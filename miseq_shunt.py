#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "1/29/2014"
__version__ = 1.0

from optparse import OptionParser, OptionGroup
from numpy import mean
import sys
import time
import os
import subprocess
import tempfile
import gzip
import zipfile

PATH_TO_PEAR = "pear"

def fetch_unzip(zipfile_obj, filename):
  # accepts zipfile object and filename contained within, returns
  # temporary file containing raw reads

  gzipped = tempfile.NamedTemporaryFile(delete=False)
  gzipped.write(zipfile_obj.open(filename, "r").read())
  gzipped.flush()

  ungzipped = tempfile.NamedTemporaryFile(delete=False)
  ungzipped.write(gzip.GzipFile(gzipped.name).read())
  ungzipped.flush()

  os.unlink(gzipped.name)

  return ungzipped

def pear(left_fname, right_fname):
  # accepts left and right paired end reads, returns temporary file containing
  # presumably overlapping reads merged into single reads and percentage of
  # reads unable to be merged (high % indicates a problem)

  outfile = tempfile.NamedTemporaryFile(delete=False)

  pear = subprocess.Popen([PATH_TO_PEAR,
                           "-f", left_fname, "-r", right_fname,
                           "-o", outfile.name,
                           "-y", options.mem_size, "-j", str(options.num_threads),
                           "-p", "0.01"],
                          stdout=subprocess.PIPE)

  while pear.poll() == None:
    time.sleep(0.1)

  if pear.returncode <> 0:
    raise ValueError("Pear exited with an error (%s)" % pear.returncode)

  stats = {}

  for line in pear.stdout:
    line = line.strip()

    if line.endswith("%)"):
      perc = float(line.split(" ")[-1][1:-2])

      if line.startswith("Assembled reads"):
        stats["assembled_reads"] = perc
      elif line.startswith("Discarded reads"):
        stats["discarded_reads"] = perc
      elif line.startswith("Not assembled reads"):
        stats["unassembled_reads"] = perc

  return outfile, stats

def find_pairs(zipfile_obj):
  # accepts zipfile object and returns tuples of paired-end reads

  filenames = zipfile_obj.namelist()
  used_filenames = []
  pairs = []

  for fname in filenames:
    mate = fname.split("_")

    if mate[-2] == "R1":
      mate[-2] = "R2"
    elif mate[-2] == "R2":
      mate[-2] = "R1"
    else:
      raise ValueError("Invalid filename format (%s)" % fname)

    mate = "_".join(mate)

    if mate in filenames:
      used_filenames.extend((fname, mate))
      pairs.append(tuple(sorted((fname, mate))))

  unused_filenames = set(filenames).difference(used_filenames)

  if unused_filenames:
    raise ValueError("Found %s samples without pairs" % len(unused_filenames))

  return list(set(pairs))

def fast_fastq(fp_in):                                                          
    buf = []                                                                    
                                                                                
    for line in fp_in:                                                          
        buf.append(line)                                                        
                                                                                
        if len(buf) == 4:                                                       
            yield fastq_record(buf, 0)                                          
            buf = []                                                            
                                                                                
class fastq_record():                                                           
    def __init__(self, lines, offset):                                          
        lines = [x.strip() for x in lines]                                      
                                                                                
        self.id = lines[0][1:]                                                  
        self.sequence = lines[1]                                                
        self.quals = lines[3]                                                   
        self.offset = offset                                                    
                                                                                
    def raw(self):                                                              
        return "\n".join(["@%s" % (self.id,), self.sequence, "+", self.quals, ""])

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] <miseq_fastq.zip>",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./miseq_out]",
                    default="./miseq_out",
                    help="write output files to this directory")

  parser.add_option("-p",
                    dest="num_threads",
                    type="int",
                    metavar="[1]",
                    default=1,
                    help="number of threads used during analysis")

  parser.add_option("-m",
                    dest="mem_size",
                    metavar="[1G]",
                    default="1G",
                    help="amount of memory to use during read merging")

  parser.add_option("-q",
                    dest="min_qual",
                    type="int",
                    metavar="[25]",
                    default=25,
                    help="discard reads with mean quality below this value")

  parser.add_option("--phred_offset",
                    dest="phred_offset",
                    type="int",
                    metavar="[33]",
                    default=33,
                    help="phred quality offset (default 33 for illumina 1.8+)")

  parser.add_option("--min-qual-perc",
                    dest="min_qual_perc",
                    type="int",
                    metavar="[95]",
                    default=95,
                    help="warn if fewer than this % of reads pass quality filter")

  parser.add_option("--min-merged-perc",
                    dest="min_merged_perc",
                    type="int",
                    metavar="[95]",
                    default=95,
                    help="warn if fewer than this % of reads can be merged")

  options, args = parser.parse_args(arguments)

  if len(args) <> 1:
    print "Error: Incorrect number of arguments"
    parser.print_help()
    sys.exit(1)

  if not (0 <= options.min_merged_perc <= 100):
    print "Error: --min-merged-perc must be in [0,100]"
    parser.print_help()
    sys.exit(1)

  if not (0 <= options.min_qual_perc <= 100):
    print "Error: --min-qual-perc must be in [0,100]"
    parser.print_help()
    sys.exit(1)

  options.min_merged_perc /= 100.0
  options.min_qual_perc /= 100.0

def main():
  parse_options(sys.argv[1:])

  if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)
  elif not os.path.isdir(options.output_dir):
    print "Error: Specific path exists and is not a directory"
    parser.print_help()
    sys.exit(1)

  # open zipfile
  miseq_zip = zipfile.ZipFile(args[0])

  # find pairs
  pairs = [("_".join(x.split("/")[-1].split(".")[0].split("_")[:-2]), (x, y)) for x, y in find_pairs(miseq_zip)]

  min_read_length = 1e10
  max_read_length = 0
  total_read_size = 0.0
  total_read_length = 0.0

  with open(os.path.join(options.output_dir, "seqs.fna"), "w") as fp:
    # for every pair
    read_lengths = []

    for pair_name, pair in pairs:
      sys.stderr.write("Processing %s\n" % pair_name)
  
      # extract read files
      left = fetch_unzip(miseq_zip, pair[0])
      right = fetch_unzip(miseq_zip, pair[1])
  
      # run pear to merge reads
      merged, stats = pear(left.name, right.name)
  
      if stats["assembled_reads"] < options.min_merged_perc:
        sys.stderr.write("  Warning: Only %.02f%% of reads assembled\n" % stats["assembled_reads"])
 
      # read through fastq
      total_reads = 0.0
      kept_reads = 0.0

      for seq_rec in fast_fastq(open("%s.assembled.fastq" % merged.name, "r")):
        total_reads += 1

        # check that read passes quality filter
        if mean([ord(x) - options.phred_offset for x in seq_rec.quals]) < \
           options.min_qual:
          continue

        # keep track of stats
        if len(seq_rec.sequence) < min_read_length:
          min_read_length = len(seq_rec.sequence)

        if len(seq_rec.sequence) > max_read_length:
          max_read_length = len(seq_rec.sequence)

        total_read_size += 1
        total_read_length += len(seq_rec.sequence)

        # rename reads and write to pooled output
        kept_reads += 1
        seq_rec.id = "%s_%s %s" % (pair_name, kept_reads, seq_rec.id.split(" ")[0])
        fp.write(">%s\n%s\n" % (seq_rec.id, seq_rec.sequence))

      # remove temporary files
      os.unlink(left.name)
      os.unlink(right.name)
      os.unlink(merged.name)
      os.unlink("%s.assembled.fastq" % merged.name)
      os.unlink("%s.discarded.fastq" % merged.name)
      os.unlink("%s.unassembled.forward.fastq" % merged.name)
      os.unlink("%s.unassembled.reverse.fastq" % merged.name)

      if total_reads > 0:
        if kept_reads / total_reads < options.min_qual_perc:
          sys.stderr.write("  Warning: Only %.02f%% of reads passed quality filter (mean(qv) > %s)\n" % (kept_reads * 100 / total_reads, options.min_qual))

      sys.stderr.write("  %d/%d reads kept\n" % (kept_reads, total_reads))

  sys.stderr.write("\nWriting mapping file\n")

  with open(os.path.join(options.output_dir, "mapping.txt"), "w") as fp:
    fp.write("\t".join(["#SampleID", "BarcodeSequence", "LinkerPrimerSequence"]) + "\n")

    for pair_name, _ in pairs:
      fp.write("\t".join([pair_name, "", ""]) + "\n")

  sys.stderr.write("\nSummary")
  sys.stderr.write("\n  Min read length:  %d" % min_read_length)
  sys.stderr.write("\n  Mean read length: %d" % (total_read_length / total_read_size))
  sys.stderr.write("\n  Max read length:  %d\n" % max_read_length)
  sys.stderr.write("\n  Total reads:      %d" % total_read_size)
  sys.stderr.write("\n  Total samples:    %d" % len(pairs))
  sys.stderr.write("\n  Reads/sample:     %d\n" % (total_read_size / len(pairs)))

if __name__ == "__main__":
  main()
