#!/usr/bin/env ruby

#
# deinterleave fastq reads
#
# checks to make sure that nothing has gone wrong
# while running khmer
#

require 'trollop'

opts = Trollop::options do
  version "v0.0.1a"
  opt :fastq, "Fastq file to deinterleave", :required => true, :type => String
  opt :output, "Output prefix", :required => true, :type => String
  opt :verbose, "Be verbose"
end
Trollop::die :fastq, "must exist" if !File.exist?(opts[:fastq]) if opts[:fastq]

fastq = File.open("#{opts.fastq}", "r")

output_left = File.open("#{opts.output}.left.fastq", "w")
output_right= File.open("#{opts.output}.right.fastq", "w")

name1 = fastq.readline rescue nil
seq1 = fastq.readline rescue nil
plus1 = fastq.readline rescue nil
qual1 = fastq.readline rescue nil
name2 = fastq.readline rescue nil
seq2 = fastq.readline rescue nil
plus2 = fastq.readline rescue nil
qual2 = fastq.readline rescue nil

while name1 != nil and name2 != nil
  if name1.split(":").first == name2.split(":").first
    if seq1.length!=qual1.length
      abort "error: seq1.length != qual1.length at #{name1}"
    end
    if seq2.length!=qual2.length
      abort "error: seq2.length != qual2.length at #{name2}"
    end
    output_left.write(name1)
    output_left.write(seq1)
    output_left.write(plus1)
    output_left.write(qual1)
    output_right.write(name2)
    output_right.write(seq2)
    output_right.write(plus2)
    output_right.write(qual2)
  else
    abort "Houston we have a problem"
  end

  name1 = fastq.readline rescue nil
  seq1 = fastq.readline rescue nil
  plus1 = fastq.readline rescue nil
  qual1 = fastq.readline rescue nil
  name2 = fastq.readline rescue nil
  seq2 = fastq.readline rescue nil
  plus2 = fastq.readline rescue nil
  qual2 = fastq.readline rescue nil

end
