#!/usr/bin/env ruby

#
# the expression output part of the pipeline rakefile as a ruby script
#
# Chris Boursnell (cmb211@cam.ac.uk) 2014-02-14
#

require 'csv'
require 'rubygems'
require 'trollop'

class Express
  attr_accessor :bundle_id,:target_id,:length,:eff_length,:tot_counts,:uniq_counts,:est_counts,:eff_counts,:ambig_distr_alpha,:ambig_distr_beta,:fpkm,:fpkm_conf_low,:fpkm_conf_high,:solvable,:tpm

  def initialize(line)
    cols = line.chomp.split("\t")
    @bundle_id = cols[0].to_f
    @target_id = cols[1]
    @length = cols[2].to_f
    @eff_length = cols[3].to_f
    @tot_counts = cols[4].to_f
    @uniq_counts = cols[5].to_f
    @est_counts = cols[6].to_f
    @eff_counts = cols[7].to_f
    @ambig_distr_alpha = cols[8].to_f
    @ambig_distr_beta = cols[9].to_f
    @fpkm = cols[10].to_f
    @fpkm_conf_low = cols[11].to_f
    @fpkm_conf_high = cols[12].to_f
    @solvable = cols[13]
    @tpm = cols[14].to_f
  end
end

opts = Trollop::options do
  version "v0.0.1a"
  opt :reads, "CSV file containing list of trimmed read fastq files (filepath, rep, section, type) type is 1=left, 2=right, 3=unpaired", :required => true, :type => String
  opt :fasta, "Fasta file of contigs to map the reads against", :required => true, :type => String
  opt :path, "Output path", :required => true, :type => String
  opt :prefix, "Prefix", :required => true, :type => String
  opt :threads, "Threads", :default => 4, :type => :int
  opt :output, "Where to save the output summary file", :required => true, :type => String
  opt :annotations, "Reciprocal hits output file from rbusearch", :type => String
  opt :all, "Align all reads"
  opt :k, "Align k reads", :type => :int
  opt :node, "What node is this?", :type => String
  opt :test, "Test"
  opt :verbose, "Be verbose"
end

Trollop::die :reads,   "must exist" if !File.exist?(opts[:reads])   if opts[:reads]
Trollop::die "Must select either -a or -k " if !opts.all and !opts.k

node=opts.node

left = Hash.new
right = Hash.new
single = Hash.new
CSV.foreach("#{opts.reads}") do |row|
  section = "#{row[1]}-#{row[2]}".to_sym
  if row[3]=="1"
    left[section]=row[0]
  elsif row[3]=="2"
    right[section]=row[0]
  elsif row[3]=="3"
    single[section]=[] if !single.has_key?(section)
    single[section] << row[0]
  else
    abort "Unexpected value in 4th column"
  end
end

bowtie2 = "/home/cmb211/apps/bowtie2/bowtie2"

# MAKE A BOWTIE INDEX
if !File.exists?("#{opts.path}/#{opts.prefix}.index.4.bt2")
  index_cmd = "bowtie2-build #{opts.fasta} #{opts.path}/#{opts.prefix}.index"
  puts index_cmd if opts.verbose
  `#{index_cmd}` if !opts.test
end

left.keys.each do |section|

  if !File.exists?("#{opts.path}/#{opts.prefix}#{section}.sam") # if the sam file doesn't exist run bowtie
    bowtie_cmd = "#{bowtie2} -t --very-sensitive -p #{opts.threads} -x #{opts.path}/#{opts.prefix}.index "
    bowtie_cmd << " -a " if opts.all
    bowtie_cmd << " -k #{opts.k} " if opts.k
    bowtie_cmd << " -1 #{left[section]} "
    bowtie_cmd << " -2 #{right[section]} "
    bowtie_cmd << " -U #{single[section].join(",")} "
    # bowtie_cmd << " -S #{opts.path}/#{opts.prefix}#{section}.sam"
    if node=="node9"
      bowtie_cmd << " -S /disk2/tmp/cmb211/#{opts.prefix}#{section}.sam"
    elsif node=="node8"
      bowtie_cmd << " -S /tmp/cmb211/#{opts.prefix}#{section}.sam" 
    else
      bowtie_cmd << " -S #{opts.path}/#{opts.prefix}#{section}.sam"
    end
    puts bowtie_cmd if opts.verbose
    `#{bowtie_cmd}` if !opts.test
    if node=="node9"
      mv_cmd = "mv /disk2/tmp/cmb211/#{opts.prefix}#{section}.sam #{opts.path}/#{opts.prefix}#{section}.sam" 
    elsif node=="node8"
      mv_cmd = "mv /tmp/cmb211/#{opts.prefix}#{section}.sam #{opts.path}/#{opts.prefix}#{section}.sam" 
    else
      mv_cmd = "echo \"do not have to copy anything\""
    end
    puts mv_cmd if opts.verbose
    `#{mv_cmd}` if !opts.test
  end

  if !File.exists?("#{opts.path}/express_#{opts.prefix}#{section}/results.xprs") # if the eXpress output doesn't exist run eXpress
    express_cmd = "express --output-align-prob "
    express_cmd << " -o #{opts.path}/express_#{opts.prefix}#{section} "
    express_cmd << " --no-update-check "   
    express_cmd << " -B 2 " 
    express_cmd << " #{opts.fasta} " # fasta file to align reads to
    express_cmd << " #{opts.path}/#{opts.prefix}#{section}.sam"
    puts express_cmd if opts.verbose
    `#{express_cmd}` if !opts.test
  end

  # if opts.remove_sam # NOOOOO don't delete the sam files yet
  #   rm_cmd = "rm #{opts.path}/#{opts.prefix}#{section}.sam"
  #   `#{rm_cmd}` if !opts.test
  # end
end

puts "Filling hash..."

expression_hash = Hash.new
left.keys.each do |section|
  puts "opening #{section}"
  if File.exists?("#{opts.path}/express_#{opts.prefix}#{section}/results.xprs")
    File.open("#{opts.path}/express_#{opts.prefix}#{section}/results.xprs", "r").each_line do |line|
      e = Express.new(line)
      if e.target_id != "target_id"
        expression_hash[e.target_id] = Hash.new if !expression_hash.has_key?(e.target_id)
        expression_hash[e.target_id][section] = e
      end
    end
  end
end

annotation_hash = Hash.new
if opts.annotations
  # if a 'reciprocal hits' annotation file has been provided then use the names in there to add
  # to the output
  File.open("#{opts.annotations}", "r").each_line do |line|
    cols = line.chomp.split("\t")
    contig = cols[0]
    agi = cols[1]
    annotation_hash[contig] = agi
  end
end

puts "writing output"
output = File.open("#{opts.path}/#{opts.output}", "w")
output.write("contig\tlength\t")
left.keys.each do |section|
  output.write("#{section}\t")
end
output.write("t1\tt2\tt3\tt3\tt4\tt5\tt6")
output.write("\n")

expression_hash.each_pair do |contig, hash|
  if opts.annotations && annotation_hash.has_key?(contig)
    output.write("#{contig}_#{annotation_hash[contig]}\t")
  else
    output.write("#{contig}\t")
  end
  output.write("#{hash[left.keys[0]].length}\t")
  # totals = [] # TODO fill 0-5 with zeros
  totals = [0,0,0,0,0,0]
  hash.each_pair do |section, data|
    (rep, leaf_section) = section.to_s.split("-")
    totals[leaf_section.to_i-1] += data.tpm
    output.write("#{data.eff_counts}\t") # TPM
  end
  totals.each do |i|
    if i!=nil
      output.write("#{i/3}\t") 
    else
      output.write("0.0\t")
    end
  end
  output.write("\n")
end
puts "Done"
