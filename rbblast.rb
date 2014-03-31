#!/usr/bin/env ruby

#
# rbblast
# 

require 'rubygems'
require 'trollop'
require 'bio'

class Hit
  # Fields: query id, subject id, % identity, alignment length, mismatches,
  # gap opens, q. start, q. end, s. start, s. end, evalue, bit score
  attr_accessor :query, :target, :id, :alnlen, :mismatches, :gaps, :qstart, :qend, :tstart, :tend, :evalue, :bitscore
  def initialize(list)
    @query      = list[0].split(/[\|\ ]/).first
    @target     = list[1].split(/[\|\ ]/).first
    @id         = list[2]
    @alnlen     = list[3].to_i
    @mismatches = list[4].to_i
    @gaps       = list[5].to_i
    @qstart     = list[6].to_i
    @qend       = list[7].to_i
    @tstart     = list[8].to_i
    @tend       = list[9].to_i
    @evalue     = list[10].to_f
    @bitscore   = list[11].to_f
  end
  def to_s
    return "#{@query}\t#{@target}\t#{@id}\t#{@alnlen}\t#{@evalue}\t#{@bitscore}\ttype:#{@type}\tpath:#{@path}"
  end
end

opts = Trollop::options do
  version "v0.0.4"
  banner <<-EOS
  rbblast - reciprocal best blast

  author: Chris Boursnell (cmb211@cam.ac.uk)
  ideas and help: Richard X Smith (rds45@cam.ac.uk)
  ideas: Steve Kelly

  EOS
  opt :query, "Query in nucleotide fasta format (required)", :required => true, :type => String
  opt :target, "Target fasta (required)", :required => true, :type => String
  opt :nucl, "Target fasta file is in nucleotide format"
  opt :prot, "Target fasta file is in protein format"
  opt :output, "Output directory", :required => true, :type => String
  opt :cores, "Number of cores", :default => 2, :type => :int
  opt :nofasta, "Don't output a fasta file"
  opt :verbose, "Turn on verbose output"
  opt :test, "Don't actually do anything"
end

Trollop::die :query, "must exist" if !File.exist?(opts[:query]) if opts[:query]
Trollop::die :target, "must exist" if !File.exist?(opts[:target]) if opts[:target]
Trollop::die "One or the other Timmy, not both" if opts.nucl and opts.prot
Trollop::die "Either nucleotide or protein format must be specified for the target fasta file" if !opts.nucl and !opts.prot

## build the database of the target

query_name = File.basename(opts.query).split(".").first    # so /data/input.fasta => input
target_name = File.basename(opts.target).split(".").first

makedb1_cmd = "makeblastdb -in #{opts.target} " 
makedb1_cmd << " -dbtype nucl" if opts.nucl
makedb1_cmd << " -dbtype prot" if opts.prot
makedb1_cmd << " -title #{target_name} -out #{target_name}"

makedb2_cmd = "makeblastdb -in #{opts.query} "
makedb2_cmd << " -dbtype nucl" 
makedb2_cmd << " -title #{query_name} -out #{query_name}"

db1 = "#{target_name}.nsq" if opts.nucl
db1 = "#{target_name}.psq" if opts.prot
db2 = "#{query_name}.nsq"
if !File.exists?("#{db1}")
  puts makedb1_cmd if opts.verbose
  `#{makedb1_cmd}` if !opts.test
end
if !File.exists?("#{db2}")
  puts makedb2_cmd if opts.verbose
  `#{makedb2_cmd}` if !opts.test
end

output1 = "#{query_name}_into_#{target_name}.1.blast"
output2 = "#{target_name}_into_#{query_name}.2.blast"

if !Dir.exists?("#{opts.output}")
  mkdir = "mkdir #{opts.output}"
  `#{mkdir}` # if !opts.test
end

if opts.nucl
  ## nucl -> nucl 
  ## blastn
  ## time blastn -query Fr_200.fasta -task blastn -db Fp -out Fr_into_Fp.blast -evalue 1e-5 -outfmt 6 -num_descriptions 50 -num_alignments 50 -num_threads 23
  cmd1 = "blastn -query #{opts.query} -task blastn -db #{target_name} -out #{opts.output}/#{output1} -evalue 1e-5 -outfmt 6 -num_descriptions 50 -num_alignments 50 -num_threads #{opts.cores}"
  cmd2 = "blastn -query #{opts.target} -task blastn -db #{query_name} -out #{opts.output}/#{output2} -evalue 1e-5 -outfmt 6 -num_descriptions 50 -num_alignments 50 -num_threads #{opts.cores}"
  if !File.exists?("#{opts.output}/#{output1}")
    puts cmd1 if opts.verbose if !File.exists?("#{opts.output}/reciprocal_hits.txt")
    `#{cmd1}` if !opts.test if !File.exists?("#{opts.output}/reciprocal_hits.txt")
  end
  if !File.exists?("#{opts.output}/#{output2}")
    puts cmd2 if opts.verbose if !File.exists?("#{opts.output}/reciprocal_hits.txt")
    `#{cmd2}` if !opts.test if !File.exists?("#{opts.output}/reciprocal_hits.txt")
  end
elsif opts.prot
  ## nucl -> protein
  ## blastx query into target
  ## tblastn target into query
  cmd1 = "blastx  -query #{opts.query} -db #{target_name} -out #{opts.output}/#{output1} -evalue 1e-5 -outfmt 6 -num_descriptions 50 -num_alignments 50 -num_threads #{opts.cores}"
  cmd2 = "tblastn -query #{opts.target} -db #{query_name} -out #{opts.output}/#{output2} -evalue 1e-5 -outfmt 6 -num_descriptions 50 -num_alignments 50 -num_threads #{opts.cores}"
  if !File.exists?("#{opts.output}/#{output1}")
    puts cmd1 if opts.verbose if (!File.exists?("#{opts.output}/#{output1}") and !File.exists?("#{opts.output}/reciprocals_hits.txt"))
    `#{cmd1}` if !opts.test if (!File.exists?("#{opts.output}/#{output1}") and !File.exists?("#{opts.output}/reciprocal_hits.txt"))
  end
  if !File.exists?("#{opts.output}/#{output2}")
    puts cmd2 if opts.verbose if (!File.exists?("#{opts.output}/#{output2}") and !File.exists?("#{opts.output}/reciprocal_hits.txt"))
    `#{cmd2}` if !opts.test if (!File.exists?("#{opts.output}/#{output2}") and !File.exists?("#{opts.output}/reciprocal_hits.txt"))
  end
end

## hopefully this next bit should be fairly similar to rbusearch...

# output = "reciprocal_hits_#{query_name}_into_#{target_name}.txt"

# hashes to store blast output
query_results = Hash.new
target_results = Hash.new
count=0
longest=0
if !File.exists?("#{opts.output}/reciprocal_hits.txt") and !opts.test
  puts "Opening #{output1}" if opts.verbose
  File.open("#{opts.output}/#{output1}").each_line do |line|
    cols = line.chomp.split("\t")
    print "." if opts.verbose and count % 100_000 == 0
    count += 1
    hit = Hit.new(cols)

    # longest = h[:length] if h[:length] > longest
    query_results[hit.query] = [] if !query_results.has_key?(hit.query)
    query_results[hit.query] << hit
  end

  puts "\nOpening #{output2}" if opts.verbose
  File.open("#{opts.output}/#{output2}").each_line do |line|
    cols = line.chomp.split("\t")
    print "." if opts.verbose and count % 100_000 == 0
    count += 1
    hit = Hit.new(cols)

    target_results[hit.query] = [] if !target_results.has_key?(hit.query)
    target_results[hit.query] << hit
  end
  puts "Done" if opts.verbose
  
  evalues = [] # e-value and length of reciprocal hits
  missed_evalues = [] # e-value and length of non-reciprocal hits 

  reciprocals = Hash.new
  missed = Hash.new # hash of best hits that weren't reciprocal

  query_results.each_pair do |query_id, list_of_hits|
    best_hit_1 = list_of_hits[0] # as the results are sorted the best one is at the top
    if target_results.has_key?(best_hit_1.target)
      list_of_hits_2 = target_results[best_hit_1.target]
      best_hit_2 = list_of_hits_2[0]
      e = best_hit_2.evalue.to_f
      e = 1e-200 if e==0
      e = -Math.log10(e)
      if best_hit_2.target == query_id # is a reciprocal hit
        reciprocals[best_hit_1.query] = best_hit_1
        longest = best_hit_1.alnlen  if best_hit_1.alnlen > longest
        evalues << {:e => e, :length => best_hit_2.alnlen} #if e<200
      else
        missed[best_hit_1.query] = best_hit_1
        missed_evalues << {:e => e, :length => best_hit_2.alnlen} if e<200
      end
    end
  end

  e_data = "evalue\tlen\n"
  evalues.each do |h|
    e_data << "#{h[:e]}\t#{h[:length]}\n"
  end

  File.open("#{opts.output}/e_data.txt", "w") { |io|  io.write(e_data)}

  length_hash = Hash.new
  
  evalues.each do |h|
    length_hash[h[:length]] = [] if !length_hash.has_key?(h[:length])
    length_hash[h[:length]] << h
  end

  # puts "Longest #{longest}"
  puts "Line fitting" if opts.verbose
  fitting = Hash.new
  (10..longest).each do |centre|
    e = 0
    count = 0
    s = centre*0.1
    s = s.to_i
    s = 5 if s < 5
    (-s..s).each do |side|
      if length_hash.has_key?(centre+side)
        length_hash[centre+side].each do |point|
          e += point[:e]
          count += 1
        end
      end
    end
    if count>0
      mean = e/count
      fitting[centre] = mean
    end
  end

  fitting_data=""
  fitting.keys.sort.each do |centre|
    fitting_data << "#{centre}\t#{fitting[centre]}\n"
  end
  File.open("#{opts.output}/fitting.txt", "w") {|io| io.write(fitting_data)}
  
  output = ""
  output_2 = ""
  reciprocals.each_pair do |id, hit|
    output << "#{hit.query}\t#{hit.target}\t#{hit.id}\t#{hit.alnlen}\t#{hit.mismatches}\t#{hit.gaps}\t#{hit.qstart}\t#{hit.qend}\t#{hit.tstart}\t#{hit.tend}\t#{hit.evalue}\t#{hit.bitscore}\t1\n"
  end

  missed.each_pair do |id, hit|
    l = hit.alnlen.to_i
    e = hit.evalue
    e = 1e-200 if e==0
    e = -Math.log10(e)

    if fitting.has_key?(l)
      if e >= fitting[l]
        output << "#{hit.query}\t#{hit.target}\t#{hit.id}\t#{hit.alnlen}\t#{hit.mismatches}\t#{hit.gaps}\t#{hit.qstart}\t#{hit.qend}\t#{hit.tstart}\t#{hit.tend}\t#{hit.evalue}\t#{hit.bitscore}\t2\n"
      end
      if !reciprocals.has_key?(id)
        reciprocals[id] = hit # adding so that these can be added to the fasta file
      end
    else
      output_2 << "#{hit.query}\t#{hit.target}\t#{hit.id}\t#{hit.alnlen}\t#{hit.mismatches}\t#{hit.gaps}\t#{hit.qstart}\t#{hit.qend}\t#{hit.tstart}\t#{hit.tend}\t#{hit.evalue}\t#{hit.bitscore}\t3\n"
    end
  end

  File.open("#{opts.output}/reciprocal_hits.txt", "w") { |io| io.write(output) }
  File.open("#{opts.output}/not_reciprocal_hits.txt", "w") { |io| io.write(output_2)}

  if !opts.nofasta
    # create new fasta file with updated contig names
    output=""
    fasta = Bio::FastaFormat.open(opts.query)
    fasta.each do |entry|
      contig = entry.definition
      if reciprocals.has_key?(contig)
        contig = "#{contig}_#{reciprocals[contig].target}"
      end
      output << ">#{contig}\n"
      output << "#{entry.seq}\n"
    end
    File.open("#{opts.output}/#{query_name}_annotated.fasta", "w") { |io|  io.write(output)}
    puts "Annotated fasta file written" if opts.verbose
  end
end