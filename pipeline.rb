#!/usr/bin/env ruby

#
# this is a non-rake version of the flaveria pipeline that can be run from anywhere
#


require 'rubygems'
require 'trollop'


opts = Trollop::options do
  version "v0.0.1a"
  opt :left, "Left reads", :required => true, :type => String
  opt :right, "Right reads", :required => true, :type => String
  opt :single, "Single reads", :required => true, :type => String
  opt :path, "Path", :required => true, :type => String
  opt :prefix, "Prefix", :required => true, :type => String
  opt :threads, "Threas", :default => 4, :type => :int
  opt :test, "Test"
  opt :verbose, "Be verbose"
end

Trollop::die :left,   "must exist" if !File.exist?(opts[:left])   if opts[:left]
Trollop::die :right,  "must exist" if !File.exist?(opts[:right])  if opts[:right]
Trollop::die :single, "must exist" if !File.exist?(opts[:single]) if opts[:single]


idba = "/home/cmb211/apps/idba_tran-1.0.13/bin/idba_tran"
sga_path = "/home/cmb211/scripts/flaveria_pipeline/sga.rb"
cd_hit_est = "/home/cmb211/apps/cd-hit-v4.6.1-2012-08-27/cd-hit-est"
gapcloser = "/home/cmb211/bin/GapCloser"

soap_scafseq = "#{opts.path}/soap/#{opts.prefix}soap.scafSeq"
soap_contigs = "#{opts.path}/soap/#{opts.prefix}soap.gapfill.fasta"
sga_contigs = "#{opts.path}/sga/#{opts.prefix}sga.assemble-contigs.fa"
idba_contigs = "#{opts.path}/idba/contig.fa"

## SOAP  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

config = "max_rd_len=20000\n"
config << "[LIB]\n"
config << "avg_ins=250\n"
config << "reverse_seq=0\n"
config << "asm_flags=3\n"
# config << "rank=2\n"
config << "q1=#{opts.left}\n"
config << "q2=#{opts.right}\n"
config << "q=#{opts.single}\n"
puts config if opts.verbose
File.open("#{opts.path}/soapdt.config", "w") {|out| out.write config} if !opts.test

if !Dir.exists?("#{opts.path}/soap")
  mkdir_cmd = "mkdir #{opts.path}/soap"
  puts mkdir_cmd if opts.verbose
  `#{mkdir_cmd}`  if !opts.test
end

if !File.exists?("#{soap_scafseq}")
  soap_cmd = "SOAPdenovo-Trans-127mer all -s #{opts.path}/soapdt.config -o #{opts.path}/soap/#{opts.prefix}soap -p #{opts.threads}"
  puts soap_cmd if opts.verbose
  `#{soap_cmd}` if !opts.test
else
  puts "soap already run #{soap_scafseq}"
end

if !File.exists?("#{opts.path}/soap/#{opts.prefix}soap.gapcloser.fasta")
  # run soap gapcloser on scafseq
  gapcloser_cmd = "#{gapcloser} -a #{opts.path}/soap/#{opts.prefix}soap.scafSeq -b #{opts.path}/soapdt.config -o #{opts.path}/soap/#{opts.prefix}soap.gapcloser.fasta" 
  puts gapcloser_cmd if opts.verbose
  `#{gapcloser_cmd}` if !opts.test
  abort "something went wrong here" if File.zero?("#{opts.path}/soap/#{opts.prefix}soap.gapcloser.fasta")
else
  puts "gapcloser already run"
end

if !File.exists?("#{soap_contigs}")
  # run sga gap filler on *.scafSeq
  pre_cmd =     "sga preprocess -p 1 -o #{opts.path}/soap/#{opts.prefix}cleaned.fastq #{opts.left} #{opts.right}" 
  index_cmd =   "sga index -p #{opts.path}/soap/#{opts.prefix} -t #{opts.threads} -a ropebwt #{opts.path}/soap/#{opts.prefix}cleaned.fastq"
  gapfill_cmd = "sga gapfill -p #{opts.path}/soap/#{opts.prefix} -e 47 -x 1 -t #{opts.threads} -o #{opts.path}/soap/#{opts.prefix}soap.gapfill.fasta #{opts.path}/soap/#{opts.prefix}soap.gapcloser.fasta"

  if !File.exists?("#{opts.path}/soap/#{opts.prefix}cleaned.fastq")
    puts pre_cmd if opts.verbose
    `#{pre_cmd}` if !opts.test
  end
  puts index_cmd if opts.verbose
  `#{index_cmd}` if !opts.test
  puts gapfill_cmd if opts.verbose
  `#{gapfill_cmd}` if !opts.test
else
  puts "sga gapfill already run #{soap_contigs}"
end

## SGA # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if !Dir.exists?("#{opts.path}/sga")
  mkdir_cmd = "mkdir #{opts.path}/sga"
  puts mkdir_cmd if opts.verbose
  `#{mkdir_cmd}` if !opts.test
end

if !File.exists?("#{sga_contigs}")
  sga_cmd = "ruby #{sga_path} --verbose --left #{opts.left} --right #{opts.right} --output #{opts.path}/sga/#{opts.prefix}sga --cores #{opts.threads}"
  puts sga_cmd if opts.verbose
  `#{sga_cmd}` if !opts.test
else
  puts "sga assembly already run #{sga}"
end

## IDBA #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if !Dir.exists?("#{opts.path}/idba")
  mkdir_cmd = "mkdir #{opts.path}/idba"
  puts mkdir_cmd if opts.verbose
  `#{mkdir_cmd}` if !opts.test
end

files = []
files << opts.left
files << opts.right
files << opts.single
# prepare reads
if !File.exists?("#{opts.path}/idba/#{opts.prefix}.fx.fa") and !opts.test
  fasta = File.open("#{opts.path}/idba/#{opts.prefix}.fx.fa", "w")
  files.each do |file|
    fastq = File.open("#{file}", "r")
    header = fastq.readline
    seq = fastq.readline
    plus = fastq.readline
    qual = fastq.readline
    while header != nil
      fasta.write(">#{header}")
      fasta.write("#{seq}")
      header = fastq.readline rescue nil
      seq = fastq.readline rescue nil
      plus = fastq.readline rescue nil
      qual = fastq.readline rescue nil
    end
  end
end

if !File.exists?("#{idba_contigs}")
  idba_cmd = "#{idba} "
  idba_cmd << "-o #{opts.path}/idba "
  idba_cmd << "-r #{opts.path}/idba/#{opts.prefix}.fx.fa "
  idba_cmd << "--num_threads #{opts.threads} "           # number of threads
  idba_cmd << "--mink 21 "                          # minimum k value (<=124)
  idba_cmd << "--maxk 77 "                          # maximum k value (<=124)
  idba_cmd << "--step 8 "                           # increment of k-mer of each iteration
  idba_cmd << "--min_count 1 "                      # minimum multiplicity for filtering k-mer when building the graph
  idba_cmd << "--no_correct "                       # do not do correction
  idba_cmd << "--max_isoforms 6 "                   # maximum number of isoforms
  idba_cmd << "--similar 0.98"                      # similarity for alignment

  puts idba_cmd if opts.verbose
  `#{idba_cmd}` if !opts.test
else
  puts "idba already run #{idba}"
end

## CD-HIT  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# concatenate all the output assemblies together
cat_cmd = "cat #{soap_contigs} #{sga_contigs} #{idba_contigs} > #{opts.path}/#{opts.prefix}combined_contigs.fa"
puts cat_cmd if opts.verbose
`#{cat_cmd}` if !opts.test

# run cd-hit-est
cd_cmd = "#{cd_hit_est} -i #{opts.path}/#{opts.prefix}combined_contigs.fa -o #{opts.path}/#{opts.prefix}cd_hit.fasta -T #{opts.threads} -c 0.99 -M 5000"
puts cd_cmd if opts.verbose
`#{cd_cmd}` if !opts.test

## RBUSEARCH   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmd = "ruby rbusearch.rb --query #{opts.path}/#{opts.prefix}contigs.fasta --target #{protein_reference} --output #{opts.path}/rbusearch --cores #{threads} --prefix #{opts.prefix} --verbose"
puts cmd if opts.verbose
`#{cmd}` if !opts.test

cmd = "mv #{opts.path}/rbusearch/#{opts.prefix}annotated.fasta #{opts.path}/#{opts.prefix}annotated.fasta"
`#{cmd}` if !opts.test

## EXPRESS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



