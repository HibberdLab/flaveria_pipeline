#!/usr/bin/env ruby

#
# script to run sga
# based on the c.elegans example in sga/src/examples 
# 

# MOL=55
# Overlap parameter used for the final assembly. 
# OL=70
#
# Turn off collapsing bubbles around indels
# MAX_GAP_DIFF=0
#
# Parameter for the small repeat resolution algorithm
# R=10
#
# The minimum length of contigs to include in a scaffold
# MIN_LENGTH=200
#
# $SGA_BIN preprocess --pe-mode 1 -o SRR065390.fastq $IN1 $IN2
# 
# $SGA_BIN index -a ropebwt -t $CPU reads.ec.k$CK.fastq
# $SGA_BIN filter -x $COV_FILTER -t $CPU --homopolymer-check --low-complexity-check reads.ec.k$CK.fastq
# $SGA_BIN fm-merge -m $MOL -t $CPU -o merged.k$CK.fa reads.ec.k$CK.filter.pass.fa
# $SGA_BIN index -d 1000000 -t $CPU merged.k$CK.fa
# $SGA_BIN rmdup -t $CPU merged.k$CK.fa
# $SGA_BIN overlap -m $MOL -t $CPU merged.k$CK.rmdup.fa
# $SGA_BIN assemble -m $OL -g $MAX_GAP_DIFF -r $R -o assemble.m$OL merged.k$CK.rmdup.asqg.gz
# 
# CTGS=assemble.m$OL-contigs.fa
# GRAPH=assemble.m$OL-graph.asqg.gz
# 
# ~/work/devel/sga/src/bin/sga-align --name celegans.pe $CTGS $IN1 $IN2
# ~/work/devel/sga/src/bin/sga-bam2de.pl -n $MIN_PAIRS --prefix libPE celegans.pe.bam
# ~/work/devel/sga/src/bin/sga-astat.py -m $MIN_LENGTH celegans.pe.refsort.bam > libPE.astat
# 
# $SGA_BIN scaffold -m $MIN_LENGTH --pe libPE.de -a libPE.astat -o scaffolds.n$MIN_PAIRS.scaf $CTGS
# $SGA_BIN scaffold2fasta -m $MIN_LENGTH -a $GRAPH -o scaffolds.n$MIN_PAIRS.fa -d $SCAFFOLD_TOLERANCE --use-overlap --write-unplaced scaffolds.n$MIN_PAIRS.scaf

require 'trollop'

opts = Trollop::options do
  version "v0.0.2"
  opt :left, "Left input fastq file ", :required => true, :type => String
  opt :right, "Right input fastq file", :required => true, :type => String
  opt :output, "Output name", :required => true, :type => String
  # opt :interleaved, "Fastq file of interleaved left and right reads", :type => String
  opt :cores, "Number of cores", :default => 2, :type => :int
  opt :coverage_filter, "The minimum k-mer coverage for the filter step.", :default => 2, :type => :int
  opt :merge_overlap, "Overlap parameter used for FM-merge", :default=>55, :type => :int
  opt :overlap, "Overlap parameter used for the final assembly", :default => 70, :type => :int
  opt :max_gap_diff, "Turn off collapsing bubbles around indels", :default => 0, :type => :int
  opt :repeat_resolution, "Parameter for the small repeat resolution algorithm", :default => 10, :type => :int
  opt :min_pairs, "The number of pairs required to link two contigs into a scaffold", :default => 3, :type => :int
  opt :min_length, "The minimum length of contigs to include in a scaffold", :default => 200, :type => :int
  opt :scaffold_tolerance, "Distance estimate tolerance when resolving scaffold sequences", :default =>1, :type => :int
  opt :sga, "Path to sga", :default => "sga", :type => String
  opt :scaffold, "Attempt to scaffold"
  opt :test, "Don't actually run anything"
  opt :verbose, "Turn on verbose output"
end

Trollop::die :left, "must exist" if !File.exist?(opts[:left]) if opts[:left]
Trollop::die :right, "must exist" if !File.exist?(opts[:right]) if opts[:right]

path = File.dirname(opts.output)
paired = File.basename(opts.output)

puts "path = #{path}"
puts "paired = #{paired}"

sga_align_path = "/home/cmb211/apps/sga/src/bin/sga-align"  
sga_bam2de_path = "/home/cmb211/apps/sga/src/bin/sga-bam2de.pl"
sga_astat_path = "/home/cmb211/apps/sga/src/bin/sga-astat.py" 

def banner string
  w=60
  puts "#"*w
  n = ((w - string.length)*0.5).to_i
  puts " "*n << string << " "*n
  puts "#"*w
end

if !File.exists?("#{path}/#{paired}.fastq")
  banner "PREPROCESS"
  preprocess = "#{opts.sga} preprocess --pe-mode 1 -o #{path}/#{paired}.fastq #{opts.left} #{opts.right}"
  puts preprocess
  `#{preprocess}` if !opts.test
end

if !File.exists?("#{path}/#{paired}.bwt")
  banner "INDEX1"
  index = "#{opts.sga} index -a ropebwt -t #{opts.cores} -p #{path}/#{paired} #{path}/#{paired}.fastq"
  puts index
  `#{index}` if !opts.test
  if !File.exists?("#{path}/#{paired}.bwt")
    abort "index put stuff in the wrong place again!"
  end
end

if !File.exists?("#{path}/#{paired}.filter.pass.fa")
  banner "FILTER"
  filter = "#{opts.sga} filter -x #{opts.coverage_filter} -t #{opts.cores} --homopolymer-check --low-complexity-check "
  filter << " -o #{path}/#{paired}.filter.pass.fa -p #{path}/#{paired} #{path}/#{paired}.fastq"
  puts filter
  `#{filter}` if !opts.test
end

if !File.exists?("#{path}/#{paired}.merged.fa")
  banner "FM MERGE"
  fm_merge = "#{opts.sga} fm-merge -m #{opts.merge_overlap} -t #{opts.cores} -p #{path}/#{paired}.filter.pass -o #{path}/#{paired}.merged.fa #{path}/#{paired}.filter.pass.fa"
  puts fm_merge
  `#{fm_merge}` if !opts.test
end

if !File.exists?("#{path}/#{paired}.merged.bwt")
  banner "INDEX2"
  index = "#{opts.sga} index -d 1000000 -t #{opts.cores} -p #{path}/#{paired}.merged #{path}/#{paired}.merged.fa"
  puts index
  `#{index}` if !opts.test
end

if !File.exists?("#{path}/#{paired}.merged.rmdup.fa") # Fb_sga.merged.rmdup.fa
  banner "REMOVE DUPLICATES"
  rm_dup = "#{opts.sga} rmdup -t #{opts.cores} -o #{path}/#{paired}.merged.rmdup.fa #{path}/#{paired}.merged.fa"
  puts rm_dup
  `#{rm_dup}` if !opts.test
end

if !File.exists?("#{path}/#{paired}.merged.rmdup.asqg.gz")
  banner "OVERLAP"
  overlap = "#{opts.sga} overlap -m #{opts.merge_overlap} -t #{opts.cores} #{path}/#{paired}.merged.rmdup.fa"
  puts overlap
  `#{overlap}` if !opts.test
end

if !File.exists?("#{path}/#{paired}.assemble-contigs.fa")
  banner "ASSEMBLE"
  assemble = "#{opts.sga} assemble -m #{opts.overlap} -g #{opts.max_gap_diff} -r #{opts.repeat_resolution} -o #{path}/#{paired}.assemble #{path}/#{paired}.merged.rmdup.asqg.gz"
  puts assemble
  `#{assemble}` if !opts.test
end

banner "FINISHED ASSEMBLY"

#
# this section doesn't seem to do anything. the scaffolder just produces a file of contig names
#

if opts.scaffold

  if !File.exists?("Fb_sga.assemble-contigs.fa.bwt")
    banner "SGA ALIGN"
    sga_align = "#{sga_align_path} --name #{paired}.pe #{paired}.assemble-contigs.fa #{opts.left} #{opts.right}"
    puts sga_align
    `#{sga_align}` if !opts.test
  end

  if !File.exists?("#{paired}.PE.de")
    banner "BAM => DE"
    bam2de = "#{sga_bam2de_path} -n #{opts.min_pairs} -m #{opts.min_length} --prefix #{paired}.PE #{paired}.pe.bam"
    puts bam2de
    `#{bam2de}` if !opts.test
  end

  if !File.exists?("#{paired}.libPE.astat")
    banner "SGA-ASTAT.PY"
    sga_astat = "#{sga_astat_path} -m #{opts.min_length} #{paired}.pe.refsort.bam > #{paired}.libPE.astat"
    puts sga_astat
    `#{sga_astat}` if !opts.test
  end

  if !File.exists?("#{paired}.scaffolds")
    banner "SCAFFOLD"
    scaffold = "#{opts.sga} scaffold -m #{opts.min_length} --pe #{paired}.PE.de -a #{paired}.libPE.astat -o #{paired}.scaffolds #{paired}.assemble-contigs.fa"
    puts scaffold
    `#{scaffold}` if !opts.test
  end

  # $SGA_BIN scaffold2fasta -m $MIN_LENGTH -a $GRAPH -o scaffolds.n$MIN_PAIRS.fa -d $SCAFFOLD_TOLERANCE --use-overlap --write-unplaced scaffolds.n$MIN_PAIRS.scaf

  # GRAPH=assemble.m$OL-graph.asqg.gz
  banner "scaffold -> fasta"
  scaffold2fasta = "#{opts.sga} scaffold2fasta -m #{opts.min_length} -a #{paired}.merged.rmdup.asqg.gz "
  scaffold2fasta << "-o #{paired}_scaffolds.fa -d #{opts.scaffold_tolerance} --use-overlap --write-unplaced #{paired}.scaffolds"
  puts scaffold2fasta
  `#{scaffold2fasta}` if !opts.test

  ## possibly add a `sga gapfill` run here if you get a lot of NNNs in the scaffolds

  # banner "GAPFILL"
  # gapfill = "#{opts.sga} gapfill -p #{paired} -s 91 -e 41 -x 1 -t #{opts.cores} -d 64"
  # puts gapfill
  # `#{gapfill}` if !opts.test

  banner "FINISHED SCAFFOLDING"

end