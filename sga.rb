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
  opt :sga, "Path to sga", :default => "sga", :type => String
  # opt :scaffold, "Attempt to scaffold"
  opt :test, "Don't actually run anything"
  opt :verbose, "Turn on verbose output"

  # filter
  opt :filter_coverage, "The minimum k-mer coverage for the filter step.", :default => 3, :type => :int
  opt :filter_kmer, "The length of the kmer to use.", :default => 27, :type => :int

  opt :filter_no_dup, "turn off duplicate removal"
  opt :filter_substring_only, "when removing duplicates, only remove substring sequences, not full-length matches"
  opt :filter_kmer, "turn off the kmer check"
  opt :filter_homopolymer, "check reads for hompolymer run length sequencing errors"
  opt :filter_low_comp, "filter out low complexity reads"

  # fm-merge
  opt :merge_overlap, "Overlap parameter used for FM-merge", :default=>45, :type => :int
  
  # rmdup
  opt :rmdup_error_rate, "the maximum error rate allowed to consider two sequences identical (default: exact matches required)", :default => 0, :type => :float # -e

  # overlap
  opt :error_rate, "the maximum error rate allowed to consider two sequences aligned ", :default => 0, :type => :float # -e
  opt :min_overlap, "minimum overlap required between two reads", :default => 45, :type => :int # -m
  opt :exhaustive, "output all overlaps, including transitive edges" # yes/no # -x
  opt :seed_length, "force the seed length to be LEN. By default, the seed length in the overlap step is calculated to guarantee all overlaps with --error-rate differences are found.  This option removes the guarantee but will be (much) faster. As SGA can tolerate some missing edges, this option may be preferable for some data sets.", :type => :int # -l

  # assemble
  opt :assembly_overlap, "Overlap parameter used for the final assembly", :default => 70, :type => :int 
  opt :transitive_reduction, "remove transitive edges from the graph. Off by default."
  opt :max_edges, "limit each vertex to a maximum of N edges. For highly repetitive regions this helps save memory by culling excessive edges around unresolvable repeats", :default => 128, :type =>:int
  opt :bubble_remove, "perform N bubble removal steps", :default => 3, :type => :int
  opt :max_diff, "only remove variation if the divergence between sequences is less than F", :default => 0.05, :type => :float
  opt :max_gap_diff, "Only remove variation if the divergence between sequences when only counting indels is less than this", :default => 0.01, :type => :float
  opt :max_indel, "do not remove variation that is an indel of length greater than this", :default => 20, :type => :int
  opt :cut_terminal, "cut off terminal branches in N rounds", :default=> 10, :type => :int
  opt :min_branch_length, "remove terminal branches only if they are less than LEN bases in length", :default => 150, :type => :int
  opt :repeat_resolution, "Parameter for the small repeat resolution algorithm", :default => 10, :type => :int
  

  # opt :min_pairs, "The number of pairs required to link two contigs into a scaffold", :default => 3, :type => :int
  # opt :min_length, "The minimum length of contigs to include in a scaffold", :default => 200, :type => :int
  # opt :scaffold_tolerance, "Distance estimate tolerance when resolving scaffold sequences", :default =>1, :type => :int
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
  filter = "#{opts.sga} filter "
  filter << "-t #{opts.cores} "
  filter << "-p #{path}/#{paired} "
  filter << "-o #{path}/#{paired}.filter.pass.fa "

  filter << "-x #{opts.filter_coverage} " if opts.filter_coverage
  filter << "-l #{opts.filter_kmer} " if opts.filter_kmer

  filter << "--no-duplicate-check " if opts.filter_no_dup 
  filter << "--substring-only " if opts.filter_substring_only 
  filter << "--no-kmer-check " if opts.filter_kmer 
  filter << "--homopolymer-check " if opts.filter_homopolymer 
  filter << "--low-complexity-check " if opts.filter_low_comp

  filter << "#{path}/#{paired}.fastq"
  puts filter
  `#{filter}` if !opts.test
end

# opt :filter_coverage, "The minimum k-mer coverage for the filter step.", :default => 3, :type => :int
# opt :filter_kmer, "The length of the kmer to use.", :default => 27, :type => :int

# opt :filter_no_dup, "turn off duplicate removal"
# opt :filter_substring_only, "when removing duplicates, only remove substring sequences, not full-length matches"
# opt :filter_kmer, "turn off the kmer check"
# opt :filter_homopolymer, "check reads for hompolymer run length sequencing errors"
# opt :filter_low_comp, "filter out low complexity reads"

if !File.exists?("#{path}/#{paired}.merged.fa")
  banner "FM MERGE"
  fm_merge = "#{opts.sga} fm-merge "
  fm_merge << "-m #{opts.merge_overlap} " if opts.merge_overlap
  fm_merge << "-t #{opts.cores} "
  fm_merge << "-p #{path}/#{paired}.filter.pass "
  fm_merge << "-o #{path}/#{paired}.merged.fa "
  fm_merge << "#{path}/#{paired}.filter.pass.fa"
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
  rm_dup = "#{opts.sga} rmdup "
  rm_dup << "-t #{opts.cores} "
  rm_dup << "-e #{opts.rmdup_error_rate} " if opts.rmdup_error_rate
  rm_dup << "-o #{path}/#{paired}.merged.rmdup.fa "
  rm_dup << " #{path}/#{paired}.merged.fa "
  puts rm_dup
  `#{rm_dup}` if !opts.test
end

if !File.exists?("#{path}/#{paired}.merged.rmdup.asqg.gz")
  banner "OVERLAP"
  overlap = "#{opts.sga} overlap "
  overlap << "-t #{opts.cores} "
  overlap << "-e #{opts.error_rate}" if opts.error_rate
  overlap << "-m #{opts.min_overlap} " if opts.min_overlap
  overlap << "-x " if opts.exhaustive
  overlap << "-l #{opts.seed_length}" if opts.seed_length
  overlap << " #{path}/#{paired}.merged.rmdup.fa"
  puts overlap
  `#{overlap}` if !opts.test
end

if !File.exists?("#{path}/#{paired}.assemble-contigs.fa")
  banner "ASSEMBLE"
  assemble = "#{opts.sga} assemble "
  assemble << "-m #{opts.assembly_overlap} " if opts.assembly_overlap
  assemble << "--transitive-reduction " if opts.transitive_reduction
  assemble << "--max-edges #{opts.max_edges}" if opts.max_edges
  assemble << "-b #{opts.bubble_remove} " if opts.bubble_remove
  assemble << "-d #{opts.max_diff} " if opts.max_diff
  assemble << "-g #{opts.max_gap_diff} " if opts.max_gap_diff
  assemble << "--max-indel=#{opts.max_indel} " if opts.max_indel
  assemble << "-x #{opts.cut_terminal} " if opts.cut_terminal
  assemble << "-l #{opts.min_branch_length} " if opts.min_branch_length
  assemble << "-r #{opts.repeat_resolution} " if opts.repeat_resolution
  assemble << "-o #{path}/#{paired}.assemble " # output
  assemble << "#{path}/#{paired}.merged.rmdup.asqg.gz" # input
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