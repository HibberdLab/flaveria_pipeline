#!/usr/bin/env ruby

# rbusearch based on steve kelly's reciprocal best blast
#

# Input:
#   Contigs that need annotating in nucleotide fasta format
#   Reference Transcriptome in aa fasta format
# Output:
#   List of contigs and their corresponding proteins from the reference
#   New fasta file with contigs renamed

require 'rubygems'
require 'trollop'
require 'bio'
require 'csv'

# a few class/method definitions:

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

class Bio::FastaFormat
  def isNucl?
    Bio::Sequence.guess(self, 0.9, 1000) == Bio::Sequence::NA
  end

  def isProt?
    Bio::Sequence.guess(self, 0.9, 1000) == Bio::Sequence::AA
  end
end

def load_query query_file, prot
  banner "LOADING QUERY"
  hash = Hash.new # a hash just to check that there are no multiple entries in the query
  file = Bio::FastaFormat.open(query_file)            # http://bioruby.org/rdoc/Bio/FastaFormat.html
  file.each do |entry|
    if entry.isNucl? || (prot && entry.isProt?)
      # query = entry.definition.scan(/[^|]+/).first.split.first
      query = entry.definition.split(/[\|\ ]/).first
      # query = entry.definition
      if hash.has_key?(query)
        abort "Found multiple entries in #{query_file} with the same name: #{query}"
      end
      hash[query] = entry
    else
      abort "Query sequence isn't in a format I expected (nucl/prot)"
    end
  end
  return hash
end

# minimum_orf_length = 65
# (1..6).each do |frame|
#   translated = Bio::Sequence::NA.new(entry.seq).translate(frame)
#   translated.split(/\*/).each do |orf|
#     if orf.length > minimum_orf_length
#       target = "#{entry.definition.scan(/[^|]+/).first.split.first}:#{frame}"            
#       hash[target] = orf
#     end
#   end
# end

def load_target target_file, nucl
  banner "LOADING TARGET"
  hash = Hash.new # a hash just to check that there are no multiple entries in the target
  file = Bio::FastaFormat.open(target_file)            # http://bioruby.org/rdoc/Bio/FastaFormat.html
  file.each do |entry|
    if entry.isProt? || (nucl && entry.isNucl?)
      # target = entry.definition.scan(/[^|]+/).first.split.first
      target = entry.definition.split(/[\|\ ]/).first
      if hash.has_key?(target)
        abort "Found multiple entries in #{target_file} with the same name: #{target}"
      end
      hash[target] = entry
    else
      abort "Query sequence isn't in a format I expected (nucl/prot)"
    end
  end
  return hash
end

def make_databases(vars, opts)
  # make database from query file # # # # # # 
  #
  query_db = ""
  if opts.databases
    # query_db = "#{opts.databases}/query_#{File.basename(opts.query)}.udb"
    vars[:query_db] = "#{opts.databases}/#{File.basename(opts.query)}.udb"
    vars[:nc_query_db] = "#{opts.databases}/nc_#{File.basename(opts.query)}.udb"
  else
    # query_db = "#{opts.output}/query_#{File.basename(opts.query)}.udb"
    vars[:query_db] = "#{opts.output}/#{File.basename(opts.query)}.udb"
    vars[:nc_query_db] = "#{opts.output}/nc_#{File.basename(opts.query)}.udb"
  end

  if opts.target_nucl and !opts.protein
    if !File.exists?(vars[:nc_query_db])
      make_query_db = "#{vars[:makeblastdb]} #{opts.query} -quiet -output #{vars[:nc_query_db]} > #{opts.output}/db_query.log"
      banner "MAKE QUERY DATABASE" if opts.verbose
      `#{make_query_db}`
    end
  else
    if !File.exists?(vars[:query_db])
      make_orfs = opts.protein ? "cp #{opts.query} #{opts.output}/query_orfs.fa" : "#{vars[:findorfs]} #{opts.query} -quiet -output #{opts.output}/query_orfs.fa -xlat -orfstyle 7"
      make_query_db = "#{vars[:makeblastdb]} #{opts.output}/query_orfs.fa -quiet -output #{vars[:query_db]} > #{opts.output}/db_query.log"
      banner "MAKE QUERY DATABASE" if opts.verbose
      `#{make_orfs}`
      `#{make_query_db}`
    end
  end

  # make database from target file # # # # # # 
  #
  if opts.databases
    # target_db = "#{opts.databases}/target_#{File.basename(opts.target)}.udb"
    vars[:target_db] = "#{opts.databases}/#{File.basename(opts.target)}.udb"
    vars[:nc_target_db] = "#{opts.databases}/nc_#{File.basename(opts.target)}.udb"
  else
    vars[:target_db] = "#{opts.output}/#{File.basename(opts.target)}.udb"
    vars[:nc_target_db] = "#{opts.output}/nc_#{File.basename(opts.target)}.udb"
  end

  if opts.target_nucl and !opts.protein
    if !File.exists?(vars[:nc_target_db])
      make_target_db = "#{vars[:makeblastdb]} #{opts.target} -quiet -output #{vars[:nc_target_db]} > #{opts.output}/db_target.log"
      banner "MAKE TARGET DATABASE" if opts.verbose
      `#{make_target_db}`
    end
  else
    if !File.exists?(vars[:target_db])
      make_orfs = "cp #{opts.target} #{opts.output}/target_orfs.fa"
      make_target_db = "#{vars[:makeblastdb]} #{opts.output}/target_orfs.fa -quiet -output #{vars[:target_db]} > #{opts.output}/db_target.log"
      `#{make_orfs}`
      `#{make_target_db}`
    end
  end
end

def banner string
  w=60
  puts "#"*w
  n = ((w - string.length)*0.5).to_i
  puts " "*n << string << " "*n
  puts "#"*w
end

if __FILE__ == $0

  puts "RBUsearch"
  puts " version 0.0.4"

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # get inputs
  opts = Trollop::options do
    version "v0.0.3"
    banner <<-EOS
    rbusearch - reciprocal best usearch

    author: Chris Boursnell (cmb211@cam.ac.uk)
    ideas and help: Richard X Smith (rds45@cam.ac.uk)
    ideas and help: Steve Kelly (rds45@cam.ac.uk)

    EOS
    opt :query, "Query in nucleotide fasta format (required)", :required => true, :type => String
    opt :target, "Target in protein fasta format (required)", :required => true, :type => String
    opt :output, "Output directory (required)", :required => true, :type => String
    opt :prefix, "String to prefix to annotated output fasta file", :default => "out_", :type => String
    # opt :protein, "Query is protein format (don't translate)"
    opt :target_nucl, "The target is nucleotide fasta file"
    opt :cores, "Number of cores", :default => 2, :type => :int
    opt :databases, "Store databases in a different place", :type => String
    opt :nofasta, "Don't output a fasta file"
    opt :verbose, "Turn on verbose output"
  end

  Trollop::die :query, "must exist" if !File.exist?(opts[:query]) if opts[:query]
  Trollop::die :target, "must exist" if !File.exist?(opts[:target]) if opts[:target]

  if !Dir.exists?("#{opts.output}")
    puts "#{opts.output} output directory not found. Creating directory"
    cmd = "mkdir #{opts.output}"
    `#{cmd}`
  else
    puts "Output directory already exists. files will be overwritten (maybe)"
  end

  usearch = `which usearch`.strip
  vars = {
    :usearch => "#{usearch}",
    :blastx => "#{usearch} -ublast",
    :tblastn => "#{usearch} -ublast",
    :blastn => "#{usearch} -ublast",
    :findorfs => "#{usearch} -findorfs",
    :makeblastdb => "#{usearch} -makeudb_ublast",
    :print_lim => 40,
    :limit_evalue => 1e-4
  }

  blastx_output = "#{opts.output}/query_#{File.basename(opts.query).split(".")[0..-2].join(".")}-#{File.basename(opts.target).split(".")[0..-2].join(".")}_result.txt"
  tblastn_output = "#{opts.output}/target_#{File.basename(opts.target).split(".")[0..-2].join(".")}-#{File.basename(opts.query).split(".")[0..-2].join(".")}_result.txt"
  
  # MAKE DATABASES # # # # # #
  make_databases(vars, opts)
  # # # # # # # # # # # # # #
  
  if !opts.protein and opts.target_nucl # 
    # puts "blastn"
    query_blastn_cmd = "#{vars[:blastn]} #{opts.query} -quiet -db #{vars[:nc_target_db]} -threads #{opts.cores} --strand both"
    query_blastn_cmd << " -blast6out #{blastx_output} -maxhits #{vars[:print_lim]} -maxrejects 12 -evalue #{vars[:limit_evalue]}"

    target_blastn_cmd = "#{vars[:blastn]} #{opts.target} -quiet -db #{vars[:nc_query_db]} -threads #{opts.cores} --strand both"
    target_blastn_cmd << " -blast6out #{tblastn_output} -maxhits #{vars[:print_lim]} -maxrejects 12 -evalue #{vars[:limit_evalue]}"

    if !File.exists?("#{blastx_output}")
      banner "BLASTN #{File.basename(opts.query)} -> #{File.basename(opts.target)}" if opts.verbose
      puts query_blastn_cmd if opts.verbose
      blastx_out = `#{query_blastn_cmd}`
      puts "Done blastn 1" if opts.verbose
    end

    if !File.exists?("#{tblastn_output}")
      banner "BLASTN #{File.basename(opts.target)} -> #{File.basename(opts.query)}" if opts.verbose
      puts target_blastn_cmd if opts.verbose
      tblastn_out = `#{target_blastn_cmd}`
      puts "Done blastn 2" if opts.verbose
    end
  else
    # puts "blastx"
    # query -> target database
    blastxCmd = "#{vars[:blastx]} #{opts.query} -quiet -db #{vars[:target_db]} -threads #{opts.cores}"
    blastxCmd << " -blast6out #{blastx_output} -maxhits #{vars[:print_lim]} -maxrejects 12 -evalue #{vars[:limit_evalue]}"

    # target -> query database
    tblastnCmd = "#{vars[:tblastn]} #{opts.target} -quiet -db #{vars[:query_db]} -threads #{opts.cores}"
    tblastnCmd << " -blast6out #{tblastn_output} -maxhits #{vars[:print_lim]} -maxrejects 12 -evalue #{vars[:limit_evalue]}"

    if !File.exists?("#{blastx_output}")
      banner "BLASTX #{File.basename(opts.query)} -> #{File.basename(opts.target)}" if opts.verbose
      puts blastxCmd if opts.verbose
      blastx_out = `#{blastxCmd}`
      puts "Done blastx" if opts.verbose
    end

    if !File.exists?("#{tblastn_output}")
      banner "TBLASTN #{File.basename(opts.target)} -> #{File.basename(opts.query)}" if opts.verbose
      puts tblastnCmd if opts.verbose
      tblastn_out = `#{tblastnCmd}`
      puts "Done tblastn" if opts.verbose
    end
  end

  # # # # # # # # # # # #
  # create output file  #
  # # # # # # # # # # # # 
  reciprocals = Hash.new
  if !File.exists?("#{opts.output}/reciprocal_hits.txt")
    puts "reading query output file" if opts.verbose
    query_results = Hash.new

    lines = `wc -l #{blastx_output}`
    lines = lines.to_i
    count=0
    File.open("#{blastx_output}").each_line do |line|
      line.chomp!
      if count % 15000==0
        print "#{100*count/lines}%.." if lines>0 if opts.verbose
      end
      cols = line.split("\t")
      # id = cols[0].scan(/[^|]+/).first.split.first
      id = cols[0].split(/[\|\ ]/).first
      if !query_results.has_key?(id)
        query_results[id] = []
      end
      query_results[id] << Hit.new(cols)
      count+=1
    end
    puts "Done" if opts.verbose
    # - - - 
    puts "reading target output file" if opts.verbose
    target_results = Hash.new

    lines = `wc -l #{blastx_output}`
    lines = lines.to_i
    count=0
    File.open("#{tblastn_output}").each_line do |line|
      line.chomp!
      if count % 15000==0
        print "#{100*count/lines}%.." if lines>0 if opts.verbose
      end
      cols = line.split("\t")
      # id = cols[0].scan(/[^|]+/).first.split.first
      id = cols[0].split(/[\|\ ]/).first
      # id = cols[0]
      if !target_results.has_key?(id)
        target_results[id] = []
      end
      target_results[id] << Hit.new(cols)
      count+=1
    end

    # p query_results
    # puts "----------------------"
    # p target_results

    evalues = []
    missed_evalues = []
    missed = Hash.new
    query_results.each_pair do |query, list_of_hits|               # for each query sequence
      # puts "query : #{query}"
      best_hit = list_of_hits[0]                                   # get the first hit
      # puts "best_hit : #{best_hit}"
      if target_results.has_key?(best_hit.target)

        list = target_results[best_hit.target]
        best_reciprocal_hit = list[0]
        # puts "best reciprocal : #{best_reciprocal_hit}"
        e = best_reciprocal_hit.evalue.to_f
        e = 1e-300 if e==0
        e = -Math.log10(e)
        if best_reciprocal_hit.target == query
          # puts "    #{query} => #{best_reciprocal_hit.target} "
          reciprocals[best_hit.query] = best_hit
          evalues << {:e => e, :length => best_reciprocal_hit.alnlen} if e<300    # store evalues and lengths of reciprocal best hits to plot
        else
          # puts "#{best_reciprocal_hit.target} != #{query}"
          missed[best_hit.query] = best_hit
          missed_evalues << {:e => e, :length => best_reciprocal_hit.alnlen}     # store evalues and lengths of reciprocal best hits to plot
        end
      end
    end
    # puts "There are : #{reciprocals.keys.length} data points"
    # puts "There are : #{missed.keys.length} missed data points"
    # exit
    e_data = "cond\tevalue\tlen\n"
    evalues.each do |h|
      e_data << "A\t#{h[:e]}\t#{h[:length]}\n"
    end

    File.open("#{opts.output}/e_data.txt", "w") { |io|  io.write(e_data)}
    
    puts "starting line fitting" if opts.verbose
    length_hash = Hash.new
    evalues.each do |h|
      if !length_hash.has_key?(h[:length])
        length_hash[h[:length]]=[]
      end
      length_hash[h[:length]] << h
    end
   
    fitting=Hash.new
    (10..15000).each do |centre|
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

    # export list of reciprocal hits
    output = ""
    output_2 = ""
    reciprocals.each_pair do |id,hit|
      output << "#{hit.query}\t#{hit.target}\t#{hit.id}\t#{hit.alnlen}\t#{hit.mismatches}\t#{hit.gaps}\t#{hit.qstart}\t#{hit.qend}\t#{hit.tstart}\t#{hit.tend}\t#{hit.evalue}\t#{hit.bitscore}\t1\n"
    end

    # use fitting to find hits that weren't that reciprocal
    missed.each_pair do |key, hit|
      l = hit.alnlen.to_i
      e = hit.evalue
      e = 1e-300 if e==0
      e = -Math.log10(e)

      if fitting.has_key?(l)
        if e > fitting[l]
          output << "#{hit.query}\t#{hit.target}\t#{hit.id}\t#{hit.alnlen}\t#{hit.mismatches}\t#{hit.gaps}\t#{hit.qstart}\t#{hit.qend}\t#{hit.tstart}\t#{hit.tend}\t#{hit.evalue}\t#{hit.bitscore}\t2\n"
          # add these ones to reciprocals so they can be printed to the txt file - oh yeah shit i didn't do this
          if !reciprocals.has_key?(key)
            reciprocals[key] = hit
          end
        else
          output_2 << "#{hit.query}\t#{hit.target}\t#{hit.id}\t#{hit.alnlen}\t#{hit.mismatches}\t#{hit.gaps}\t#{hit.qstart}\t#{hit.qend}\t#{hit.tstart}\t#{hit.tend}\t#{hit.evalue}\t#{hit.bitscore}\t3\n"
        end
      end
    end

    File.open("#{opts.output}/reciprocal_hits.txt", "w") { |io| io.write(output) }
    File.open("#{opts.output}/not_reciprocal_hits.txt", "w") { |io| io.write(output_2)}

  end
  if !opts.nofasta
    # export new fasta file of contigs with updated names
    puts "creating new fasta file" if opts.verbose
    fasta=""
    File.open(opts.query, "r").each_line do |line|
      if line=~/^>(\S+)/
        name=$1
        if reciprocals.has_key?(name)
          name = "#{name}_#{reciprocals[name].target}"
        end
        fasta << ">#{name}\n"
      else
        fasta << line
      end
    end
    puts "writing new fasta file" if opts.verbose
    File.open("#{opts.output}/#{opts.prefix}annotated.fasta", "w") { |io|  io.write(fasta)}
  end
  # write a log about the command run for this annotation
  #
  log = "rbusearch\n"
  log << "query\t#{opts.query}\n"
  log << "target\t#{opts.target}\n"
  log << "output\t#{opts.output}\n"
  log << "prefix\t#{opts.prefix}\n"

  File.open("#{opts.output}/#{opts.prefix}rbusearch.log", "w") {|io| io.write log}
 

end # if __FILE__ == $0