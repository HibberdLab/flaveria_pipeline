#!/usr/bin/env ruby

#
# recursive round robin (blast version)
# 
#

class Node 
  attr_accessor :name, :agi, :bitscore

  def initialize(name, agi, bitscore)
    @name = name
    @agi = agi
    @bitscore = bitscore
    # TODO add a count
    # TODO add a original of agi thingy
  end

  def to_s
    "#{@name} #{@agi} #{@bitscore}"
  end

  def <=> other
    return 0 if @name==other.name
    return 1 if @name>other.name
    return -1 if @name<other.name
  end

  def == other
    if @name == other.name
      return true
    else
      return false
    end
  end
end

class Hit
  # Fields: query id, subject id, % identity, alignment length, mismatches,
  # gap opens, q. start, q. end, s. start, s. end, evalue, bit score
  attr_accessor :query, :target, :id, :alnlen, :mismatches, :gaps, :qstart, :qend, :tstart, :tend, :evalue, :bitscore
  def initialize(list, query_name, target_name)
    @query      = "#{query_name}:#{list[0].split(/\|/).first}"
    @target     = "#{target_name}:#{list[1].split(/\|/).first}"
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
    @type       = list[12].to_i
  end

  def to_s
    return "#{@query}\t#{@target}\t#{@id}\t#{@alnlen}\t#{@evalue}\t#{@bitscore}"
  end
end

require 'bio'
require 'trollop'
require 'rgl/adjacency'
require 'rgl/bidirectional'

opts = Trollop::options do
  version "v0.0.2"
  banner <<-EOS
round robin annotation: 

the idea is the more you have to annotate the better it gets

author: Chris Boursnell (cmb211@cam.ac.uk)
ideas and help: Richard X Smith-Unna (rds45@cam.ac.uk)

EOS
  opt :reference, "Annotated reference protein fasta file", :required => true, :type => String
  opt :list, "File containing list of nucleotide fasta files to annotate", :required => true, :type => String
  opt :threads, "Threads", :required => true, :type => :int
  opt :output, "Final annotation output file", :required => true, :type => String
  opt :test, "Don't actually do anything"
  opt :verbose, "Be Verbose"
end

Trollop::die :reference, "must exist" if !File.exist?(opts[:reference]) if opts[:reference]
Trollop::die :list, "must exist" if !File.exist?(opts[:list]) if opts[:list]

list = File.readlines("#{opts.list}")
list.each do |file|
  file.chomp!
  if !File.exists?(file)
    abort "Can't find #{file}"
  end
end

cpu = opts.threads
nodes = Hash.new
g = RGL::DirectedAdjacencyGraph.new
# rbusearch = "/home/cmb211/scripts/flaveria_pipeline/rbusearch.rb"
rbblast = "/home/cmb211/scripts/flaveria_pipeline/rbblast.rb"

# run rbusearch
list.each do |file_query|
  $stderr.puts "Doing rbblast on #{file_query} against #{opts.reference}"
  dir_output = "#{(File.basename(file_query)).split(".").first}_into_#{(File.basename(opts.reference)).split(".").first}"

  query = File.basename(file_query).split(".").first
  rbu_cmd = "ruby #{rbblast} "
  rbu_cmd << "--query #{file_query} --target #{opts.reference} --output #{dir_output} "
  rbu_cmd << " --prot --cores #{cpu} --nofasta --verbose"
  puts rbu_cmd if opts.verbose
  `#{rbu_cmd}` if !opts.test

  list.each do |file_target|
    if file_query != file_target
      $stderr.puts "Doing rbblast of #{file_query} against #{file_target}"
      dir_output = "#{(File.basename(file_query)).split(".").first}_into_#{(File.basename(file_target)).split(".").first}"
      query = File.basename(file_query).split(".").first
      rbu_cmd = "ruby #{rbblast} "
      rbu_cmd << "--query #{file_query} --target #{file_target} --nucl --output #{dir_output} "
      rbu_cmd << " --cores #{cpu} --nofasta --verbose"
      puts rbu_cmd if opts.verbose
      `#{rbu_cmd}` if !opts.test
    end
  end
end

# edges = File.open("edges.txt", "w")
# scan output reciprocal hits files
list.each do |file_query|
  query_name = File.basename(file_query).split(".").first
  target_name = File.basename(opts.reference).split(".").first
  dir_output = "#{query_name}_into_#{target_name}"
  Dir.chdir(dir_output) do |dir|
    if File.exists?("reciprocal_hits.txt")
      puts "opening #{dir_output}" if opts.verbose
      File.open("reciprocal_hits.txt").each_line do |line|

        cols = line.split("\t")
        hit = Hit.new(cols,query_name, target_name)

        # here we create a node based on the query name
        # and specify the AGI of the node based on the
        # target_name
        name = "#{hit.query}"

        if !nodes.has_key?(name)
          node = Node.new(name, hit.target, hit.bitscore)
          # if line =~ /AT5G49560/ or line=~/spinosa_00056/
          #   puts "line: #{line}"
          #   puts "hit:  #{hit}"
          #   puts "name: #{name}"
          #   puts "node: #{node}"
          # end
          nodes[name]=node
        else
          if nodes[name].agi == nil
            nodes[name].agi = hit.target
            nodes[name].bitscore = hit.bitscore
            # puts "setting node #{name} agi to #{hit.target} and bitscore to #{hit.bitscore}" if name =~ /spinosa_00056/
          else
            puts "this shouldn't happen"
          end
        end
      end
    end
  end

  # scan output reciprocal hit files 

  list.each do |file_target|
    if file_query != file_target
      query_name = File.basename(file_query).split(".").first
      target_name = File.basename(file_target).split(".").first
      dir_output = "#{query_name}_into_#{target_name}"
      Dir.chdir(dir_output) do |dir|
        if File.exists?("reciprocal_hits.txt")
          puts "opening #{dir}" if opts.verbose
          File.open("reciprocal_hits.txt").each_line do |line|
            cols = line.split("\t")
            hit = Hit.new(cols,query_name, target_name)

            name1 = "#{hit.query}"
            name2 = "#{hit.target}"
            # if name1 =~ /spinosa_00056/
            #   puts "name1:"
            #   puts "line: #{line}"
            # end
            # if name2 =~ /spinosa_00056/
            #   puts "name2:"
            #   puts "line: #{line}"
            # end

            node1 = Node.new(name1, nil, nil)
            node2 = Node.new(name2, nil, nil)
            if nodes.has_key?(name1)
              node1 = nodes[name1]
            else
              nodes[name1] = node1
            end

            if nodes.has_key?(name2)
              node2 = nodes[name2]
            else
              nodes[name2] = node2
            end

            g.add_edge(node1, node2)
            # TODO add these edges to a csv file for visualisation
            # if node1.name == "C.spinosa_00056" or node2.name == "C.spinosa_00056" or node1.name == "AT5G49560.1" or node2.name == "AT5G49560.1" or node1.name == "AT4G36630.1" or node2.name == "AT4G36630.1"
            # edges.write("#{node1.name}\t#{node1.bitscore}\t#{node1.agi}\t#{node2.name}\t#{node2.bitscore}\t#{node2.agi}\n")
            # end
          end
        end
      end
    end
  end
end

# cascade best annotations through the graph
puts "cascading best annotations through the graph" if opts.verbose
(1..list.length).each do |iteration|
  File.open("nodes-#{iteration}.txt", "w") do |io|
    nodes.each_pair do |name, node|
      io.write("#{node.name}\t#{node.agi}\t#{node.bitscore}\n")
    end
  end
  nodes.each_pair do |name,node|
    if node.agi == nil
      neighbours = g.adjacent_vertices(node)
      bitscore=0
      agi=nil
      neighbours.each do |n|
        if n.bitscore!=nil and n.bitscore > bitscore
          agi=n.agi
          bitscore=n.bitscore
        end
      end
      node.agi = agi 
      node.bitscore = bitscore
    end
  end
end

annotation=Hash.new
nodes.each_pair do |name, node|
  (species, contig) = name.split(":")
  if !annotation.has_key?(species)
    annotation[species]=Hash.new
  end
  agi = node.agi.split(":").last if node.agi != nil
  annotation[species][contig]=agi if agi != "" and agi!=nil
end

# need to make the output a bit better - TODO [ ]

output=""
annotation.each_pair do |species, hash2|
  hash2.each_pair do |contig, agi|
    output << "#{species}\t#{contig}\t#{agi}\n"
  end
end
puts "writing output to #{opts.output}" if opts.verbose
File.open("#{opts.output}", "w") {|io| io.write(output)}

