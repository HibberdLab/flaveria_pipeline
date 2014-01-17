#!/usr/bin/env ruby

#
# recursive round robin
#
# this is different from the previous version in the way that the hits are stored
#
# now the data from 'reciprocal hits.txt' files is stored in a tree
#
# 
#

class Node 
  attr_accessor :name, :agi, :bitscore

  def initialize(name, agi, bitscore)
    @name = name
    @agi = agi
    @bitscore = bitscore
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
    return "#{@query}\t#{@target}\t#{@id}\t#{@alnlen}\t#{@evalue}\t#{@bitscore}\ttype:#{@type}\tpath:#{@path}"
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
ideas and help: Richard X Smith (rds45@cam.ac.uk)


EOS
  opt :reference, "Annotated reference protein fasta file", :required => true, :type => String
  opt :list, "List of nucleotide fasta files to annotate", :required => true, :type => String
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

# p list

cpu = 7
nodes = Hash.new
g = RGL::DirectedAdjacencyGraph.new

# run rbusearch
list.each do |file_query|
  $stderr.puts "Doing rbusearch on #{file_query} against #{opts.reference}"
  dir_output = "#{(File.basename(file_query)).split(".").first}_into_#{(File.basename(opts.reference)).split(".").first}"

  query = File.basename(file_query).split(".").first    
  rbu_cmd = "ruby rbusearch.rb "
  rbu_cmd << "--query #{file_query} --target #{opts.reference} --output #{dir_output} "
  rbu_cmd << "--databases . --prefix #{query} --cores #{cpu} --nofasta"
  puts rbu_cmd if opts.verbose
  # `#{rbu_cmd}`

  list.each do |file_target|
    if file_query != file_target
      $stderr.puts "Doing rbusearch of #{file_query} against #{file_target}"
      dir_output = "#{(File.basename(file_query)).split(".").first}_into_#{(File.basename(file_target)).split(".").first}"
      query = File.basename(file_query).split(".").first
      rbu_cmd = "ruby rbusearch.rb "
      rbu_cmd << "--query #{file_query} --target #{file_target} --target-nucl --output #{dir_output} "
      rbu_cmd << "--databases . --prefix #{query} --cores #{cpu} --nofasta"
      puts rbu_cmd if opts.verbose
      # `#{rbu_cmd}`
    end
  end
end

# scan output reciprocal hits files
list.each do |file_query|
  query_name = File.basename(file_query).split(".").first
  target_name = File.basename(opts.reference).split(".").first
  dir_output = "#{query_name}_into_#{target_name}"
  Dir.chdir(dir_output) do |dir|
    File.open("reciprocal_hits.txt").each_line do |line|
      cols = line.split("\t")
      hit = Hit.new(cols,query_name, target_name)

      # here we create a node based on the query name
      # and specify the AGI of the node based on the
      # target_name
      name = "#{hit.query}"
      if !nodes.has_key?(name)
        node = Node.new(name, hit.target, hit.bitscore)
        nodes[name]=node
      end
    end
  end

  list.each do |file_target|
    if file_query != file_target
      query_name = File.basename(file_query).split(".").first
      target_name = File.basename(file_target).split(".").first
      puts "Doing rbusearch of #{file_query} against #{file_target}"
      dir_output = "#{query_name}_into_#{target_name}"
      Dir.chdir(dir_output) do |dir|
        File.open("reciprocal_hits.txt").each_line do |line|
          cols = line.split("\t")
          hit = Hit.new(cols,query_name, target_name)
        
          # here we create 2 nodes based on the query_name
          # and the target_name, and then create an edge
          # # going from query to target
          # name1 = "#{hit.query}"
          # name2 = "#{hit.target}"
          # node1 = Node.new(name1, nil, nil)
          # node2 = Node.new(name2, nil, nil)
          # nodes << node1
          # nodes << node2
          # g.add_edge(node1, node2)


          name1 = "#{hit.query}"
          name2 = "#{hit.target}"
          node1 = Node.new(name1, nil, nil)
          node2 = Node.new(name2, nil, nil)
          if nodes.has_key?(name1)
            node1 = nodes[name1]
          else
            nodes[name1]=node1
          end

          if nodes.has_key?(name2)
            node2 = nodes[name2]
          else
            nodes[name2]=node2
          end

          g.add_edge(node1, node2)

        end
      end
    end
  end
end

(1..list.length).each do |iteration|
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
  agi = node.agi.split(":").last if node.agi!=nil
  annotation[species][contig]=agi
end

annotation.each_pair do |species, hash2|
  hash2.each_pair do |contig, agi|
    puts "#{species}\t#{contig}\t#{agi}"
  end
end
