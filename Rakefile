require 'bio'

required = {
  :input_reads => "raw_data",
  :quality => "quality_check",
  :trimmed_reads => "trimmed_reads",
  :yaml => "datasets",
  :corrected_reads => "corrected_reads",
  :single_corrected_reads => "single_corrected_reads",
  :hammer_log => "hammer.log",
  :khmered_reads => "khmered_reads",
  :config => "soapdt.config",
  :soap_output => "soap_output",
  :sga_output => "sga_output",
  :idba_output => "idba_output",
  :combined_assembly => "combined_assembly",
  :annotation_output => "annotation_summary",
  :expression_output => "expression_output",
  :bowtie_index => "bowtie_index"
}

threads = 22
memory = 92
maximum_files_to_hammer_at_a_time = 5
fastqc_path = "/applications/fastqc_v0.10.1/FastQC/fastqc"
hammer_path = "~/apps/SPAdes-2.5.1-Linux/bin/spades.py"
trimmomatic_path = "/home/cmb211/apps/Trimmomatic-0.32/trimmomatic-0.32.jar"
khmer_path = "/home/cmb211/.local/bin/normalize-by-median.py"
protein_reference = "/home/cmb211/flaveria/at_p.faa"
idba = "/home/cmb211/apps/idba_tran-1.0.13/bin/idba_tran"
cd_hit_est = "/home/cmb211/apps/cd-hit-v4.6.1-2012-08-27/cd-hit-est"
gapcloser = "/home/cmb211/bin/GapCloser"
lcs = ""
path = ""

task :input do
  puts "checking input"
  a=[]
  File.open(required[:input_reads], "r").each_line do |line|
    line.chomp!
    if File.exists?("#{line}")
      a << File.basename(line.gsub(".fastq", ""))
      path = File.dirname(line)
    else
      abort "Can't find #{line}"
    end
  end
  abort "Can't find protein reference for annotation" if !File.exists?(protein_reference)
  s = a.min_by(&:size)
  lcs = catch(:hit) {  s.size.downto(1) { |i| (0..(s.size - i)).each { |l| throw :hit, s[l, i] if a.all? { |item| item.include?(s[l, i]) } } } }
  lcs = "out" if lcs.length == 0
  #lcs=lcs[0..-2] if lcs[-2]=~/\_\.\,\-/
end

file required[:quality] => required[:input_reads] do
  puts "running fastqc reads..."
  if !Dir.exists?("#{path}/fastqc_output")
    `mkdir #{path}/fastqc_output`
  end

  input = []
  File.open("#{required[:input_reads]}", "r").each_line do |line|
    line.chomp!
    input << line
  end

  files = ""
  Dir.chdir("#{path}/fastqc_output") do
    input.each do |file|
      filename = File.basename(file)
      fastqc_line = filename.split(".")[0..-2].join(".") + "_fastqc"
      if File.exists?("#{fastqc_line}")
        #puts "Found #{fastqc_line}"
      else
        #puts "Didn't find #{fastqc_line}"
        files << " #{file} "
      end
    end
  end
  
  quality_threads = threads > 5 ? 5 : threads # limits to 5 threads due to io
  `#{fastqc_path} --kmers 5 --threads #{quality_threads} --outdir #{path}/fastqc_output #{files}` if files.length>0
  
  fastqc_summary = Hash.new
  Dir.chdir("#{path}/fastqc_output") do
    Dir["*fastqc"].each do |fastqc_dir|
      Dir.chdir("#{fastqc_dir}") do
        File.open("fastqc_data.txt", "r").each_line do |line|
          if line=~/Overrepresented\s+sequences\s+(\S+)/
            fastqc_summary[fastqc_dir] = $1
          end
        end
      end
    end
  end
  warn=0
  File.open("#{required[:quality]}", "w") do |out|
    fastqc_summary.each_pair do |key, value|
      out.write "#{key}\t#{value}\n"
      warn+= 1 if value=="warn"
    end
  end
  
  puts "FastQC: There were #{warn} warnings about overrepresented sequences"
end

file required[:trimmed_reads] => required[:input_reads] do
  puts "creating trimmed reads..."
  trim_threads = threads > 4 ? 4 : threads
  trim_batch_cmd = "ruby trim-batch.rb "
  trim_batch_cmd << "--jar #{trimmomatic_path} "
  trim_batch_cmd << "--pairedfile #{required[:input_reads]} "
  #trim_batch_cmd << "--singlefile #{required[:single_input_reads]} "
  trim_batch_cmd << "--threads #{trim_threads} " # due to io limitations this is capped at 4
  trim_batch_cmd << "--quality 15 "
  # puts trim_batch_cmd
  `#{trim_batch_cmd}`
  list_of_trimmed_reads = ""
  File.open("#{required[:input_reads]}", "r").each_line do |line|
    line.chomp!
    
    filename = File.basename(line)
    if File.exists?("#{path}/t.#{filename}")
      list_of_trimmed_reads << "#{path}/t.#{filename}\n"
    end
    if File.exists?("#{path}/t.#{filename}U")
      # have to rename the files because spades only works with files with known extensions eg fasta, fa, fq, fastq...
      rename =  "mv #{path}/t.#{filename}U #{path}/tU.#{filename}"
      # puts rename
      `#{rename}`
      list_of_trimmed_reads << "#{path}/tU.#{filename}\n"
    end
  end
  File.open("#{required[:trimmed_reads]}", "w") do |out|
    out.write list_of_trimmed_reads
  end
end

file required[:yaml] do # construct dataset.yaml file for bayeshammer input
  hash = Hash.new
  hash[:left]=[]
  hash[:right]=[]
  hash[:single]=[]
  File.open("#{required[:trimmed_reads]}", "r").each_line do |line|
    line.chomp!
    
    filename = File.basename(line)
    if filename=~/^t\..*R1.*/
      hash[:left] << line
    elsif filename=~/^t\..*R2.*/
      hash[:right] << line
    else
      hash[:single] << line
    end
  end
  datasets=[]
  #puts "files = #{hash[:left].length}"
  #puts "max = #{maximum_files_to_hammer_at_a_time}"
  #puts (hash[:left].length.to_f / maximum_files_to_hammer_at_a_time.to_f).ceil
  chunk = (hash[:left].length.to_f / maximum_files_to_hammer_at_a_time.to_f).ceil
  
  (1..chunk).each do |i|
    left=[]
    right=[]
    single=[]
    (1..maximum_files_to_hammer_at_a_time).each do
      left << hash[:left].shift
      right << hash[:right].shift
      single << hash[:single].shift
      single << hash[:single].shift
    end
    left = left.delete_if {|x| x==nil }
    right = right.delete_if {|x| x==nil }
    single = single.delete_if {|x| x==nil}
    yaml = "[\n"
    yaml << "  {\n"
    yaml << "    orientation: \"fr\",\n"
    yaml << "    type: \"paired-end\",\n"
    yaml << "    left reads: [\n"
    left.each_with_index do |left_read, j|
      yaml << "      \"#{left_read}\""
      yaml << "," if j < left.length-1
      yaml << "\n"
    end
    yaml << "    ],\n"
    yaml << "    right reads: [\n"
    right.each_with_index do |right_read, j|
      yaml << "      \"#{right_read}\""
      yaml << "," if j < right.length-1
      yaml << "\n"
    end
    yaml << "    ],\n"
    yaml << "  },\n"
    yaml << "  {\n"
    yaml << "    type: \"single\",\n"
    yaml << "    single reads: [\n"
    single.each_with_index do |single_read, j|
      yaml << "      \"#{single_read}\""
      yaml << "," if j < single.length-1
      yaml << "\n"
    end
    yaml << "    ]\n"
    yaml << "  }\n"
    yaml << "]\n"
    datasets << "dataset_#{i}.yaml"
    File.open("dataset_#{i}.yaml", "w") do |out|
      out.write yaml
    end
  end # end of chunk
  File.open("#{required[:yaml]}", "w") do |yaml_out|
    datasets.each do |yaml_file|
      yaml_out.write "#{yaml_file}\n" 
    end
  end
end

file required[:corrected_reads] => required[:trimmed_reads] do
  puts "running bayeshammer to correct reads..."
  
  count=0
  output_directories=[]

  write_hammer_to_tmp = true
  if write_hammer_to_tmp
    node=`hostname`
    if node=~/node9/
      hammer_out_path = "/disk2/tmp/cmb211" # node9 tmp cd /disk2/tmp/cmb211
    elsif node=~/node8/
      hammer_out_path = "/tmp/cmb211" # node8 tmp
    end
  else
    hammer_out_path = path
  end
  puts "hammer path = #{hammer_out_path}"
  File.open("#{required[:yaml]}", "r").each_line do |dataset_line|
    puts "running hammer on #{dataset_line}"
    dataset_line.chomp!
    cmd = "python #{hammer_path} --dataset #{dataset_line} --only-error-correction --disable-gzip-output -m #{memory} -t #{threads} -o #{hammer_out_path}/output_#{count}.spades"
    output_directories << "#{hammer_out_path}/output_#{count}.spades"
    if !File.exists?("#{hammer_out_path}/output_#{count}.spades")
      puts cmd
      hammer_log = `#{cmd}`
      File.open("#{path}/hammer_#{dataset_line}.log", "w") {|out| out.write hammer_log}
    else
      puts "output_#{count}.spades already exists, not running hammer on this batch."
      puts "if you wanted it to run, move or delete this directory"
    end
    count+=1
  end

  paired = []
  single = []

  output_directories.each do |dir|
    Dir.chdir(dir) do
      Dir.chdir("corrected") do
        fastq_files = Dir["*fastq"]
        abort "Something went wrong with BayesHammer and no corrected reads were created in #{dir}" if fastq_files.length ==0
        fastq_files.each do |fastq|
          if fastq =~ /t\..*R[12].*fastq/
            paired << "#{dir}/corrected/#{fastq}" # #{path}/
          elsif fastq =~ /tU.*fastq/
            single << "#{dir}/corrected/#{fastq}" # #{path}/
          end
        end
      end
    end
  end
  paired.sort!
  File.open("#{required[:corrected_reads]}", "w") do |out|
    paired.each do |pe|
      out.write "#{pe}\n"
    end
  end
  File.open("single_#{required[:corrected_reads]}", "w") do |out|
    single.each do |sng|
      out.write "#{sng}\n"
    end
  end
end

file required[:khmered_reads] => required[:corrected_reads] do
  puts "running khmer to reduce coverage of reads..."

  # make a list of the paired corrected reads
  filelist=[]
  File.open("#{required[:corrected_reads]}", "r").each_line do |line|
    line.chomp!
    abort "couldn't find file" if !File.exists?(line)
    filelist << line
  end
  
  # interleave the corrected fastq files together into #{path} with a .in suffix
  interleaved_files=[]
  interleave_out = path
  filelist.each_slice(2) do |pair|
    left = pair[0]
    right = pair[1]
    puts "interleaving #{File.basename(left)} and #{File.basename(right)}"
    base = File.basename(left)
    cmd = "paste #{left} #{right} | paste - - - - | "
    cmd << " awk -v FS=\"\t\" -v OFS=\"\n\" \'{print(\"@read\"NR\":1\",$3,$5,$7,\"@read\"NR\":2\",$4,$6,$8)}\' > #{interleave_out}/#{base}.in"
    if !File.exists?("#{interleave_out}/#{base}.in")
      # puts cmd
      `#{cmd}`
    end
    interleaved_files << "#{interleave_out}/#{base}.in"
    # rm1 = "rm #{File.basename(left)}"
    # rm2 = "rm #{File.basename(right)}"
    # `#{rm1}`
    # `#{rm2}`
  end

  # settings for khmer
  first = true
  kmer_size = 23
  n = 4
  khmer_memory = 48
  cutoff = 20
  x = (khmer_memory/n*1e9).to_i

  # run khmer on the interleaved files and export .keep files into #{path}
  interleaved_files.each do |file|
    Dir.chdir(File.dirname(file)) do
      if first
        cmd = "#{khmer_path} -p -k #{kmer_size} -C #{cutoff} -N #{n} -x #{x} --savehash table.kh #{file}"
        puts "#{cmd}"
        puts `#{cmd}`
        first = false
      else
        cmd = "#{khmer_path} -p -k #{kmer_size} -C #{cutoff} -N #{n} -x #{x} --loadhash table.kh --savehash table2.kh #{file}"
        puts "#{cmd}"
        puts `#{cmd}`
        `mv table2.kh table.kh`
      end
    end
  end

  # cat all the .keep files together from the paired khmer run
  cat_cmd = "cat "
  Dir["#{path}/*keep"].each do |file|
    cat_cmd << "#{file} "
  end
  cat_cmd << " > #{path}/#{lcs}.khmered.fastq"
  # puts cat_cmd
  `#{cat_cmd}`

  # cat all the single corrected reads into a single fastq file to run khmer on them
  cat_cmd = "cat "
  corrected_path=""
  File.open("single_#{required[:corrected_reads]}", "r").each_line do |line|
    line.chomp!
    cat_cmd << "#{line} "
  end
  cat_cmd << " > #{path}/#{lcs}.single.fastq"
  # puts cat_cmd
  `#{cat_cmd}`

  # run khmer on all the single reads
  Dir.chdir(path) do
    cmd = "#{khmer_path} -k #{kmer_size} -N #{n} -x #{x} #{path}/#{lcs}.single.fastq"
    puts "#{cmd}"
    puts `#{cmd}`
  end
  # remove lcs.single.fastq as it's just a cat of other files and it's big
  rm_cmd = "rm #{path}/#{lcs}.single.fastq"
  `#{rm_cmd}`
  
  # deinterleave the output from adds left.fastq and right.fastq to the end of the filename
  cmd = "ruby smart_deinterleave.rb -f #{path}/#{lcs}.khmered.fastq -o #{path}/#{lcs}" 
  puts cmd
  `#{cmd}`

  # write names of files to khmered_reads
  File.open("#{required[:khmered_reads]}", "w") do |out|
    out.write "#{path}/#{lcs}.left.fastq\n"
    out.write "#{path}/#{lcs}.right.fastq\n"
    out.write "#{path}/#{lcs}.single.fastq.keep\n"
  end
end

file required[:config] do # construct dataset.yaml file for bayeshammer input
  reads  = File.open("#{required[:khmered_reads]}", "r")
  left  = reads.readline.chomp!
  right = reads.readline.chomp!
  single = reads.readline.chomp!
  config = "max_rd_len=20000\n"
  config << "[LIB]\n"
  config << "avg_ins=250\n"
  config << "reverse_seq=0\n"
  config << "asm_flags=3\n"
  # config << "rank=2\n"
  config << "q1=#{left}\n"
  config << "q2=#{right}\n"
  config << "q=#{single}\n"
  File.open("soapdt.config", "w") {|out| out.write config}
end

file required[:soap_output] => required[:khmered_reads] do
  puts "running soap on khmered reads..."

  if !Dir.exists?("#{path}/soap")
    mkdir_cmd = "mkdir #{path}/soap"
    `#{mkdir_cmd}`
  end

  soap_cmd = "SOAPdenovo-Trans-127mer all -s #{required[:config]} -o #{path}/soap/#{lcs}soap -p #{threads}"
  puts soap_cmd
  `#{soap_cmd}`

  # run soap gapcloser on scafseq
  gapcloser_cmd = "#{gapcloser} -a #{path}/soap/#{lcs}soap.scafSeq -b soapdt.config -o #{path}/soap/#{lcs}soap.gapcloser.fasta" # GapCloser -a Ft_soap.scafSeq -b soapdt.config -o Ft_soap.filled.fasta -l 101
  puts gapcloser_cmd
  `#{gapcloser_cmd}`

  # run sga gap filler on *.scafSeq
  pre_cmd =     "sga preprocess -p 1 -o #{path}/soap/#{lcs}cleaned.fastq #{path}/#{lcs}.left.fastq #{path}/#{lcs}.right.fastq" 
  index_cmd =   "sga index -p #{path}/soap/#{lcs} -t #{threads} -a ropebwt #{path}/soap/#{lcs}cleaned.fastq"
  gapfill_cmd = "sga gapfill -p #{path}/soap/#{lcs} -e 47 -x 1 -t #{threads} -o #{path}/soap/#{lcs}soap.gapfill.fasta #{path}/soap/#{lcs}soap.gapcloser.fasta"

  puts pre_cmd
  `#{pre_cmd}`

  puts index_cmd
  `#{index_cmd}`

  puts gapfill_cmd
  `#{gapfill_cmd}`

  # make a histogram of the lengths of the contigs in the output assembly file
  #
  if File.exists?("#{path}/soap/#{lcs}soap.gapfill.fasta") # soapdtgraph.scafSeq
    abort "soap assembly file is empty!" if File.zero?("#{path}/soap/#{lcs}soap.gapfill.fasta")
    contigs = Bio::FastaFormat.open("#{path}/soap/#{lcs}soap.gapfill.fasta")
    histogram = {}
    contigs.each do |entry|
      bucket = (entry.seq.length*0.01).round
      if !histogram.has_key?(bucket)
        histogram[bucket]=0
      end
      histogram[bucket]+=1
    end
    File.open("#{required[:soap_output]}", "w") do |out|
      keys = histogram.keys.sort
      keys.each do |key|
        out.write "#{100*key}\t#{histogram[key]}\n"
      end
    end
  else
    abort "couldn't find soap output file"
  end
end

file required[:sga_output] => required[:khmered_reads] do
  puts "running sga on khmered reads..."

  if !Dir.exists?("#{path}/sga")
    mkdir_cmd = "mkdir #{path}/sga"
    `#{mkdir_cmd}`
  end

  files = File.readlines("#{required[:khmered_reads]}").map{|n| n.chomp!}.compact
  left = files[0]
  right = files[1]
  sga_cmd = "ruby sga.rb --verbose --left #{left} --right #{right} --output #{path}/sga/#{lcs}sga --cores #{threads}"
  puts sga_cmd
  `#{sga_cmd}`

  # make a histogram of the lengths of the contigs in the output assembly file
  #
  if File.exists?("#{path}/sga/#{lcs}sga.assemble-contigs.fa") # Fb_sga.assemble-contigs.fa
    abort "sga assembly file is empty!" if File.zero?("#{path}/sga/#{lcs}sga.assemble-contigs.fa")
    contigs = Bio::FastaFormat.open("#{path}/sga/#{lcs}sga.assemble-contigs.fa")
    histogram = {}
    contigs.each do |entry|
      bucket = (entry.seq.length*0.01).round
      if !histogram.has_key?(bucket)
        histogram[bucket]=0
      end
      histogram[bucket]+=1
    end
    File.open("#{required[:sga_output]}", "w") do |out|
      keys = histogram.keys.sort
      keys.each do |key|
        out.write "#{100*key}\t#{histogram[key]}\n"
      end
    end
  else
    abort "can't find #{path}/sga/#{lcs}sga.assemble-contigs.fa. sga must've failed"
  end
end

file required[:idba_output] => required[:khmered_reads] do
  puts "running idba trans on khmered reads..."

  if !Dir.exists?("#{path}/idba")
    mkdir_cmd = "mkdir #{path}/idba"
    `#{mkdir_cmd}`
  end

  # prepare reads
  if !File.exists?("#{path}/idba/#{lcs}.fx.fa")
    files = File.readlines("#{required[:khmered_reads]}").map{|n| n.chomp!}.compact
    fasta = File.open("#{path}/idba/#{lcs}.fx.fa", "w")
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

  idba_cmd = "#{idba} "
  idba_cmd << "-o #{path}/idba "
  idba_cmd << "-r #{path}/idba/#{lcs}.fx.fa "
  idba_cmd << "--num_threads #{threads} "           # number of threads
  idba_cmd << "--mink 21 "                          # minimum k value (<=124)
  idba_cmd << "--maxk 77 "                          # maximum k value (<=124)
  idba_cmd << "--step 4 "                           # increment of k-mer of each iteration
  idba_cmd << "--min_count 1 "                      # minimum multiplicity for filtering k-mer when building the graph
  idba_cmd << "--no_correct "                       # do not do correction
  idba_cmd << "--max_isoforms 6 "                   # maximum number of isoforms
  idba_cmd << "--similar 0.98"                      # similarity for alignment

  puts idba_cmd
  `#{idba_cmd}`

  # rm_cmd = "rm #{path}/#{lcs}.fx.fa"
  # `#{rm_cmd}`

  if File.exists?("#{path}/idba/contig.fa") # soapdtgraph.scafSeq
    abort "idba assembly file is empty!" if File.zero?("#{path}/idba/contig.fa")
    contigs = Bio::FastaFormat.open("#{path}/idba/contig.fa")
    histogram = {}
    contigs.each do |entry|
      bucket = (entry.seq.length*0.01).round
      if !histogram.has_key?(bucket)
        histogram[bucket]=0
      end
      histogram[bucket]+=1
    end
    File.open("#{required[:idba_output]}", "w") do |out|
      keys = histogram.keys.sort
      keys.each do |key|
        out.write "#{100*key}\t#{histogram[key]}\n"
      end
    end
  else
    abort "couldn't find idba output file contig"
  end

end

file required[:combined_assembly] do
  puts "combining outputs from soap, sga and idba and running cd-hit-est"

  soap = "#{path}/soap/#{lcs}soap.gapfill.fasta" 
  sga = "#{path}/sga/#{lcs}sga.assemble-contigs.fa"
  idba = "#{path}/idba/contig.fa"

  # concatenate all the output assemblies together
  cat_cmd = "cat #{soap} #{sga} #{idba} > #{path}/#{lcs}combined_contigs.fa"
  puts cat_cmd
  `#{cat_cmd}`

  # run cd-hit-est
  cd_cmd = "#{cd_hit_est} -i #{path}/#{lcs}combined_contigs.fa -o #{path}/#{lcs}cd_hit.fasta -T #{threads} -c 0.99 -M 5000"
  puts cd_cmd
  `#{cd_cmd}`

  output = "#{path}/#{lcs}cd_hit.fasta"
  File.open("#{required[:combined_assembly]}", "w") {|io| io.write output}
  # abort "stopping in combining assembly files"
end

file required[:bowtie_index] do
  puts "making bowtie2 index..."
  abort "Something went wrong with annotation. Can't find annotated fasta file" if !File.exists?("#{path}/#{lcs}annotated.fasta")
  index_cmd = "bowtie2-build #{path}/#{lcs}annotated.fasta #{path}/#{lcs}.index"
  puts index_cmd
  `#{index_cmd}`
  File.open("#{required[:bowtie_index]}", "w") {|out| out.write("#{index_cmd}\n")}
end

file required[:expression_output] => required[:annotation_output] do
  puts "running eXpress with trimmed reads against transcripts"

  # construct list of reads to align to transcripts
  left=Hash.new
  right=Hash.new
  single=Hash.new
  File.open("#{required[:trimmed_reads]}").each_line do |line|
    line.chomp!
    filename=File.basename(line)
    filepath=File.dirname(line)
    if filename=~/^t\..*(.)_(.)_R1\.fastq/ # $1 = replicate, $2 = section
      replicate = $1
      section = $2
      key = "#{section}-#{repliace}".to_sym
      left[key] = [] if !left.has_key?(key)
      left[key] << line
    elsif filename=~/^t\..*(.)_(.)_R2\.fastq/
      replicate = $1
      section = $2
      key = "#{section}-#{repliace}".to_sym
      right[key] = [] if !right.has_key?(key)
      right[key] << line
    elsif filename=~/^tU\..*(.)_(.)_R1\.fastq/
      replicate = $1
      section = $2
      key = "#{section}-#{repliace}".to_sym
      single[key] = [] if !single.has_key?(key)
      single[key] << line
    end
  end

  directories = []
  # run bowtie2 streamed into express
  # bowtie2 2.1.0 and express 1.5.1 were used to test this
  left.keys.each do |section|
    #make an output directory for this section
    mkdir_cmd = "mkdir #{path}/express_#{lcs}#{section}"
    puts mkdir_cmd
    directories << "#{path}/express_#{lcs}#{section}"
    if !File.exists?("#{path}/express_#{lcs}#{section}")
      `#{mkdir_cmd}`
    else
      puts "directory already exists"
    end

    express_cmd = "bowtie2 -t -a --very-sensitive -p #{threads} -x #{path}/#{lcs}.index "
    express_cmd << "-1 #{left[section].join(",")} "
    express_cmd << "-2 #{right[section].join(",")} "
    express_cmd << "-U #{single[section].join(",")} "
    express_cmd << " | express --output-align-prob -o #{path}/express_#{lcs}#{section} "
    express_cmd << " --no-update-check "   # extra batch rounds -B 2 can't do with streaming output
    express_cmd << " #{path}/#{lcs}annotated.fasta " # fasta file to align reads to
    puts express_cmd
    if !File.exists?("#{path}/express_#{lcs}#{section}/results.xprs")
      `#{express_cmd}`
    else
      puts "results.xprs already exists for section #{section}"
    end
  end

  count=0
  directories.each do |dir|
    File.open("#{dir}/results.xprs", "r").each_line do |line|
      line.chomp!
      if line.split(/\t+/)[14].to_f > 0
        count+=1
      end
    end
  end
  #
  # TODO 
  # make a summary of the expression counts of the annotated genes in one file
  # eg
  # AT1G01010.1   1.2e2   2.4e2   3.6e2   3.7e2   3.8e2   3.8e2
  # ...
  #
  File.open("#{required[:expression_output]}", "w") {|out| out.write "#{count} expressed transcripts"}
end

file required[:annotation_output] => required[:combined_assembly] do
  abort "ABORT: Something went wrong with #{required[:combined_assembly]} and the output file is empty!"  if File.size(required[:combined_assembly]) < 10
  puts "running RBUsearch to annotate transcripts..."
  if !File.exists?("#{path}/#{lcs}contigs.fasta")
    puts "removing short contigs"
    # remove contigs with length less than  X
    length_cutoff = 179
    count=0
    contigs = Bio::FastaFormat.open("#{path}/#{lcs}cd_hit.fasta")
    File.open("#{path}/#{lcs}contigs.fasta", "w") do |out|
      contigs.each do |entry|
        if entry.seq.length >= length_cutoff
          name = "contig_#{count}"
          out.write ">#{name}\n#{entry.seq}\n"
          count+=1
        end
      end
    end
  end
  
  cmd = "ruby rbusearch.rb --query #{path}/#{lcs}contigs.fasta --target #{protein_reference} --output #{path}/rbusearch --cores #{threads} --prefix #{lcs} --verbose"
  puts cmd
  log = `#{cmd}`
  File.open("#{path}/rbusearch.log", "w") {|out| out.write log}

  cmd = "mv #{path}/rbusearch/#{lcs}annotated.fasta #{path}/#{lcs}annotated.fasta"
  `#{cmd}`

  cmd = "cp #{path}/rbusearch/reciprocal_hits.txt #{required[:annotation_output]}"
  `#{cmd}`
end

task :default => :build

# task :build => [:expression]
task :build => [:annotation]

task :index => [required[:bowtie_index]]

task :expression => [:annotation, :index, required[:expression_output]]

task :annotation => [:assemble, required[:annotation_output]]

# task :scaffold => [:assemble, required[:scaffold_output]] # TODO [ ]

task :assemble => [:khmer, :soap, :sga, :idba, :combine]

task :combine => [required[:combined_assembly]]

task :idba => [required[:idba_output]]

task :sga => [required[:sga_output]]

task :soap => [:config, required[:soap_output]]

task :config => [required[:config]]

task :khmer => [:correct, required[:khmered_reads]]

task :yaml => [required[:yaml]]

task :correct => [:trim, :yaml, required[:corrected_reads]]

task :trim => [:fastqc, required[:trimmed_reads]]

task :fastqc => [:input, required[:quality]]

task :clean do
  files = required.values.delete_if {|a| a=="#{required[:input_reads]}"} # don't delete the input files :P
  files.each do |file|
    if File.exists?(file)
      `rm #{file}`
    end
  end
  `rm dataset*yaml`
  `rm datasets`
end

task :test do
   
end
