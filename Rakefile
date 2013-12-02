require 'bio'

required = {
  :assembly_output => "assembly_stats",
  :quality => "quality_check_output",
  :input_reads => "raw_data",
  :trimmed_reads => "trimmed_reads",
  :corrected_reads => "corrected_reads",
  :hammer_log => "hammer.log",
  :single_corrected_reads => "single_corrected_reads",
  :yaml => "datasets",
  :config => "soapdt.config",
  :khmered_reads => "khmered_reads",
  :annotation_output => "annotation_summary.csv",
  :expression_output => "expression_output",
  :bowtie_index => "bowtie_index"
}

threads = 22
memory = 90
maximum_files_to_hammer_at_a_time = 6
fastqc_path = "/applications/fastqc_v0.10.1/FastQC/fastqc"
hammer_path = "~/apps/SPAdes-2.5.1-Linux/bin/spades.py"
trimmomatic_path = "/home/cmb211/apps/Trimmomatic-0.32/trimmomatic-0.32.jar"
khmer_path = "/home/cmb211/.local/bin/normalize-by-median.py"
lcs = ""
path=""

task :input do
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
  s = a.min_by(&:size)
  lcs = catch(:hit) {  s.size.downto(1) { |i| (0..(s.size - i)).each { |l| throw :hit, s[l, i] if a.all? { |item| item.include?(s[l, i]) } } } }
  # lcs=lcs[0..-2] if lcs[-2]=~/\_\.\,\-/
  # lcs = File.basename(lcs)
  puts "lcs = #{lcs} path=#{path}"
end

file required[:quality] => required[:input_reads] do
  puts "running fastqc reads..."
  if !Dir.exists?("#{path}/fastqc_output")
    sh "mkdir #{path}/fastqc_output"
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

  sh "#{fastqc_path} --kmers 5 --threads #{threads} --outdir #{path}/fastqc_output #{files}" if files.length>0
  
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

  trim_batch_cmd = "ruby trim-batch.rb "
  trim_batch_cmd << "--jar #{trimmomatic_path} "
  trim_batch_cmd << "--pairedfile #{required[:input_reads]} "
  #trim_batch_cmd << "--singlefile #{required[:single_input_reads]} "
  trim_batch_cmd << "--threads #{threads} "
  trim_batch_cmd << "--quality 15 "
  # puts trim_batch_cmd
  sh "#{trim_batch_cmd}"
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
      sh "#{rename}"
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
      yaml << "," if j < right_.length-1
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
  File.open("#{required[:yaml]}", "r").each_line do |dataset_line|
    puts "running hammer on #{dataset_line}"
    dataset_line.chomp!
    cmd = "python #{hammer_path} --dataset #{dataset_line} --only-error-correction --disable-gzip-output -m #{memory} -t #{threads} -o #{path}/output_#{count}.spades"
    output_directories << "#{path}/output_#{count}.spades"
    if !File.exists("output_#{count}.spades")
      puts cmd
      hammer_log = `#{cmd}`
      File.open("#{path}/hammer_#{dataset_line}.log", "w") {|out| out.write hammer_log}
    end
    count+=1
  end

  paired = []
  single = []

  output_directories.each do |dir|
    Dir.chdir(dir) do
      Dir.chdir("corrected") do
        Dir["*fastq"].each do |fastq|
          if fastq =~ /t\..*R[12].*fastq/
            paired << "#{path}/#{dir}/corrected/#{fastq}"
          elsif fastq =~ /tU.*fastq/
            single << "#{path}/#{dir}/corrected/#{fastq}"
          end
        end
      end
    end
  end
  if paired.length == 0 and single.length == 0
    abort "Something went wrong with BayesHammer and no corrected reads were created"
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

  # run khmer-batch on paired corrected reads
  khmer_cmd = "ruby khmer-batch.rb "
  khmer_cmd << "--script #{khmer_path} "
  khmer_cmd << "--input #{required[:corrected_reads]} "
  khmer_cmd << "--paired "
  khmer_cmd << "--interleave "
  khmer_cmd << "--verbose "
  puts khmer_cmd
  `#{khmer_cmd}` 
  
  # cat all the single corrected reads into a single fastq file to run khmer on them
  cat_cmd = "cat "
  corrected_path=""
  File.open("single_#{required[:corrected_reads]}", "r").each_line do |line|
    line.chomp!
    corrected_path = File.dirname(line)
    cat_cmd << " #{line} "
  end
  cat_cmd << " > #{path}/#{lcs}.single.fastq"
  puts cat_cmd
  `#{cat_cmd}`

  # cat all the .keep files together from the paired khmer run
  cat_cmd = "cat "
  Dir.chdir(corrected_path) do |dir|
    Dir["*keep"].each do |file|
      cat_cmd << " #{corrected_path}/#{file} "
    end
  end
  cat_cmd << " > #{path}/#{lcs}.khmered.fastq"
  puts cat_cmd
  `#{cat_cmd}`

  # run khmer on the catted single corrected fastq files
  khmer_cmd = "ruby khmer-batch.rb "
  khmer_cmd << "--files #{path}/#{lcs}.single.fastq "
  khmer_cmd << "--verbose "
  puts khmer_cmd
  `#{khmer_cmd}`

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

file required[:assembly_output] => required[:khmered_reads] do
  puts "running soap on khmered reads..."

  soap_cmd = "SOAPdenovo-Trans-127mer all -s #{required[:config]} -o #{path}/#{lcs}soap -p #{threads}"
  sh "#{soap_cmd}"
  
  if File.exists?("#{path}/#{lcs}soap.contig") # soapdtgraph.scafSeq
    # calculate n50 and stats etc
    File.open("#{required[:assembly_output]}", "w") do |out|
      contigs = Bio::FastaFormat.open("#{path}/#{lcs}soap.contig")
      contigs.each do |entry|
        out.write "#{entry.definition} #{entry.seq.length}\n"
      end
    end
  else
    abort "can't find #{lcs}soap.contig. soap must've failed or something"
  end
end

file required[:bowtie_index] do
  puts "making bowtie2 index..."
  # construct a bowtie2 index  
  index_cmd = "bowtie2-build #{path}/#{lcs}soap.contig #{path}/#{lcs}.index"
  puts index_cmd
  sh "#{index_cmd}"
  File.open("#{required[:bowtie_index]}", "w") {|out| out.write("#{index_cmd}\n")}
end

file required[:expression_output] => required[:assembly_output] do
  puts "running eXpress with trimmed reads against transcripts"

  # construct list of reads to align to transcripts
  left=[]
  right=[]
  single=[] # can you add single reads as well as paired reads to bowtie2? YES
  File.open("#{required[:trimmed_reads]}").each_line do |line|
    line.chomp!
    filename=File.basename(line)
    filepath=File.dirname(line)
    if filename=~/^t\..*R1.*fastq/
      left << line
    elsif filename=~/^t\..*R2.*fastq/
      right << line
    elsif filename=~/^tU\..*fastq/
      single << line
    end
  end

  # run bowtie2 streamed into express
  # bowtie2 2.1.0 and express 1.5.0 were used to test this
  express_cmd = "bowtie2 -t -a --very-sensitive -p #{threads} -x #{path}/#{lcs}.index "
  express_cmd << "-1 #{left.join(",")} "
  express_cmd << "-2 #{right.join(",")} "
  express_cmd << "-U #{single.join(",")} "
  express_cmd << " | express --output-align-prob -o #{path} #{path}/#{lcs}soap.contig "
  puts express_cmd
  sh "#{express_cmd}"

  count=0
  File.open("#{path}/results.xprs", "r").each_line do |line|
    line.chomp!
    if line.split(/\t+/)[14].to_f > 0
      count+=1
    end
  end
  File.open("#{required[:expression_output]}", "w") {|out| out.write "#{count} expressed transcripts"}
end

file required[:annotation_output] => required[:assembly_output] do
  if File.size(required[:assembly_output]) > 0
    puts "running RBUsearch and round-robin to annotate transcripts..."
    sh "touch #{required[:annotation_output]}"
  else
    abort "ABORT: Something went wrong with #{required[:assembly_output]} and the output file is empty!"
  end
end

task :default => :build

task :build => [:expression, :annotation]

task :index => [required[:bowtie_index]]

task :expression => [:assemble, :index, required[:expression_output]]

task :annotation => [:assemble, required[:annotation_output]]

task :assemble => [:khmer, :config, required[:assembly_output]]

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
      sh "rm #{file}"
    end
  end
end

