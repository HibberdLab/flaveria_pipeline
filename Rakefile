require 'bio'

required = {
  :assembly_output => "assembly_stats",
  :quality => "quality_check_output",
  :input_reads => "raw_data",
  :trimmed_reads => "trimmed_reads",
  :corrected_reads => "corrected_reads",
  :hammer_log => "hammer.log",
  :single_corrected_reads => "single_corrected_reads",
  :yaml => "dataset.yaml",
  :config => "soapdt.config",
  :khmered_reads => "khmered_reads",
  :annotation_output => "annotation_summary.csv",
  :expression_output => "expression_quantification_output.csv"
}

threads = 12
fastqc_path = "/applications/fastqc_v0.10.1/FastQC/fastqc"
lcs = ""
path=""

task :input do
  a=[]
  File.open(required[:input_reads], "r").each_line do |line|
    line.chomp!
    if File.exists?("#{line}")
      a << line
      path = File.dirname(line)
    else
      abort "Can't find #{line}"
    end
  end
  a.map {|f| f.gsub!(".fastq","")}
  s = a.min_by(&:size)
  lcs = catch(:hit) {  s.size.downto(1) { |i| (0..(s.size - i)).each { |l| throw :hit, s[l, i] if a.all? { |item| item.include?(s[l, i]) } } } }
  # lcs=lcs[0..-2] if lcs[-2]=~/\_\.\,\-/
  lcs = File.basename(lcs)
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
        puts "Found #{fastqc_line}"
      else
        puts "Didn't find #{fastqc_line}"
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
  File.open("#{required[:quality]}", "w") do |out|
    fastqc_summary.each_pair do |key, value|
      out.write "#{key}\t#{value}\n"
    end
  end
  
end

file required[:trimmed_reads] => required[:input_reads] do
  puts "creating trimmed reads..."

  trim_batch_cmd = "ruby trim-batch.rb "
  trim_batch_cmd += "--jar /home/cmb211/apps/Trimmomatic-0.32/trimmomatic-0.32.jar "
  trim_batch_cmd += "--pairedfile #{required[:input_reads]} "
  #trim_batch_cmd += "--singlefile #{required[:single_input_reads]} "
  trim_batch_cmd += "--threads #{threads} "
  trim_batch_cmd += "--quality 15 "
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
      rename =  "mv #{path}/t.#{filename}U #{path}/tU.#{filename}"
      puts rename
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
  yaml = "[\n"
  yaml << "  {\n"
  yaml << "    orientation: \"fr\",\n"
  yaml << "    type: \"paired-end\",\n"
  yaml << "    left reads: [\n"
  hash[:left].each_with_index do |left_read, i|
    yaml << "      \"#{left_read}\""
    yaml << "," if i < hash[:left].length-1
    yaml << "\n"
  end
  yaml << "    ],\n"
  yaml << "    right reads: [\n"
  hash[:right].each_with_index do |right_read, i|
    yaml << "      \"#{right_read}\""
    yaml << "," if i < hash[:right].length-1
    yaml << "\n"
  end
  yaml << "    ],\n"
  yaml << "  },\n"
  yaml << "  {\n"
  yaml << "    type: \"single\",\n"
  yaml << "    single reads: [\n"
  hash[:single].each_with_index do |single_read, i|
    yaml << "      \"#{single_read}\""
    yaml << "," if i < hash[:single].length-1
    yaml << "\n"
  end
  yaml << "    ]\n"
  yaml << "  }\n"
  yaml << "]\n"
  File.open("#{required[:yaml]}", "w") do |out|
    out.write yaml
  end

end

file required[:corrected_reads] => required[:trimmed_reads] do
  puts "running bayeshammer to correct reads..."
  
  cmd = "python ~/apps/SPAdes-2.5.1-Linux/bin/spades.py --dataset #{required[:yaml]} --only-error-correction --disable-gzip-output -m 90 -t #{threads} -o #{path}/output.spades"
  puts cmd
  hammer_log = `#{cmd}`
  File.open("hammer.log", "w") {|out| out.write hammer_log}
  paired = []
  single = []
  Dir.chdir("#{path}/output.spades") do
    Dir.chdir("corrected") do
      Dir["*fastq"].each do |fastq|
        if fastq =~ /t\..*R[12].*fastq/
          paired << fastq
        elsif fastq =~ /tU.*fastq/
          single << fastq
        end
      end
    end
  end
  paired.sort!
  File.open("#{required[:corrected_reads]}", "w") do |out|
    paired.each do |pe|
      out.write "#{path}/output.spades/corrected/#{pe}\n"
    end
  end
  File.open("single_#{required[:corrected_reads]}", "w") do |out|
    single.each do |sng|
      out.write "#{path}/output.spades/corrected/#{sng}\n"
    end
  end
end

file required[:khmered_reads] => required[:corrected_reads] do
  puts "running khmer to reduce coverage of reads..."

  # puts "path = #{path}"
  # puts "lcs = #{lcs}"

  # run khmer-batch on paired corrected reads
  khmer_cmd = "ruby khmer-batch.rb "
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
  reads = File.open("#{required[:khmered_reads]}", "r")
  left = reads.readline.chomp!
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

file required[:expression_output] => required[:assembly_output] do
  puts "running eXpress with trimmed reads against transcripts"
  sh "touch #{required[:expression_output]}"
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

task :expression => [:assemble, required[:expression_output]]

task :annotation => [:assemble, required[:annotation_output]]

task :assemble => [:khmer, :config, required[:assembly_output]]

task :config => [required[:config]]

task :khmer => [:correct, required[:khmered_reads]]

task :yaml => [required[:yaml]]

task :correct => [:trim, :yaml, required[:corrected_reads]]

task :trim => [:fastqc, required[:trimmed_reads]]

task :fastqc => [:input, required[:quality]]

task :clean do
  files = required.values.delete_if {|a| a=="raw_data"} # don't delete the input files :P
  files.each do |file|
    if File.exists?(file)
      sh "rm #{file}"
    end
  end
end

