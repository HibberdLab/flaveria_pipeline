required = {
  :assembly_output => "assembly.fasta",
  :quality => "quality_check_output",
  :input_reads => "raw_data",
  :trimmed_reads => "trimmed_reads",
  :corrected_reads => "corrected_reads",
  :yaml => "dataset.yaml",
  :khmered_reads => "khmered_reads",
  :annotation_output => "annotation_summary.csv",
  :expression_output => "expression_quantification_output.csv"
}

threads = 1
fastqc_path = "/applications/fastqc_v0.10.1/FastQC/fastqc"

task :input do
  File.open(required[:input_reads], "r").each_line do |line|
    line.chomp!
    if !File.exists?("#{line}")
      abort "Can't find #{line}"
    end
  end
end

file required[:quality] => required[:input_reads] do
  if !Dir.exists?("fastqc_output")
    sh "mkdir fastqc_output"
  end
  files = ""
  Dir.chdir("fastqc_output") do
    File.open("../#{required[:input_reads]}", "r").each_line do |line|
      line.chomp!
      fastqc_line = line.split(".")[0..-2].join(".") + "_fastqc"
      if File.exists?("#{fastqc_line}")
        puts "Found #{fastqc_line}"
      else
        puts "Didn't find #{fastqc_line}"
        files << " #{line} "
      end
    end
  end

  # maybe add a thing here that checks the contents of "fastqc_output" and sees if the fastqc has already been done?
  sh "#{fastqc_path} --kmers 5 --threads #{threads} --outdir fastqc_output #{files}" if files.length>0
  
  fastqc_summary = Hash.new
  Dir.chdir("fastqc_output") do
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
    if File.exists?("t.#{line}")
      list_of_trimmed_reads << "t.#{line}\n"
    end
    if File.exists?("t.#{line}U")
      rename =  "mv t.#{line}U tU.#{line}"
      puts rename
      sh "#{rename}"
      list_of_trimmed_reads << "tU.#{line}\n"
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
    if line=~/^t\..*R1.*/
      hash[:left] << line
    elsif line=~/^t\..*R2.*/
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
  
  cmd = "python ~/apps/SPAdes-2.5.1-Linux/bin/spades.py --dataset #{required[:yaml]} --only-error-correction --disable-gzip-output -m 90 -t #{threads} -o output.spades"
  puts cmd
  hammer_log = `#{cmd}`
  File.open("hammer.log", "w") {|out| out.write hammer_log}
  paired = []
  single = []
  Dir.chdir("output.spades") do
    Dir.chdir("corrected") do
      Dir["*fastq"].each do |fastq|
        if fastq =~ /t\..*R[12].*fastq/
          #paired
          paired << fastq
        elsif fastq =~ /tU.*fastq/
          #single
          single << fastq
        end
      end
    end
  end
  File.open("#{required[:corrected_reads]}", "w") do |out|
    paired.each do |pe|
      out.write "output.spades/corrected/#{pe}\n"
    end
  end
  File.open("single_#{required[:corrected_reads]}", "w") do |out|
    single.each do |sng|
      out.write "output.spades/corrected/#{sng}\n"
    end
  end
end

file required[:khmered_reads] => required[:corrected_reads] do
  puts "running khmer to reduce coverage of reads..."
  khmer_cmd = "ruby khmer-batch.rb "
  khmer_cmd << "--input #{required[:corrected_reads]} "
  khmer_cmd << "--paired "
  khmer_cmd << "--interleave "
  khmer_cmd << "--verbose "
  puts khmer_cmd
  `#{khmer_cmd}`

  #sh "touch #{required[:khmered_reads]}"
end

file required[:assembly_output] => required[:khmered_reads] do
  puts "running assemblotron on khmered reads..."
  sh "touch #{required[:assembly_output]}"
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

task :assemble => [:khmer, required[:assembly_output]]

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

