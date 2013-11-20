required_files = {
  :assembly_output => "assembly.fasta",
  :quality => "quality_check_output",
  :input_reads => "input_reads.fastq",
  :trimmed_reads => "trimmed_reads.fastq",
  :khmered_reads => "khmered_reads.fastq",
  :annotation_output => "annotation_summary.csv",
  :expression_output => "expression_quantification_output.csv"
}

file required_files[:quality] => required_files[:input_reads] do
  puts "running quality check..."
  sh "touch #{required_files[:quality]}"
end

file required_files[:trimmed_reads] => required_files[:input_reads] do
  puts "creating trimmed reads..."
  sh "touch #{required_files[:trimmed_reads]}"
end

file required_files[:assembly_output] => required_files[:khmered_reads] do
  puts "running assemblotron on khmered reads..."
  sh "touch #{required_files[:assembly_output]}"
end

file required_files[:khmered_reads] => required_files[:trimmed_reads] do
  puts "running khmer to reduce coverage of reads..."
  sh "touch #{required_files[:khmered_reads]}"
end

file required_files[:expression_output] => required_files[:assembly_output] do
  puts "running eXpress with trimmed reads against transcripts"
  sh "touch #{required_files[:expression_output]}"
end

file required_files[:annotation_output] => required_files[:assembly_output] do
  if File.size(required_files[:assembly_output]) > 0
    puts "running RBUsearch and round-robin to annotate transcripts..."
    sh "touch #{required_files[:annotation_output]}"
  else
    abort "ABORT: Something went wrong with #{required_files[:assembly_output]} and the output file is empty!"
  end
end

task :default => :build

task :build => [:expression, :annotation]

task :fastqc => required_files[:quality]

task :trim => [:fastqc, required_files[:trimmed_reads]]

task :khmer => [:trim, required_files[:khmered_reads]]

task :assemble => [:khmer, required_files[:assembly_output]]

task :annotation => [:assemble, required_files[:annotation_output]]

task :expression => [:assemble, required_files[:expression_output]]

task :clean do
  files = required_files.values.delete_if {|a| a=="input_reads.fastq"} # don't delete the input files :P
  files.each do |file|
    if File.exists?(file)
      sh "rm #{file}"
    end
  end
end

