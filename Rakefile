task :default => :build

file "trimmed_reads.fastq" => "input_reads.fastq" do
  puts "creating trimmed reads..."
  sh "touch trimmed_reads.fastq"
end

file "quality_check_output" => "input_reads.fastq" do
  puts "running quality check..."
  sh "touch quality_check_output"
end

file "expression_quantification_output.csv" => "assembly.fasta" do
  puts "running eXpress with trimmed reads against transcripts"
  sh "touch expression_quantification_output.csv"
end

file "assembly.fasta" => "khmered_reads.fastq" do
  puts "running assemblotron on khmered reads..."
  sh "touch assembly.fasta"
end

file "khmered_reads.fastq" => "trimmed_reads.fastq" do
  puts "running khmer to reduce coverage of reads..."
  sh "touch khmered_reads.fastq"
end

file "annotation_summary.csv" => "assembly.fasta" do
  puts "running RBUsearch and round-robin to annotate transcripts..."
  sh "touch annotation_summary.csv"
end

task :build => [:expression, :annotation]

task :fastqc => "quality_check_output"

task :trim => [:fastqc, "trimmed_reads.fastq"]

task :khmer => [:trim, "khmered_reads.fastq"]

task :assemble => [:khmer, "assembly.fasta"]

task :annotation => [:assemble, "annotation_summary.csv"]

task :expression => [:assemble, "expression_quantification_output.csv"]

task :clean do
  files = ["trimmed_reads.fastq", 
    "quality_check_output", 
    "expression_quantification_output.csv",
    "annotation_summary.csv",
    "assembly.fasta",
    "khmered_reads.fastq",
    "trimmed_reads.fastq"
  ]
  files.each do |file|
    if File.exists?(file)
      sh "rm #{file}"
    end
  end
end

