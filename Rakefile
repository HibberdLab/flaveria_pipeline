task :default => [:build]

task :fastqc  do
  puts "Running quality checking..."
end

task :trim => [:fastqc]do
  puts "Trimming reads..."
end

task :khmer => [:trim] do
  puts "Reducing coverage..."
end

task :assemble  => [:khmer] do
  puts "Assembling reads..."
end

task :annotation => [:assemble] do
  puts "Annotating transcripts..."
end

task :expression => [:annotation] do
  puts "Calculating expression quantification..."
end

task :build => [:expression] do
  puts "Running pipeline..."
end