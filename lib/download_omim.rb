#!/usr/bin/env ruby
require 'biotcm'
require 'ruby-progressbar'

# Load mim2gene
if File.exist?('../tmp/mim2gene.txt')
  mim2gene = File.read('../tmp/mim2gene.txt')
else
  mim2gene = BioTCM.curl('https://omim.org/static/omim/data/mim2gene.txt')
  File.open('../tmp/mim2gene.txt', 'w').puts mim2gene
end

# Select MIMs
mims = mim2gene.split("\n").map do |line|
  col = line.split("\t")

  if col[0] =~ /^#/
    nil
  elsif %w[gene moved/removed].include?(col[1])
    nil
  else
    col[0]
  end
end
mims.compact!

# Cache MIMs locally
pb = ProgressBar.create(total: mims.size, format: '%t: |%B| %a %e')
mims.each do |mim|
  pb.title = "Caching MIM \##{mim}"
  pb.increment
  BioTCM::Databases::OMIM.get(mim) rescue retry
end
