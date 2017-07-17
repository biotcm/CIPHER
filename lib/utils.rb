#!/usr/bin/env ruby
require 'biotcm'
require 'ruby-progressbar'

# Get a list of OMIM phenotypes
def omim_phenotypes
  return @omim_phenotypes if @omim_phenotypes

  # Load mim2gene
  if File.exist?('tmp/mim2gene.txt')
    mim2gene = File.read('tmp/mim2gene.txt')
  else
    mim2gene = BioTCM.curl('https://omim.org/static/omim/data/mim2gene.txt')
    File.open('tmp/mim2gene.txt', 'w').puts mim2gene
  end

  # Select MIMs
  mims = mim2gene.split("\n").map do |line|
    col = line.split("\t")

    if col[0] =~ /^#/
      nil # throw comments
    elsif %w[gene moved/removed].include?(col[1])
      nil # throw non-phenotypes
    else
      col[0]
    end
  end

  @omim_phenotypes = mims.compact
end
