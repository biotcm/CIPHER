#!/usr/bin/env ruby
require 'yaml'
require 'biotcm'
require 'ruby-progressbar'

# Load OMIM mim2gene
def omim_mim2gene
  if File.exist?('tmp/mim2gene.txt')
    mim2gene = File.read('tmp/mim2gene.txt')
  else
    mim2gene = BioTCM.curl('https://omim.org/static/omim/data/mim2gene.txt')
    File.open('tmp/mim2gene.txt', 'w').puts mim2gene
  end

  mim2gene
end

# Load OMIM phenotypes
def omim_phenotypes
  @omim_phenotypes ||= omim_mim2gene.split("\n").map do |line|
    col = line.split("\t")

    if col[0] =~ /^#/
      nil # throw comments
    elsif %w[gene moved/removed].include?(col[1])
      nil # throw non-phenotypes
    else
      col[0]
    end
  end.compact
end
