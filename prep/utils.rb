#!/usr/bin/env ruby
require 'yaml'
require 'biotcm'
require 'ruby-progressbar'

# Load OMIM mim2gene
def omim_mim2gene
  if File.exist?('../temp/mim2gene.txt')
    @omim_mim2gene ||= File.read('../temp/mim2gene.txt')
  else
    @omim_mim2gene = BioTCM.curl('https://omim.org/static/omim/data/mim2gene.txt')
    File.open('../temp/mim2gene.txt', 'w').puts @omim_mim2gene
  end

  @omim_mim2gene
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

def mesh_tree_table
  if File.exist?('../temp/mtrees2017.tab')
    @mesh_tree_table ||= BioTCM::Table.load('../temp/mtrees2017.tab')
  else
    @mesh_tree_table = BioTCM::Table.new(primary_key: 'mindex', col_keys: %w[term parent])

    File.open('../data/mtrees2017.bin').each do |line|
      term, mindex = line.chomp.split(';')
      next unless mindex =~ /^[AC]/

      @mesh_tree_table.row(
        mindex,
        'term' => term.downcase,
        'parent' => mindex =~ /\./ ? mindex.split('.')[0...-1].join('.') : mindex[0]
      )
    end

    @mesh_tree_table.save('../temp/mtrees2017.tab')
  end

  @mesh_tree_table
end
