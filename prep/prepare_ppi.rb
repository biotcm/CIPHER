#!/usr/bin/env ruby
require_relative './utils'
BioTCM.logger.level = Logger::ERROR
BioTCM::Databases::HGNC.ensure

# Change this value to change the source of PPIs
SOURCE = 'old_extended'.freeze

File.open('../temp/ppi.txt', 'w') do |fout|
  case SOURCE
  when /^(old_extended|old_hprd)$/
    dict = {}

    File.open("../data/#{SOURCE}/inter_gene.txt").each do |line|
      col = line.chomp.split("\t", -1)
      sym = [
        col[4].to_formal_symbol,
        col[3].refseq2symbol,
        col[2].uniprot2symbol
      ].reject(&:empty?).uniq.first
      dict[col[0]] = sym if sym
    end

    edges = File.open("../data/#{SOURCE}/inner_ppi.txt").map do |line|
      col = line.chomp.split("\t").map { |i| dict[i] }.compact
      if col.first == col.last
        nil
      elsif col.size < 2
        nil
      else
        col.sort.join("\t")
      end
    end

    fout.puts edges.compact.sort
  end
end
