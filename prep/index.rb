#!/usr/bin/env ruby
require_relative './utils'

gene = {}
File.open('../temp/ppi.txt').each do |line|
  line.chomp.split("\t").each { |g| gene[g] ||= gene.size + 1 }
end

mesh = {
  'A' => { index: 1, mindex: 'A', term: 'A', parent: nil },
  'C' => { index: 2, mindex: 'C', term: 'C', parent: nil }
}
BioTCM::Table.load('../temp/mtrees2017.tab').each_row do |rk, row|
  mesh[rk] = {
    index: mesh.size + 1,
    mindex: rk,
    term: row['term'],
    parent: row['parent']
  }
end

pheno = omim_phenotypes.zip(1..omim_phenotypes.size).to_h

File.open('../temp/inter_gene.txt', 'w').puts(gene.map { |k, v| "#{v}\t#{k}" })
File.open('../temp/inter_mesh.txt', 'w').puts(mesh.map { |_, h| h.values.join("\t") })
File.open('../temp/inter_pheno.txt', 'w').puts(pheno.map { |k, v| "#{v}\t#{k}" })

File.open('../temp/inner_ppi.txt', 'w') do |fout|
  File.open('../temp/ppi.txt').each do |line|
    fout.puts line.chomp.split("\t").map { |g| gene[g] }.join("\t")
  end
end

%w[direct extend curation].each do |scope|
  File.open("../temp/inner_pheno2gene_#{scope}.txt", 'w') do |fout|
    File.open("../temp/pheno2gene_#{scope}.txt").each do |line|
      col = line.chomp.split("\t").map.with_index do |c, i|
        i.zero? ? pheno[c] : gene[c]
      end.compact
      fout.puts col.join("\t") if col.size > 1
    end
  end
end

File.open('../temp/inner_pheno2mesh_freq.txt', 'w') do |fout|
  File.open('../temp/pheno2mesh_freq.txt').each do |line|
    col = line.chomp.split("\t")
    fout.puts [pheno[col[0]], mesh[col[1]][:index], col[2]].join("\t")
  end
end
