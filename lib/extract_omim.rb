#!/usr/bin/env ruby
require_relative './utils'
BioTCM.logger.level = Logger::ERROR
BioTCM::Databases::HGNC.ensure
GD = BioTCM::Apps::GeneDetector.new
PB = ProgressBar.create(total: omim_phenotypes.size, format: '%t: |%B| %a %e')

# Extract phenotype-gene relationships from an OMIM entry
def extract_phenotype_gene_relationships(entry)
  genes = { direct: [], extend: [] }

  if entry['phenotypeMapExists']
    genes[:direct] = entry['phenotypeMapList']
      .map { |item| item['phenotypeMap']['geneSymbols'].split(', ') }
      .flatten.to_formal_symbols.uniq.reject(&:empty?)
  end

  lines = (entry['textSectionList'] || [])
    .map { |ts| ts['textSection']['textSectionContent'] }.join("\n")
  genes[:extend] = (genes[:direct] + GD.detect(lines)).uniq

  genes.tap { |g| yield(g) }
end

# Extract phenotype-MeSH relationships from an OMIM entry
def extract_phenotype_mesh_relationships(_entry)
  raise NotImplementedError # TODO
end

fout = {
  pheno2gene: {
    direct: File.open('tmp/pheno2gene_direct.txt', 'w'),
    extend: File.open('tmp/pheno2gene_extend.txt', 'w')
  }
}

omim_phenotypes.each do |mim|
  PB.increment
  PB.title = "Extracting MIM \##{mim}"
  entry = BioTCM::Databases::OMIM.get(mim)

  extract_phenotype_gene_relationships(entry) do |genes|
    %i[direct extend].each do |scope|
      fout[:pheno2gene][scope].puts genes[scope].unshift(mim).join("\t")
    end
  end

  # extract_phenotype_mesh_relationships(entry)
end
