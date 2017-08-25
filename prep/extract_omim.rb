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
def extract_phenotype_mesh_relationships(entry)
  meshs = {}
  content = ''

  entry['textSectionList'].each do |section|
    content += section['textSection']['textSectionContent']
  end

  if entry['clinicalSynopsisExists']
    if entry['clinicalSynopsis']['oldFormatExists']
      entry['clinicalSynopsis']['oldFormat'].values.each { |v| content += v }
    else
      entry['clinicalSynopsis'].each do |key, value|
        next unless entry['clinicalSynopsis']["#{key}Exists"] == true
        content += value
      end
    end
  end

  # Ignore letter cases
  content = content.downcase

  mesh_tree_table.col('term').each do |mindex, term|
    counts = term.split(/,\s*/).map do |segment|
      content.scan(Regexp.new("\\W#{segment}(s|es)?\\W")).size
    end
    meshs[mindex] = counts.max if counts.max > 0
  end

  meshs.tap { |m| yield(m) }
end

fout = {
  pheno2gene: {
    direct: File.open('../temp/pheno2gene_direct.txt', 'w'),
    extend: File.open('../temp/pheno2gene_extend.txt', 'w')
  },
  pheno2mesh: File.open('../temp/pheno2mesh_freq.txt', 'w')
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

  extract_phenotype_mesh_relationships(entry) do |meshs|
    meshs.each do |mindex, count|
      fout[:pheno2mesh].puts [mim, mindex, count].join("\t")
    end
  end
end
