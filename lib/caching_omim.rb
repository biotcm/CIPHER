#!/usr/bin/env ruby
require_relative './utils'

# Cache MIMs
pb = ProgressBar.create(total: omim_phenotypes.size, format: '%t: |%B| %a %e')
omim_phenotypes.each do |mim|
  pb.title = "Caching MIM \##{mim}"
  pb.increment
  BioTCM::Databases::OMIM.get(mim) rescue retry
end
