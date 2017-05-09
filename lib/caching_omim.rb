#!/usr/bin/env ruby
require_relative './utils'

# Get the list
mims = get_OMIM_phenotypes

# Cache MIMs
pb = ProgressBar.create(total: mims.size, format: '%t: |%B| %a %e')
mims.each do |mim|
  pb.title = "Caching MIM \##{mim}"
  pb.increment
  BioTCM::Databases::OMIM.get(mim) rescue retry
end
