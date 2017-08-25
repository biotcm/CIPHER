#!/usr/bin/env ruby
#
# CIPHER2 workflow - prepare all input data files
#
# This script describes the workflow of CIPHER2. It may take days to finish this
# pipeline. So it is recommended to run each step separately.
#
require 'logger'
logger = Logger.new(STDOUT)

# Setup the environment
if `gem list --local --installed ^bundler$` == 'false'
  logger.info('Installing Bundler...')
  system('gem install bundler')
end
logger.info('Checking dependencies...')
system('bundle install')

# Process OMIM content
logger.info('Caching OMIM content...')
system('ruby caching_omim.rb')
logger.info('Extracting OMIM content...')
system('ruby extract_omim.rb')

# Prepare the PPIs
logger.info('Preparing PPIs...')
system('ruby prepare_ppi.rb')

# Index all inputs
logger.info('Indexing all inputs...')
system('ruby index.rb')
