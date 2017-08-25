#!/usr/bin/env ruby
#
# iCIPHER workflow - prepare all input data files
#
# This script describes the workflow of iCIPHER. It may take days to finish this
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
