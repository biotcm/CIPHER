#!/usr/bin/env ruby
#
# iCIPHER workflow
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
unless File.exist?('temp/omim_content.bin')
  logger.info('Caching OMIM content...')
  system('ruby lib/caching_omim.rb')
  logger.info('Extracting OMIM content...')
  system('ruby lib/extract_omim.rb')
end
