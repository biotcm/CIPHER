#!/usr/bin/env ruby
# This script aims at parsing CIPHER's logs

File.open(ARGV[0][0...-3] + 'txt', 'w') do |fout|
  File.open(ARGV[0]).each do |line|
    next unless line.index("\t")
    col = line.chomp.split("\t")

    case col[1]
    when /info/
    when /data/
      fout.puts col[2..-1].join("\t")
    end
  end
end
