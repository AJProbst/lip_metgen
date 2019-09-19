#!/usr/bin/ruby -w

# The MIT License (MIT)
# Copyright (c) 2019 Alexander J Probst

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# input files
spc2lip = ("species2lipids_best_hit.txt")
lip2lip = ("first_second_lipid_connections.txt")
puts ("reading input data ...")

# read primary connections for each species
total_prim_lipids =[]
species = Hash.new{|h,k| h[k]=[]}
File.open(spc2lip).each do |line|
    line.chomp
    cols = line.split(/\t/)
    s = cols[0]
    l = cols[1]
    p = species[s]
    p.push(l)
    total_prim_lipids.push(l)
end

# read secondary connections and unify for each species
lipids_nested = Hash.new{|h,k| h[k]=[]}
total_2nd_lipids =[]
File.open(lip2lip).each do |line|
    line.chomp
    cols = line.split(/\t/)
    l1 = cols[0]
    l2 = cols[1]
    pl = lipids_nested[l1]
    unless lipids_nested[l1].include?(l2)
        pl.push(l2)
    end
    total_2nd_lipids.push(l2)
end

# determine all uniq secondary connections
# 1. check if sec connection is part of any primary connection; if so, remove
puts ("removing secondary connections that are primaries of others ...")
lipids_nested_cleaned_1 = {}
lipids_nested.each do |l1, l2|
    lipids_nested_cleaned_1[l1] = l2 - total_prim_lipids
end

# 2. check if sec connection is part of any other secondary connection; if so, remove
# create a new hash, where species are the key and the secondary lipids are the values
puts ("removing shared secondary connections ...")
species2 = Hash.new{|h,k| h[k]=[]}
species.each do | spc, lip1 |
  p = species2[spc]
  lip1.each do |l11|
    lipids_nested_cleaned_1.each do |l1, l2|
      if l1 == l11 && l2.length > 0
        p.push(l2)
      end
    end
  end
end


# new clean hash with species
species2_cleaned = Hash.new{|h,k| h[k]=[]}
# hash with counts of each lipid across all species
hash_count_lipids = species2.values.flatten.inject(Hash.new(0)) { |total, e| total[e] += 1; total }
#puts hash_count_lipids
# check how often each secondary connection occurs in each species and compare to the counts across all species (if bigger or the same, keep!)
species2.each do |l1, l2|
  lip_sec = l2.flatten.inject(Hash.new(0)) { |total, e| total[e] += 1; total }
  lip_sec.each do | entry, amount |
    unless hash_count_lipids[entry] > amount
      val1 = l1
      p = species2_cleaned[l1]
      p.push(entry)
    end
  end
end

# create a network hash where species are the key and the values is an array of 2 arrays, one for primary lipids, one for secondary lipids
network = Hash.new{|h,k| h[k]=[[]]}
species.each do | spc, lip1|
    p = network[spc]
    p[0].push(lip1)
    network[spc][1]=species2_cleaned[spc]
end

puts ("... this amount of species were assigned lipids: #{species.keys.length}")
puts ("... this amount of lipids were primary: #{species.values.flatten.length}")
puts ("... this amount of lipids were secondary: #{species2_cleaned.flatten.length}")

# write network to file
File.open("org-lipids_network.txt", "w+") do | file |
    file.puts("organism\tprimary_lipid\tseconday_lipid") #\tscore")
    network.each do | spc, lipids |
        file.puts("#{spc}\t#{lipids[0].to_s}\t#{lipids[1].to_s}") #\t#{lipids[0].length/lipids[1].length}")
    end
end

# write primary lipids to file
File.open("primary_lipids.txt", "w+") do | file |
  file.puts total_prim_lipids
end

# write secondary lipids to file
File.open("secondary_lipids.txt", "w+") do | file |
  species2_cleaned.each do | l1, l2|
    l2.each do | secondlipid |
        file.puts secondlipid
    end
  end
end

puts ("finished!")
