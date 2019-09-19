# The MIT License (MIT)
# Copyright (c) 2019 Alexander J Probst

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

print ("... reading input files ...")

#input
library(gplots)
#setwd("/Users/ajp/Documents/Projects/UC_Berkeley/CG_data/CG_lipids/correlation_rpS3_lipids_2018Mar14")
m1 <- c("lipids.txt")
m2 <- c("rpS3.txt") # SORTED BY MAX VALUE OF EACH SPECIES ACROSS ALL SAMPLES ! ! !
dat1=t(read.table(m1, sep="\t", row.names=1, header=T))
dat2=t(read.table(m2, sep="\t", row.names=1, header=T))
#cormethod=c("pearson")

lipids=c()
species=c()
cor_vals=c()
q_vals=c()
p_vals=c()
scores=c()
max_abds=c()

print ("... correlating lipids with rpS3 ...")

for (spc1 in 1:length(dat1[1,])) {
  for (spc2 in 1:length(dat2[1,])) {
    correlation = cor.test(dat1[,spc1],dat2[,spc2])
    cor_val = as.vector(correlation$estimate)
    p_val = correlation$p.value
    q_val = p_val*length(dat2[1,]) # Bonferroni correction
    max_abd = max(dat2[,spc2])
    score = q_val/(max_abd)
    #print (q_val)
    #print (max(dat2[,spc2]))
    #print (score)
                          if ( (q_val < 0.05) && (cor_val > 0) ) {
                              #print(q_val)
                            lipids[length(lipids)+1] = colnames(dat1)[spc1]
                            species[length(species)+1] = colnames(dat2)[spc2]
                            cor_vals[length(cor_vals)+1] = cor_val
                            q_vals[length(q_vals)+1] = q_val
                            p_vals[length(p_vals)+1] = p_val
                            scores[length(scores)+1] = score
                            max_abds[length(max_abds)+1] = max_abd
                          }
  }
}

final_matrix=t(rbind(species,lipids,cor_vals,scores,q_vals,max_abds,p_vals))
write.table(final_matrix, "species2lipids_correlation.txt", quote=F, sep='\t', row.names=F)
#second_matrix=t(rbind(species,lipids,cor_vals))
#write.table(second_matrix, "species2lipids_cor.txt", quote=F, sep='\t', row.names=F)


print ("... correlating lipids with lipids ...")

# do correlation of lipids to themselves

lipids1=c()
lipids2=c()
cor_vals=c()
q_vals=c()
p_vals=c()

for (spc1 in 1:length(dat1[1,])) {
  for (spc2 in 1:length(dat1[1,])) {
    correlation = cor.test(dat1[,spc1],dat1[,spc2])
    cor_val = as.vector(correlation$estimate)
    p_val = correlation$p.value
    q_val = p_val*length(dat1[1,]) # Bonferroni correction
                          if (spc1 == spc2) {  } #  print ("identical hit") # remove hits to identical lipids
                          else if ( (q_val < 0.05) && (cor_val > 0) ) {
                              #print(q_val)
                            lipids1[length(lipids1)+1] = colnames(dat1)[spc1]
                            lipids2[length(lipids2)+1] = colnames(dat1)[spc2]
                            cor_vals[length(cor_vals)+1] = cor_val
                            q_vals[length(q_vals)+1] = q_val
                            p_vals[length(p_vals)+1] = p_val
                          }
  }
}

final_matrix=t(rbind(lipids1,lipids2,cor_vals,q_vals,p_vals))
write.table(final_matrix, "lipids2lipids_correlation.txt", quote=F, sep='\t', row.names=F)
#second_matrix=t(rbind(lipids2,lipids1,cor_vals))
#write.table(second_matrix, "lipids2lipids_cor.txt", quote=F, sep='\t', row.names=F)


print ("... correlations finished!")

print ("... parsing out best score table and corresponding first and second lipid connections for network analysis ...")

# some quick parsing done in bash
system ("cat species2lipids_correlation.txt | sed 1d | awk '{print$2}' | sort -u > uniqLipids.tmp")
system ("cat species2lipids_correlation.txt | sed 1d | sort -k4,4n > sorted_species2lipids.tmp")
system ("for lipid in $(cat uniqLipids.tmp); do grep -w -m 1 $lipid sorted_species2lipids.tmp >> best_hit_species2lipids.tmp; done")
system ("for lipid in $(cat uniqLipids.tmp); do grep -w \"^${lipid}\" lipids2lipids_correlation.txt >> first_second_lipids.tmp; done ")
system ("cut -f 1,2,3 species2lipids_correlation.txt | head -n 1 > file1.tmp")
system ("cut -f 1,2,3 best_hit_species2lipids.tmp > file2.tmp")
system ("cut -f 1,2,3 first_second_lipids.tmp > file3.tmp")
system ("cat file1.tmp file2.tmp file3.tmp > network_analysis_table.txt")
system ("mv file2.tmp species2lipids_best_hit.txt")
system ("mv file3.tmp first_second_lipid_connections.txt")
system ("rm uniqLipids.tmp sorted_species2lipids.tmp best_hit_species2lipids.tmp first_second_lipids.tmp file1.tmp")

print ("FINISHED!")
print ("Species against lipids are stored here: species2lipids_correlation.txt")
print ("Lipids against lipids are stored here: lipids2lipids_correlation.txt")
print ("File for network analysis of species to first and second connections of lipids is stored here: network_analysis_table.txt")

q()
