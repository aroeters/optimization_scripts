input_file = "E:/Data/files/f_statistic_scores.tsv"

# header for the files is this (not really in the file)
# file_name, f_score, cutoff, jaccard, adjacency index, domain similarity score, anchor boost
f_statistic_data = read.table(input_file, sep="\t")
colnames(f_statistic_data) <- c("file", "f_score", "cutoff", "Jaccard", "AI", "DSS", "AB")

# Separate data based on group
all_mix = f_statistic_data[grepl("all_mix", f_statistic_data[,1]), ]
RiPPs = f_statistic_data[grepl("allRiPPs", f_statistic_data[,1]), ]
PKS_NRP = f_statistic_data[grepl("allPKS-NRP", f_statistic_data[,1]), ]
PKSI = f_statistic_data[grepl("allPKSI", f_statistic_data[,1]), ]
PKSother = f_statistic_data[grepl("allPKSother", f_statistic_data[,1]), ]
Saccharides = f_statistic_data[grepl("allSaccharides", f_statistic_data[,1]), ]
Terpene = f_statistic_data[grepl("allTerpene", f_statistic_data[,1]), ]
NRPS = f_statistic_data[grepl("allNRPS", f_statistic_data[,1]), ]
Other = f_statistic_data[grepl("allOther", f_statistic_data[,1]), ]

# Order from high to low scores based on their cutoff and f-score
all_mix = all_mix[order(all_mix$cutoff, all_mix$f_score, decreasing = TRUE),]
RiPPs = RiPPs[order(RiPPs$cutoff, RiPPs$f_score, decreasing = TRUE),]
PKS_NRP = PKS_NRP[order(PKS_NRP$cutoff, PKS_NRP$f_score, decreasing = TRUE),]
PKSI = PKSI[order(PKSI$cutoff, PKSI$f_score, decreasing = TRUE),]
PKSother = PKSother[order(PKSother$cutoff, PKSother$f_score, decreasing = TRUE),]
Saccharides = Saccharides[order(Saccharides$cutoff, Saccharides$f_score, decreasing = TRUE),]
Terpene = Terpene[order(Terpene$cutoff, Terpene$f_score, decreasing = TRUE),]
NRPS = NRPS[order(NRPS$cutoff, NRPS$f_score, decreasing = TRUE),]
Other = Other[order(Other$cutoff, Other$f_score, decreasing = TRUE),]

# Write to files
write.table(all_mix, file="E:/f_score_results/all_mix.tsv", sep=",\t", append = F, row.names=F)
write.table(RiPPs, file="E:/f_score_results/RiPPs.tsv", sep=",\t", append = F, row.names=F)
write.table(PKS_NRP, file="E:/f_score_results/PKS_NRP.tsv", sep=",\t", append = F, row.names=F)
write.table(PKSI, file="E:/f_score_results/PKSI.tsv", sep=",\t", append = F, row.names=F)
write.table(PKSother, file="E:/f_score_results/PKSother.tsv", sep=",\t", append = F, row.names=F)
write.table(Saccharides, file="E:/f_score_results/Saccharides.tsv", sep=",\t", append = F, row.names=F)
write.table(Terpene, file="E:/f_score_results/Terpene.tsv", sep=",\t", append = F, row.names=F)
write.table(NRPS, file="E:/f_score_results/NRPS.tsv", sep=",\t", append = F, row.names=F)
write.table(Other, file="E:/f_score_results/Other.tsv", sep=",\t", append = F, row.names=F)

# Process the results
par(mfrow=c(3,3))
plot(all_mix$f_score, main = "all_mix")
plot(RiPPs$f_score, main = "RiPPs")
plot(PKS_NRP$f_score, main = "PKS_NRP")
plot(PKSI$f_score, main = "PKSI")
plot(PKSother$f_score, main = "PKSother")
plot(Saccharides$f_score, main = "Saccharides")
plot(Terpene$f_score, main = "Terpene")
plot(NRPS$f_score, main = "NRPS")
plot(Other$f_score, main = "Other")
par(mfrow=c(1,1))
# currently used weights in big-scape
JI =  c(0.20, 0.28, 0.00, 0.22, 0.00, 0.00, 0.20, 0.00, 0.01)
DSS = c(0.75, 0.71, 0.78, 0.76, 0.32, 0.00, 0.75, 1.00, 0.97)
AI =  c(0.05, 0.01, 0.22, 0.02, 0.68, 1.00, 0.05, 0.00, 0.02)
AB =  c(2.0, 1.0, 1.0, 1.0, 4.0, 1.0, 2.00, 4.0, 4.0)
groups_temp <- c("all-mix", "RiPPs", "PKS-NRP", "PKSI", "PKSother", "Saccharides", "Terpene", "NRPS", "Other")
big_scape_weights = data.frame(groups_temp, JI, DSS, AI, AB)
big_scape_weights = big_scape_weights[order(big_scape_weights$groups_temp),]
rm(JI,DSS,AI,AB,groups_temp)







