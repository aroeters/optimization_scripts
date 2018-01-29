###############
## FUNCTIONS ##
###############
avg_weight_plot <- function(averages1, sd1, max_scores, big_scape, grouping, title, ylabel) {
  # This is made for 9 different classes, alter code if more classes are added
  colorCodes = c("#e41a1c", "#377eb8", "#4daf4a")
  for (i in 1:9) {
      plot(averages1[,i+1] ~ averages1$weight, xlim=c(0,1), ylim=c(min(averages1-sd1), max(averages1+sd1)), col=colorCodes[1], lwd=3, type="l", main=paste(title, grouping[i], "vs avg normalized score" , sep=" "), xlab="weight", ylab=ylabel, cex=1.5)
      lines(max_scores[,i+1] ~ max_scores$weight, lwd=3, lty="longdash", col=colorCodes[1])
      try({
        # this is made for 9 different classes, alter code if more classes are added
        lines(averages1[,i+10] ~ averages1$weight, col=colorCodes[2], lwd=3, lty="dotdash")
        lines(max_scores[,i+10] ~ max_scores$weight, lwd=3, lty="dotted", col=colorCodes[2])
      }, silent=T)
      if (title == "JI") {
        abline(v=big_scape$JI[big_scape$groups_temp == grouping[i]], lwd=2)
      } else if (title == "DSS") {
        abline(v=big_scape$DSS[big_scape$groups_temp == grouping[i]], lwd=2)
      } else if (title == "AI") {
        abline(v=big_scape$AI[big_scape$groups_temp == grouping[i]], lwd=2)
      }
      polygon(border=NA, x=c(sd1$weight,rev(sd1$weight)), y=c(averages1[,i+1]-sd1[,i+1],rev(averages1[,i+1]+sd1[,i+1])), col=adjustcolor(colorCodes[3], alpha.f = .5))
      legend("topleft", legend = c("average", "max score", "new AP average", "new AP max"), lty = c(1, 5, 4, 3), lwd=3, col=c(colorCodes[1], colorCodes[1], colorCodes[2], colorCodes[colorCodes[1]]))
      }
}
########
## Calculates the avarage per class per weight
########
# data = all data, weight_class = JI, DSS or AI
calculate_avg_for_weight_and_class <- function(data, metric, grouping) {
  x = c("weight", grouping)
  mean_frame = as.data.frame(matrix(,0, length(x)))
  for (i in 0:100) {
    current_weight = i/100
    average_per_weight = c(current_weight)
    for (i in 1:length(grouping)) {
      current_data = data[data$grouping==grouping[i],]
      if (metric == "JI") {
        average_per_weight = c(average_per_weight, mean(current_data$normalized_score[current_data$JI==current_weight]))
      } else if (metric == "DSS") {
        average_per_weight = c(average_per_weight, mean(current_data$normalized_score[current_data$DSS==current_weight]))
      } else if (metric == "AI") {
        average_per_weight = c(average_per_weight, mean(current_data$normalized_score[current_data$AI==current_weight]))
      }
    }
    mean_frame = rbind(mean_frame, average_per_weight)
  }
  (names(mean_frame) <- x)
  mean_frame = mean_frame[complete.cases(mean_frame),]
  mean_frame <- mean_frame[!is.infinite(rowSums(mean_frame)),]
  rm(current_data, current_weight)
  return(mean_frame)
}
########
## Calculates the standard deviation per class per weight
########
calculate_sd_for_weight_and_class <- function(data, metric, grouping) {
  x = c("weight", grouping)
  sd_frame = as.data.frame(matrix(,0, length(x)))
  for (i in 0:100) {
    current_weight = i/100
    sd_per_weight = c(current_weight)
    for (i in 1:length(grouping)) {
      current_data = data[data$grouping==grouping[i],]
      if (metric == "JI") {
        sd_per_weight = c(sd_per_weight, sd(current_data$normalized_score[current_data$JI==current_weight]))
      } else if (metric == "DSS") {
        sd_per_weight = c(sd_per_weight, sd(current_data$normalized_score[current_data$DSS==current_weight]))
      } else if (metric == "AI") {
        sd_per_weight = c(sd_per_weight, sd(current_data$normalized_score[current_data$AI==current_weight]))
      }
    }
    sd_frame = rbind(sd_frame, sd_per_weight)
  }
  names(sd_frame) <- x
  sd_frame = sd_frame[complete.cases(sd_frame),]
  sd_frame <- sd_frame[!is.infinite(rowSums(sd_frame)),]
  rm(current_data, current_weight)
  return(sd_frame)
}
########
## Retrieves the max score per class per weight
########
get_max_scores <- function(data, metric, grouping) {
  x = c("weight", grouping)
  max_score_frame = as.data.frame(matrix(,0, length(x)))
  for (i in 0:100) {
    current_weight = i/100
    max_score_per_weight = c(current_weight)
    for (i in 1:length(grouping)) {
      current_data = data[data$grouping==grouping[i],]
      if (metric == "JI") {
        max_score_per_weight = c(max_score_per_weight, max(current_data$normalized_score[current_data$JI==current_weight]))
      } else if (metric == "DSS") {
        max_score_per_weight = c(max_score_per_weight, max(current_data$normalized_score[current_data$DSS==current_weight]))
      } else if (metric == "AI") {
        max_score_per_weight = c(max_score_per_weight, max(current_data$normalized_score[current_data$AI==current_weight]))
      }
    }
    max_score_frame = rbind(max_score_frame, max_score_per_weight)
  }
  names(max_score_frame) <- x
  max_score_frame = max_score_frame[complete.cases(max_score_frame),]
  max_score_frame <- max_score_frame[!is.infinite(rowSums(max_score_frame)),]
  rm(current_data, current_weight)
  return(max_score_frame)
}
########
## Retrieves the avg max score per class per weight window
########
get_avg_max_scores_in_window <- function(max_scores, window_size) {
  x = c("weight_min","weight_max", colnames(max_scores)[2:ncol(max_scores)])
  avg_max_scores = as.data.frame(matrix(,0,length(x)))
  for (i in 0:100) {
    current_weight = i/100
    current_data = max_scores[max_scores$weight >= current_weight & max_scores$weight <= current_weight+window_size,2:ncol(max_scores)]
    if (length(current_data[,1]) >= 3 && (current_weight+window_size) <= 1.0) {
      new_row = c(current_weight, current_weight+window_size, apply(current_data, 2, mean))
      avg_max_scores = rbind(avg_max_scores, new_row)
    }
  }
  names(avg_max_scores) <- x
  return(avg_max_scores)
}
########
## Retrieves the max scores per class
########
get_max <- function(x, top_x_of_best_scores=1) {
  max_scores_per_class <- data.frame(matrix(ncol=2,nrow=0))
  names(max_scores_per_class) <- c("group", "score")
  groups_temp <- unique(x$grouping)
    
  for (i in 1:length(groups_temp)) {
    for (z in 1:top_x_of_best_scores) {
      
      score <- max(x$normalized_score[x$grouping==groups_temp[i]])
      x = x[x$normalized_score!=score,]
      max_scores_per_class = rbind(max_scores_per_class, data.frame(group=groups_temp[i], score=as.numeric(score)))
    }
  }
    return(max_scores_per_class)
}
########
# Retrieves the best weights
########
get_best_weigths <- function(x, best_scores, nr_of_entries=1) {
  weights_max_score = data.frame(matrix(ncol=3,nrow=0))
  names(weights_max_score) <- c("group", "best_weights", "score")
  for (i in 1:length(best_scores$group)) {
    weights <- as.character(x$file[x$normalized==best_scores$score[i] & x$grouping==best_scores$group[i]][1:nr_of_entries])
    weights_max_score = rbind(weights_max_score, data.frame(group=best_scores$group[i], best_weights=weights,score=best_scores$score[i]))
  }
  return(weights_max_score)
}
########
## Retrieves the best scores per cutoff value per class
########
get_best_scores_per_cutoff <- function(x, top_x_of_best_scores=1, entries_per_score=1) {
  temp_data = x
  max_scores_per_class_per_cutoff <- data.frame(matrix(ncol=3,nrow=0))
  names(max_scores_per_class_per_cutoff) <- c("group", "cutoff", "score")
  groups_temp <- unique(x$grouping)
  cutoffs <- unique(temp_data$cutoff)
  for (i in 1:length(groups_temp)) {
    for (j in 1:length(cutoffs)) {
      for (k in 1:top_x_of_best_scores) {
      score <- max(temp_data$normalized_score[temp_data$grouping==groups_temp[i] & temp_data$cutoff==cutoffs[j]])
      if (score != 0) {
        temp_data = temp_data[temp_data$normalized_score!=score,]
        max_scores_per_class_per_cutoff = rbind(max_scores_per_class_per_cutoff, data.frame(group=groups_temp[i], cutoff=as.numeric(cutoffs[j]), score=as.numeric(score)))
        }
      }
    }
  }
  rm(temp_data)
  weights_max_score_per_cutoff = data.frame(matrix(ncol=4,nrow=0))
  names(weights_max_score_per_cutoff) <- c("group", "cutoff", "best_weights", "score")
  
  for (i in 1:length(max_scores_per_class_per_cutoff$group)) {
    for (j in 1:length(cutoffs)) {
    weights <- as.character(x$file[x$normalized==max_scores_per_class_per_cutoff$score[i] & x$grouping==max_scores_per_class_per_cutoff$group[i] & x$cutoff==cutoffs[j]][1:entries_per_score])
    weights_max_score_per_cutoff = rbind(weights_max_score_per_cutoff, data.frame(group=max_scores_per_class_per_cutoff$group[i],cutoff=cutoffs[j], best_weights=weights,score=max_scores_per_class_per_cutoff$score[i]))
    }
  }
  
  weights_max_score_per_cutoff = na.omit(weights_max_score_per_cutoff)

  weights_max_score_per_cutoff = weights_max_score_per_cutoff[order(weights_max_score_per_cutoff$group,weights_max_score_per_cutoff$cutoff),]
  return(weights_max_score_per_cutoff)
  
}


##########
## DATA ##
##########

base_groups = c("all-mix", "RiPPs", "PKS-NRP", "PKSI", "PKSother", "Saccharides", "Terpene", "NRPS", "Other")
base_groups = base_groups[order(base_groups)]
colors = rainbow(length(base_groups))

# Load in the data
infile = "E:/Data/files/results_scoring_min.txt" # MIBiG sparse min, scikit default
# infile = "E:Data/files/results_scoring_median.txt" # MIBiG sparse median, scikit default

data = read.table(infile, sep="\t", header=T)
# Get the groups of files (all the same i.e. RiPPs)
data[,15] <- vector(length=length(data$file))
colnames(data)[15] <- "grouping"
for (i in 1:length(base_groups)) {
  data[grepl(base_groups[i], data$file), 15] <- base_groups[i]
  data[which(data$grouping == base_groups[i] & grepl("new_", data$file)), 15] <- paste(base_groups[i], "_new", sep="")
}
groups = as.vector(unique(data$grouping))

# currently used weights in big-scape
JI =  c(0.20, 0.28, 0.00, 0.22, 0.00, 0.00, 0.20, 0.00, 0.01)
DSS = c(0.75, 0.71, 0.78, 0.76, 0.32, 0.00, 0.75, 1.00, 0.97)
AI =  c(0.05, 0.01, 0.22, 0.02, 0.68, 1.00, 0.05, 0.00, 0.02)
AB =  c(2.0, 1.0, 1.0, 1.0, 4.0, 1.0, 2.00, 4.0, 4.0)
groups_temp <- c("all-mix", "RiPPs", "PKS-NRP", "PKSI", "PKSother", "Saccharides", "Terpene", "NRPS", "Other")
big_scape_weights = data.frame(groups_temp, JI, DSS, AI, AB)
big_scape_weights = big_scape_weights[order(big_scape_weights$groups_temp),]
rm(JI,DSS,AI,AB)

###########
## PLOTS ##
###########
# plot the data (shows how stable the  normalized scores are of each class)
boxplot(data$normalized_score ~ data$grouping, col=colors, cex=.5, cex.axis=.7)

# plot in a different way (shows how the normalized score behaves compared to the unnormalized score)
for (i in 1:length(groups)) {
  if (i == 1) {
    plot(data$normalized_score[data$grouping==groups[i]] ~ data$score[data$grouping==groups[i]],xlab="score", ylab="Normalized score", col=colors[i], pch=16, xlim=c(min(data$score), max(data$score)), ylim=c(min(data$normalized_score),max(data$normalized_score)))
  } else {
    points(data$normalized_score[data$grouping==groups[i]] ~ data$score[data$grouping==groups[i]], col=colors[i], pch=16)
  }
}
legend("topright", inset=c(-0.15,0), legend=groups, title="Classes", cex=.75, col=colors, lty=c(1,1), lwd=3)

# JI
avg_jaccard = calculate_avg_for_weight_and_class(data, "JI", groups)
sd_jaccard = calculate_sd_for_weight_and_class(data, "JI", groups)
max_scores_jaccard =  get_max_scores(data, "JI", groups)
avg_weight_plot(avg_jaccard, sd_jaccard, max_scores_jaccard, big_scape_weights, groups, "JI", "avg normalized score")

# AI
avg_adjcency = calculate_avg_for_weight_and_class(data, "AI", groups)
sd_adjacency = calculate_sd_for_weight_and_class(data, "AI", groups)
max_scores_adjacency =  get_max_scores(data, "AI", groups)
avg_weight_plot(avg_adjcency, sd_adjacency, max_scores_adjacency, big_scape_weights, groups, "AI", "avg normalized score")

# DSS
avg_domain_similarity = calculate_avg_for_weight_and_class(data, "DSS", groups)
sd_domain_similarity = calculate_sd_for_weight_and_class(data, "DSS", groups)
max_scores_domain_similarity =  get_max_scores(data, "DSS", groups)
avg_weight_plot(avg_domain_similarity, sd_domain_similarity, max_scores_domain_similarity, big_scape_weights, groups, "DSS", "avg normalized score")

# Get the maximum score per class
max_score_per_class = c()
for (i in 1:length(groups)) {
  max_score_per_class = c(max_score_per_class, max(data$normalized_score[data$grouping==groups[i]]))
}

# and the best weights corresponding to the scores
weights_max_score = c()
for (i in 1:length(groups)) {
  weights_max_score = c(weights_max_score, as.character(data$file[data$normalized==max(data$normalized_score[data$grouping==groups[i]]) & data$grouping==groups[i]][1:5]))
}
weights_max_score

##########################################################
## The three different methods of AP results comparison ##
##########################################################
## sparse current
infile_sparse = "E:/R/sparse_current_results/results_scoring.txt" # sparse current, min, damping 0.8
data = read.table(infile_sparse, sep="\t", header=T)
# Get the groups of files (all the same i.e. RiPPs)
data[,16] <- vector(length=length(data$file))
colnames(data)[16] <- "grouping"
for (i in 1:length(base_groups)) {
  data[grepl(base_groups[i], data$file), 16] <- base_groups[i]
  # data[which(data$grouping == base_groups[i] & grepl("new_", data$file)), 15] <- paste(base_groups[i], "_new", sep="")
}
groups = as.vector(unique(data$grouping))

data_sparse <- data


####################################################################
## scikit
infile_scikit = "E:/R/scikit_results/results_scoring.txt" # scikit, default
data = read.table(infile_scikit, sep="\t", header=T)
# Get the groups of files (all the same i.e. RiPPs)
data[,16] <- vector(length=length(data$file))
colnames(data)[16] <- "grouping"
for (i in 1:length(base_groups)) {
  data[grepl(base_groups[i], data$file), 16] <- base_groups[i]
  data[which(data$grouping == base_groups[i] & grepl("new_", data$file)), 16] <- paste(base_groups[i], "_new", sep="")
}
groups = as.vector(unique(data$grouping))

data_scikit <- data


####################################################################
## sparse new
infile_sparse_new = "E:/R/sparse_new_results/results_scoring.txt" # sparse new, median, damping 0.5
data = read.table(infile_sparse_new, sep="\t", header=T)
# Get the groups of files (all the same i.e. RiPPs)
data[,16] <- vector(length=length(data$file))
colnames(data)[16] <- "grouping"
for (i in 1:length(base_groups)) {
  data[grepl(base_groups[i], data$file), 16] <- base_groups[i]
  }
groups = as.vector(unique(data$grouping))

data_sparse_new <- data

# get all best scores per class
max_score_sparse_current <- get_max(data_sparse, 2)
max_score_scikit <- get_max(data_scikit, 2)
max_score_sparse_new <- get_max(data_sparse_new, 2)

# put them all in a single matrix
scores <- max_score_sparse_current
scores <- cbind(scores, max_score_sparse_new$score)
scores <- cbind(scores, max_score_scikit$score)
names(scores) <- c("group", "current", "sparse_new", "scikit")
View(scores)

# get all best weights
best_weights_sparse_current <- get_best_weigths(data_sparse, max_score_sparse_current,20)
best_weights_scikit <- get_best_weigths(data_scikit, max_score_scikit,20)
best_weights_sparse_new <- get_best_weigths(data_sparse_new, max_score_sparse_new, 5)

# put all weights in a single matrix
best_weights <- data.frame(best_weights_sparse_current, best_weights_scikit$best_weights, best_weights_scikit$score)
names(best_weights) <- c("group", "current","score_current", "scikit", "score_scikit")
View(best_weights)

best_scores_per_cutoff_sparse <- get_best_scores_per_cutoff(data_sparse, 1, 1)
best_scores_per_cutoff_sparse_new <- get_best_scores_per_cutoff(data_sparse_new, 1, 1)
best_scores_per_cutoff_scikit <- get_best_scores_per_cutoff(data_scikit, 1, 1)
write.table(best_scores_per_cutoff_sparse, file="E:/Data/images/best_weights_sparse.tsv", sep=",\t", append = F, row.names=F)
write.table(best_scores_per_cutoff_sparse_new, file="E:/Data/images/best_weights_sparse_new.tsv", sep=",\t", append = F, row.names=F)
write.table(best_scores_per_cutoff_scikit, file="E:/Data/images/best_weights_scikit.tsv", sep=",\t", append = F, row.names=F)








##################################
## sparse current lower penalty ############################################################################################
##################################
infile_sparse_lower_penalty = "M:/Thesis/R/sparse_current_lower_penalty/results_scoring_lower_penalty.txt" # sparse current, min, damping 0.8
data_lower_penalty = read.table(infile_sparse_lower_penalty, sep="\t", header=T)
# Get the groups of files (all the same i.e. RiPPs)
data_lower_penalty[,16] <- vector(length=length(data_lower_penalty$file))
colnames(data_lower_penalty)[16] <- "grouping"
for (i in 1:length(base_groups)) {
  data_lower_penalty[grepl(base_groups[i], data_lower_penalty$file), 16] <- base_groups[i]
}
groups = as.vector(unique(data_lower_penalty$grouping))
best_scores_per_cutoff_lower_penalty <- get_best_scores_per_cutoff(data_lower_penalty, 3, 2)
write.table(best_scores_per_cutoff_lower_penalty, file="M:/Thesis/R/sparse_current_lower_penalty/best_weights_sparse_lower_penalty.tsv", sep=",\t", append = F, row.names=F)

############################
## sparse current entropy ##################################################################################################
############################
infile_sparse_entropy = "M:/Thesis/R/sparse_current_entropy/results_scoring_entropy.txt" # sparse current, min, damping 0.8
data_entropy = read.table(infile_sparse_entropy, sep="\t", header=T)
# Get the groups of files (all the same i.e. RiPPs)
data_entropy[,16] <- vector(length=length(data_entropy$file))
colnames(data_entropy)[16] <- "grouping"
for (i in 1:length(base_groups)) {
  data_entropy[grepl(base_groups[i], data_entropy$file), 16] <- base_groups[i]
}
groups = as.vector(unique(data_entropy$grouping))
best_scores_per_cutoff_entropy <- get_best_scores_per_cutoff(data_entropy, 3, 2)
write.table(best_scores_per_cutoff_entropy, file="M:/Thesis/R/sparse_current_entropy/best_weights_sparse_entropy.tsv", sep=",\t", append = F, row.names=F)

