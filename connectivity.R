data = read.csv("M:/Thesis/connectivity_results.tsv", sep="\t")
#plot(data[,3:10])

data = data[data$nr_of_clusters > 20,]

# Filter out the networks with a percentage of all BGCs in a single cluster
data = data[data$biggest_cluster_size < data$total_bgcs*.6, ] 
# Filter all networks with only a single cluster
data = data[data$smallest_cluster_size != data$biggest_cluster_size, ]
# Filter out the data with too much singletons
data = data[data$nr_of_singletons <= data$total_bgcs*(0.9-data$cutoff), ]

