import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# Load VCF file using pandas
vcf_file = 'example.vcf'
vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)

# Extract relevant columns for clustering analysis
vcf_df = vcf_df.iloc[:, 9:]
vcf_df.columns = ['sample' + str(i) for i in range(1, vcf_df.shape[1]+1)]

# Convert genotypes to numerical values
vcf_df = vcf_df.replace(['0/0', '0/1', '1/0', '1/1'], [0, 1, 1, 2])

# Calculate pairwise distance matrix
distance_matrix = np.zeros((vcf_df.shape[1], vcf_df.shape[1]))
for i in range(vcf_df.shape[1]):
    for j in range(i+1, vcf_df.shape[1]):
        distance = np.mean(np.abs(vcf_df.iloc[:, i] - vcf_df.iloc[:, j]))
        distance_matrix[i, j] = distance
        distance_matrix[j, i] = distance

# Perform clustering analysis using KMeans
num_clusters = 3
kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(distance_matrix)

# Plot the results
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
for i in range(num_clusters):
    cluster_samples = [vcf_df.columns[j] for j in range(len(kmeans.labels_)) if kmeans.labels_[j] == i]
    cluster_distances = distance_matrix[[j for j in range(len(kmeans.labels_)) if kmeans.labels_[j] == i], :]
    cluster_distances = cluster_distances[:, [j for j in range(len(kmeans.labels_)) if kmeans.labels_[j] == i]]
    plt.scatter(cluster_distances.flatten(), np.zeros(len(cluster_distances.flatten())), color=colors[i], alpha=0.5, label='Cluster ' + str(i+1))
plt.legend(loc='best')
plt.show()
