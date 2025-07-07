import sys
import os
import logging
logging.getLogger('tensorflow').disabled = True
from logger import *

from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS
# from sklearn.manifold import MDS
from hdbscan import HDBSCAN
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
np.random.seed(42)

import seaborn as sns
from umap import UMAP
import networkx as nx

import pandas as pd
# from sklearn.preprocessing import MinMaxScaler
pd.set_option('display.max_rows', None)  # Displays all rows
pd.set_option('display.max_columns', None)  # Displays all columns

import warnings
warnings.filterwarnings("ignore", message="n_jobs value .* overridden to 1 by setting random_state")
warnings.filterwarnings("ignore", message="using precomputed metric; inverse_transform will be unavailable")

from utils import get_family_id, create_directory, get_known_motif_shortcode


def combine_all_data(pdb_chainwise_loops, pdb_chainwise_centroids, spatial_proximity_data):
	combined_loops = []
	combined_centroids = []
	pdb_chain_by_loop = {}
	for pdb_id in pdb_chainwise_loops:
		for chain_id in pdb_chainwise_loops[pdb_id]:
			if len(pdb_chainwise_loops[pdb_id][chain_id]) >= 20:
				# print(len(pdb_chainwise_loops[pdb_id][chain_id]))
				# continue
				for loop in pdb_chainwise_loops[pdb_id][chain_id]:
					combined_loops.append(loop)
					combined_centroids.append(pdb_chainwise_centroids[pdb_id][chain_id][loop])
					pdb_chain_by_loop[loop] = (pdb_id, chain_id)
	# sys.exit()
	total_loops = len(combined_loops)
	combined_distance_matrix = np.full((total_loops, total_loops), 9999999.99)
	for i in range(total_loops):
		for j in range(i+1, total_loops):
			if i == j:
				combined_distance_matrix[i, j] = combined_distance_matrix[j, i] = 0.0
			else:
				pdb_id1, chain_id1 = pdb_chain_by_loop[combined_loops[i]]
				pdb_id2, chain_id2 = pdb_chain_by_loop[combined_loops[j]]
				if pdb_id1 == pdb_id2 and chain_id1 == chain_id2:
					combined_distance_matrix[i, j] = combined_distance_matrix[j, i] = spatial_proximity_data[pdb_id1][chain_id1][combined_loops[i]][combined_loops[j]]

	return combined_distance_matrix, combined_centroids, combined_loops


def precompute_distance_matrices(pdb_chainwise_loops, spatial_proximity_data):
	distance_matrices = {}
	for pdb_id in pdb_chainwise_loops:
		distance_matrices[pdb_id] = {}
		for chain_id in pdb_chainwise_loops[pdb_id]:
			loops = pdb_chainwise_loops[pdb_id][chain_id]
			n_motifs = len(loops)
			distances = np.zeros((n_motifs, n_motifs))
			for i in range(n_motifs):
				for j in range(i+1, n_motifs):
					if i == j:
						distances[i,j] = distances[j,i] = 0.0
					else:
						distances[i,j] = distances[j,i] = spatial_proximity_data[pdb_id][chain_id][loops[i]][loops[j]]
			distance_matrices[pdb_id][chain_id] = distances

	return distance_matrices

def get_centroid(coordinates):
	coordinates_array = np.array(coordinates)
	centroid = np.mean(coordinates_array, axis=0)
	return list(centroid)

def calculate_centroids(pdb_chainwise_residue_data):
	pdb_chainwise_centroids = {}

	for pdb_id in pdb_chainwise_residue_data:
		pdb_chainwise_centroids[pdb_id] = {}

		for chain_id in pdb_chainwise_residue_data[pdb_id]:
			pdb_chainwise_centroids[pdb_id][chain_id] = {}

			(loop_info_dict, min_coords, max_coords) = pdb_chainwise_residue_data[pdb_id][chain_id]
			for loop in loop_info_dict:
				loop_coords = loop_info_dict[loop]
				# print(loop_coords)
				res_centroid_coords = []
				for res in loop_coords:
					res_centroid_coords.append(get_centroid(res))
				centroid = get_centroid(res_centroid_coords)
				# print(centroid)
				# sys.exit()
				pdb_chainwise_centroids[pdb_id][chain_id][loop] = centroid

	return pdb_chainwise_centroids



def visualize_cluster_3DPlot(families, pdb_id, chain_id, pdb_chainwise_loops, pdb_chainwise_centroids, cluster_labels):
	loops = pdb_chainwise_loops[pdb_id][chain_id]
	motif_centroids = []
	motif_families = []
	for i in range(len(loops)):
		motif_centroids.append(pdb_chainwise_centroids[pdb_id][chain_id][loops[i]])
		motif_families.append(get_family_id(loops[i], families))

	motif_centroids = np.array(motif_centroids)


	# loops = pdb_chainwise_loops[pdb_id][chain_id]

	unique_families = list(families.keys())
	family_to_color = {family: plt.cm.tab20(i/len(unique_families))
	for i, family in enumerate(unique_families)}

	fig = plt.figure(figsize=(12, 10))
	ax = fig.add_subplot(111, projection='3d')

	for family in unique_families:
		# print('showing')
		# print(family)
		# Get indices of motifs belonging to this family

		# for i, loop in enumerate(loops):
		# 	get_family_id(loop, families)
		# family_indices = [i for i, f in enumerate(motif_families) if f == family]
		family_indices = [i for i, loop in enumerate(loops) if get_family_id(loop, families) == family]
		# print(family_indices)
		# print([loops[i] for i in family_indices])
		
		# Get their coordinates
		family_points = []
		for i in family_indices:
			family_points.append(pdb_chainwise_centroids[pdb_id][chain_id][loops[i]])
		# family_points = [pdb_chainwise_centroids[pdb_id][chain_id][loops[i]] for i in family_indices]
		family_points = np.array(family_points)
		# print(family_points)
		
		# Get their cluster assignments
		family_clusters = [cluster_labels[i] for i in family_indices]
		# print(family_clusters)
		
		# Use different markers for noise (-1) vs clustered points
		clustered_indices = [i for i, c in enumerate(family_clusters) if c != -1]
		noise_indices = [i for i, c in enumerate(family_clusters) if c == -1]

		# print(clustered_indices)
		# print(noise_indices)
		
		# Plot clustered points with circles
		if clustered_indices:
			ax.scatter(
				family_points[clustered_indices, 0],
				family_points[clustered_indices, 1],
				family_points[clustered_indices, 2],
				color=family_to_color[family],
				marker='o',
				label=f'{family}',
				s=80,
				alpha=0.7
			)
		
		# Plot noise points with x markers
		if noise_indices:
			ax.scatter(
				family_points[noise_indices, 0],
				family_points[noise_indices, 1],
				family_points[noise_indices, 2],
				color=family_to_color[family],
				marker='x',
				s=80,
				alpha=0.4
			)

	for cluster_id in sorted(set(cluster_labels)):
		if cluster_id == -1:  # Skip noise
			continue
		
		# Get points in this cluster
		cluster_indices = np.where(cluster_labels == cluster_id)[0]
		cluster_points = motif_centroids[cluster_indices]
		
		# Add a light boundary or annotation
		ax.text(
			np.mean(cluster_points[:, 0]),
			np.mean(cluster_points[:, 1]),
			np.mean(cluster_points[:, 2]),
			f'Cluster {cluster_id}',
			fontsize=10,
			bbox=dict(facecolor='white', alpha=0.5)
		)

	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.set_title('RNA Motif Clustering by Family')
	plt.legend(title="Motif Family")

	# Create a secondary legend for cluster/noise distinction
	# from matplotlib.lines import Line2D
	# legend_elements = [
	# 	Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', label='In cluster', markersize=10),
	# 	Line2D([0], [0], marker='x', color='gray', label='Noise', markersize=10)
	# ]
	# ax.legend(handles=legend_elements, title="Cluster status", loc='upper right')

	plt.tight_layout()
	plt.show()

def visualize_cluster_PCA(families, pdb_id, chain_id, pdb_chainwise_loops, pdb_chainwise_centroids, cluster_labels):
	loops = pdb_chainwise_loops[pdb_id][chain_id]
	motif_centroids = []
	motif_families = []
	for i in range(len(loops)):
		motif_centroids.append(pdb_chainwise_centroids[pdb_id][chain_id][loops[i]])
		motif_families.append(get_family_id(loops[i], families))

	motif_centroids = np.array(motif_centroids)

	# Apply PCA
	pca = PCA(n_components=2)
	pca_result = pca.fit_transform(motif_centroids)

	# Plot with family colors and cluster markers
	plt.figure(figsize=(10, 8))

	unique_families = list(families.keys())
	family_to_color = {family: plt.cm.tab20(i/len(unique_families))
	for i, family in enumerate(unique_families)}

	print(family_to_color)

	for family in unique_families:
		family_indices = [i for i, f in enumerate(motif_families) if f == family]
		
		# Split into clustered vs noise
		clustered = [i for i in family_indices if cluster_labels[i] != -1]
		noise = [i for i in family_indices if cluster_labels[i] == -1]
		
		# Plot clustered points
		plt.scatter(
			pca_result[clustered, 0], pca_result[clustered, 1],
			color=family_to_color[family],
			marker='o',
			label=family,
			s=80,
			alpha=0.7
		)
		
		# Plot noise points
		plt.scatter(
			pca_result[noise, 0], pca_result[noise, 1],
			color=family_to_color[family],
			marker='x',
			s=80,
			alpha=0.4
		)

	# Add cluster annotations
	for cluster_id in sorted(set(cluster_labels)):
		if cluster_id == -1:
			continue
		
		cluster_indices = np.where(cluster_labels == cluster_id)[0]
		centroid = np.mean(pca_result[cluster_indices], axis=0)
		
		plt.annotate(
			f'Cluster {cluster_id}',
			centroid,
			bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.5)
		)

	plt.title('PCA Visualization of RNA Motif Clusters')
	plt.xlabel('Principal Component 1')
	plt.ylabel('Principal Component 2')
	plt.legend(title="Motif Family")
	plt.tight_layout()
	plt.show()

def clustering_analysis(families, pdb_chainwise_loops, pdb_chainwise_residue_data, spatial_proximity_data, distance_threshold_to_be_nearest_residue, output_dir):

	print('')
	logger.info('Clustering motifs in each RNA chain using DBSCAN and HDBSCAN based on distance to get an insight from clustering point of view.')
	logger.info('PCA (both 2D and 3D), UMAP, and t-SNE are being used for visualizing the clusters. Moreover, a heatmap will be generated for each chain to show the family distribution in the clustering output.')
	print('')

	distance_matrices = precompute_distance_matrices(pdb_chainwise_loops, spatial_proximity_data)
	pdb_chainwise_centroids = calculate_centroids(pdb_chainwise_residue_data)
	

	# combined test:
	# outcome and decision: combining chains makes them appear as separate clusters in the clustering output, therefore, they should not be combined while generating the clusters.
	# combined_distance_matrix, combined_centroids, combined_loops = combine_all_data(pdb_chainwise_loops, pdb_chainwise_centroids, spatial_proximity_data)
	# combined_clustering2 = HDBSCAN(min_cluster_size=2, min_samples=None, cluster_selection_epsilon=distance_threshold_to_be_nearest_residue, metric='precomputed').fit(combined_distance_matrix)
	# combined_cluster_labels2 = combined_clustering2.labels_
	# combined_motif_families, combined_unique_families, combined_unique_clusters, combined_family_color_map, combined_cluster_marker_map = generate_preliminary_objects(combined_loops, families, combined_cluster_labels2)
	# plot_umap_with_family_colors(combined_motif_families, combined_unique_families, combined_unique_clusters, combined_distance_matrix, combined_cluster_labels2, combined_cluster_marker_map, combined_family_color_map)
	# plot_tsne_with_family_colors(combined_motif_families, combined_unique_families, combined_unique_clusters, combined_distance_matrix, combined_cluster_labels2, combined_cluster_marker_map, combined_family_color_map)
	# # plot_network_with_family_colors_cluster_markers(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map, distance_threshold_to_be_nearest_residue)
	# # plot_circular_layout_with_family_colors_cluster_markers(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map)
	# plot_cluster_family_distribution(combined_motif_families, combined_unique_families, combined_unique_clusters, combined_distance_matrix, combined_cluster_labels2, combined_cluster_marker_map, combined_family_color_map)
	# plot_detailed_umap(combined_motif_families, combined_unique_families, combined_unique_clusters, combined_distance_matrix, combined_cluster_labels2, combined_cluster_marker_map, combined_family_color_map)
	# sys.exit()



	plots_dir = os.path.join(output_dir, 'plots')
	create_directory(plots_dir)
	cluster_dir = os.path.join(plots_dir, 'unsupervised_clustering')
	create_directory(cluster_dir)
	dbscan_dir = os.path.join(cluster_dir, 'DBSCAN')
	create_directory(dbscan_dir)
	hdbscan_dir = os.path.join(cluster_dir, 'HDBSCAN')
	create_directory(hdbscan_dir)

	for pdb_id in distance_matrices:
		for chain_id in distance_matrices[pdb_id]:
			pdb_chain = str(pdb_id) + '_' + str(chain_id)

			distances = distance_matrices[pdb_id][chain_id]
			if len(pdb_chainwise_loops[pdb_id][chain_id]) < 5:
				continue
			# print(pdb_id, chain_id)
			# print(pdb_chainwise_loops[pdb_id][chain_id])
			# print(len(pdb_chainwise_loops[pdb_id][chain_id]))

			# DBSCAN test
			clustering_method = 'DBSCAN'
			clustering1 = DBSCAN(eps=distance_threshold_to_be_nearest_residue, min_samples=2, metric='precomputed').fit(distances)
			cluster_labels1 = clustering1.labels_

			# print(cluster_labels1)
			motif_families, unique_families, unique_clusters, family_color_map, cluster_marker_map = generate_preliminary_objects(pdb_chainwise_loops[pdb_id][chain_id], families, cluster_labels1)

			cluster_output_fname = os.path.join(dbscan_dir, pdb_chain + '_DBSCAN.tsv')
			fp = open(cluster_output_fname, 'w')
			fp.write('Motif\tFamilyName\tClusterLabel\n')
			for i, loop in enumerate(pdb_chainwise_loops[pdb_id][chain_id]):
				fp.write(str(loop) + '\t' + str(motif_families[i]) + '\t' + str(cluster_labels1[i]) + '\n')
			fp.close()

			# visualize_cluster_3DPlot(families, pdb_id, chain_id, pdb_chainwise_loops, pdb_chainwise_centroids, cluster_labels1)

			plot_output_fname = os.path.join(dbscan_dir, pdb_chain + '_pca.png')
			plot_PCA_with_family_colors(motif_families, unique_families, unique_clusters, distances, cluster_labels1, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			# visualize_cluster_PCA(families, pdb_id, chain_id, pdb_chainwise_loops, pdb_chainwise_centroids, cluster_labels1)
			plot_output_fname = os.path.join(dbscan_dir, pdb_chain + '_pca3D.png')
			plot_PCA3D_with_family_colors(motif_families, unique_families, unique_clusters, distances, cluster_labels1, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)

			if len(pdb_chainwise_loops[pdb_id][chain_id]) >= 15:
				plot_output_fname = os.path.join(dbscan_dir, pdb_chain + '_umap.png')
				plot_umap_with_family_colors(motif_families, unique_families, unique_clusters, distances, cluster_labels1, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			else:
				print('Skipping ' + clustering_method + ' clustering output visualisation with UMAP for ' + pdb_chain + ' due to having too few (<15) motifs. This chain contains ' + str(len(pdb_chainwise_loops[pdb_id][chain_id])) + ' motifs.')
			plot_output_fname = os.path.join(dbscan_dir, pdb_chain + '_tsne.png')
			plot_tsne_with_family_colors(motif_families, unique_families, unique_clusters, distances, cluster_labels1, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			plot_output_fname = os.path.join(dbscan_dir, pdb_chain + '_DBSCAN_cluster_family_distribution.png')
			plot_cluster_family_distribution(motif_families, unique_families, unique_clusters, distances, cluster_labels1, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			# plot_output_fname = os.path.join(dbscan_dir, pdb_chain + '_umap_detailed.png')
			# plot_detailed_umap(motif_families, unique_families, unique_clusters, distances, cluster_labels1, cluster_marker_map, family_color_map, plot_output_fname)
			# sys.exit()


			# HDBSCAN test
			clustering_method = 'HDBSCAN'
			clustering2 = HDBSCAN(min_cluster_size=2, min_samples=None, cluster_selection_epsilon=distance_threshold_to_be_nearest_residue, metric='precomputed').fit(distances)
			# print(clustering)
			cluster_labels2 = clustering2.labels_
			# print(cluster_labels2)
			# print(type(cluster_labels1))
			# visualize_cluster_3DPlot(families, pdb_id, chain_id, pdb_chainwise_loops, pdb_chainwise_centroids, cluster_labels1)
			# visualize_cluster_PCA(families, pdb_id, chain_id, pdb_chainwise_loops, pdb_chainwise_centroids, cluster_labels1)
			motif_families, unique_families, unique_clusters, family_color_map, cluster_marker_map = generate_preliminary_objects(pdb_chainwise_loops[pdb_id][chain_id], families, cluster_labels2)
			
			cluster_output_fname = os.path.join(hdbscan_dir, pdb_chain + '_HDBSCAN.tsv')
			fp = open(cluster_output_fname, 'w')
			fp.write('Motif\tFamilyName\tClusterLabel\n')
			for i, loop in enumerate(pdb_chainwise_loops[pdb_id][chain_id]):
				fp.write(str(loop) + '\t' + str(motif_families[i]) + '\t' + str(cluster_labels1[i]) + '\n')
			fp.close()

			# np.random.seed(42)
			plot_output_fname = os.path.join(hdbscan_dir, pdb_chain + '_pca.png')
			plot_PCA_with_family_colors(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			plot_output_fname = os.path.join(hdbscan_dir, pdb_chain + '_pca3D.png')
			plot_PCA3D_with_family_colors(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			# plot_umap_with_family_colors(families, pdb_chainwise_loops[pdb_id][chain_id], distances, cluster_labels2)

			if len(pdb_chainwise_loops[pdb_id][chain_id]) >= 15:
				plot_output_fname = os.path.join(hdbscan_dir, pdb_chain + '_umap.png')
				plot_umap_with_family_colors(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			else:
				print('Skipping ' + clustering_method + ' clustering output visualisation with UMAP for ' + pdb_chain + ' due to having too few (<15) motifs. This chain contains ' + str(len(pdb_chainwise_loops[pdb_id][chain_id])) + ' motifs.')
			plot_output_fname = os.path.join(hdbscan_dir, pdb_chain + '_tsne.png')
			plot_tsne_with_family_colors(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			# plot_network_with_family_colors_cluster_markers(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map, distance_threshold_to_be_nearest_residue)
			# plot_circular_layout_with_family_colors_cluster_markers(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map)
			plot_output_fname = os.path.join(hdbscan_dir, pdb_chain + '_HDBSCAN_cluster_family_distribution.png')
			plot_cluster_family_distribution(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain)
			# plot_output_fname = os.path.join(hdbscan_dir, pdb_chain + '_umap_detailed.png')
			# plot_detailed_umap(motif_families, unique_families, unique_clusters, distances, cluster_labels2, cluster_marker_map, family_color_map, plot_output_fname)
			# sys.exit()
	# test DBSCAN
	# test HDBSCAN

def generate_preliminary_objects(loops, families, cluster_labels):
	motif_families = []
	for i in range(len(loops)):
		motif_families.append(get_family_id(loops[i], families))
	unique_families = list(set(motif_families))

	family_color_map = dict(zip(unique_families, 
							  sns.color_palette("husl", len(unique_families))))

	marker_styles = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'H']#, '+']
	unique_clusters = np.unique(cluster_labels)
	cluster_marker_map = {}
	for i, cluster in enumerate(unique_clusters):
		if cluster == -1:
			continue
		cluster_marker_map[cluster] = marker_styles[i % len(marker_styles)]
	cluster_marker_map[-1] = 'x'

	motif_families = np.array(motif_families)

	return motif_families, unique_families, unique_clusters, family_color_map, cluster_marker_map









def plot_PCA3D_with_family_colors(motif_families, unique_families, unique_clusters, distance_matrix, labels, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain):
	n_samples = len(labels)
	mds = MDS(n_components=n_samples-1, dissimilarity='precomputed', random_state=42, normalized_stress='auto')
	coordinates = mds.fit_transform(distance_matrix)

	# Step 3: Apply PCA to the MDS coordinates, but now with 3 components
	# -----------------------------------------------------------------
	pca = PCA(n_components=3)  # Use 3 components for 3D visualization
	X_pca = pca.fit_transform(coordinates)

	# Step 4: Create a DataFrame with all the information
	# -------------------------------------------------
	df = pd.DataFrame({
		'PC1': X_pca[:, 0],
		'PC2': X_pca[:, 1],
		'PC3': X_pca[:, 2],  # Now we have a third component
		'Cluster': [f'Cluster {l}' if l >= 0 else 'Noise' for l in labels],
		'Family': motif_families,
		'ClusterNum': labels
	})

	# Step 5: Create the 3D visualization with two legend sets
	# ------------------------------------------------------
	fig = plt.figure(figsize=(14, 12))
	ax = fig.add_subplot(111, projection='3d')  # Create a 3D subplot

	# Plot each point with appropriate color (family) and marker (cluster)
	for family in unique_families:
		for cluster in unique_clusters:
			# Get points that belong to this family and cluster
			mask = (df['Family'] == family) & (df['ClusterNum'] == cluster)
			if sum(mask) > 0:  # Only plot if we have points in this combination
				ax.scatter(
					df.loc[mask, 'PC1'],
					df.loc[mask, 'PC2'],
					df.loc[mask, 'PC3'],  # Use the third dimension
					color=family_color_map[family],
					marker=cluster_marker_map[cluster],
					s=100,
					alpha=0.8,
					# edgecolors='k',
					linewidths=0.5,
					label=f"{family}-{cluster}"  # We won't use this label directly
				)

	# Add point labels if needed
	for i in range(len(df)):
		ax.text(
			df.iloc[i]['PC1'], 
			df.iloc[i]['PC2'], 
			df.iloc[i]['PC3'],
			str(i),  # You can replace this with actual motif IDs if available
			size=8, 
			zorder=1
		)

	# Create custom legends - one for families, one for clusters
	# Family legend (colors)
	family_legend_elements = [
		Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=family)
		for family, color in family_color_map.items() if family in df['Family'].values
	]

	# Cluster legend (markers)
	cluster_legend_elements = [
		Line2D([0], [0], marker=marker, color='none', markeredgecolor='black', 
           markerfacecolor='black', markersize=10, label=f'Cluster {cluster}' if cluster >= 0 else 'Noise')
		for cluster, marker in cluster_marker_map.items() if cluster in df['ClusterNum'].values
	]

	# Add legends to the plot
	family_legend = ax.legend(
		handles=family_legend_elements,
		title="RNA Families",
		loc='upper left',
		bbox_to_anchor=(1.05, 1)
	)

	# Add the second legend
	ax.add_artist(family_legend)  # Add the first legend to the plot
	ax.legend(
		handles=cluster_legend_elements,
		title=clustering_method + " Clusters",
		loc='upper left',
		bbox_to_anchor=(1.05, 0.5)
	)

	# Add explained variance information
	explained_var = pca.explained_variance_ratio_
	ax.set_xlabel(f'PC1 ({explained_var[0]:.2%} variance)')
	ax.set_ylabel(f'PC2 ({explained_var[1]:.2%} variance)')
	ax.set_zlabel(f'PC3 ({explained_var[2]:.2%} variance)')  # Z-axis label with explained variance

	# ax.set_title('RNA Structural Motifs Clustering Visualization using ' + clustering_method + ' (' + pdb_chain + ')', fontsize=16)
	plt.title('Clustering motifs in ' + pdb_chain + ' with ' + clustering_method + ' and visualizing output with 3D PCA', fontsize=16)
	ax.grid(True)

	# Add total explained variance as text annotation
	total_var = sum(explained_var)
	ax.text2D(0.05, 0.95, f'Total variance explained: {total_var:.2%}', 
			  transform=ax.transAxes, fontsize=12)

	plt.tight_layout()
	plt.subplots_adjust(right=0.75)  # Make room for the legends

	# plt.show()
	plt.savefig(plot_output_fname)
	plt.close()


def plot_PCA_with_family_colors(motif_families, unique_families, unique_clusters, distance_matrix, labels, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain):
	n_samples = len(labels)
	mds = MDS(n_components=n_samples-1, dissimilarity='precomputed', random_state=42, normalized_stress='auto')
	coordinates = mds.fit_transform(distance_matrix)

	# Step 3: Apply PCA to the MDS coordinates
	pca = PCA(n_components=2)
	X_pca = pca.fit_transform(coordinates)

	# Create a DataFrame for easier plotting
	df = pd.DataFrame({
		'PC1': X_pca[:, 0],
		'PC2': X_pca[:, 1],
		'Cluster': [f'Cluster {l}' if l >= 0 else 'Noise' for l in labels],
		'Family': motif_families,
		'ClusterNum': labels
	})

	# Create the visualization
	# fig, ax = plt.subplots(figsize=(10, 8))
	plt.figure(figsize=(12, 10))

	
	# Plot each point with appropriate color (family) and marker (cluster)
	for family in unique_families:
		for cluster in unique_clusters:
			# Get points that belong to this family and cluster
			mask = (df['Family'] == family) & (df['ClusterNum'] == cluster)
			if sum(mask) > 0:  # Only plot if we have points in this combination
				plt.scatter(
					df.loc[mask, 'PC1'],
					df.loc[mask, 'PC2'],
					color=family_color_map[family],
					marker=cluster_marker_map[cluster],
					s=100,
					alpha=0.8,
					# edgecolors='k',
					linewidths=0.5,
					label=f"{family}-{cluster}"  # We won't use this label directly
				)

	# Add point labels if needed
	for i in range(len(df)):
		plt.annotate(
			i,  # You can replace this with actual motif IDs if available
			(df.iloc[i]['PC1'], df.iloc[i]['PC2']),
			textcoords="offset points",
			xytext=(0, 5),
			ha='center',
			fontsize=8
		)

	# Create custom legends - one for families, one for clusters
	# Family legend (colors)
	family_legend_elements = [
		Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=family)
		for family, color in family_color_map.items() if family in df['Family'].values
	]

	# Cluster legend (markers)
	cluster_legend_elements = [
		Line2D([0], [0], marker=marker, color='none', markeredgecolor='black', 
           markerfacecolor='black', markersize=10, label=f'Cluster {cluster}' if cluster >= 0 else 'Noise')
		for cluster, marker in cluster_marker_map.items() if cluster in df['ClusterNum'].values
	]

	# Add legends to the plot
	family_legend = plt.legend(
		handles=family_legend_elements,
		title="RNA Families",
		loc='upper left',
		bbox_to_anchor=(1.05, 1)
	)

	# Add the second legend
	plt.gca().add_artist(family_legend)  # Add the first legend to the plot
	plt.legend(
		handles=cluster_legend_elements,
		title=clustering_method + " Clusters",
		loc='upper left',
		bbox_to_anchor=(1.05, 0.5)
	)

	# Add explained variance information
	explained_var = pca.explained_variance_ratio_
	plt.xlabel(f'PC1 ({explained_var[0]:.2%} variance)')
	plt.ylabel(f'PC2 ({explained_var[1]:.2%} variance)')

	# plt.title('RNA Structural Motifs Clustering Visualization using ' + clustering_method + ' (' + pdb_chain + ')', fontsize=16)
	plt.title('Clustering motifs in ' + pdb_chain + ' with ' + clustering_method + ' and visualizing output with 2D PCA', fontsize=16)
	plt.grid(True, linestyle='--', alpha=0.3)
	plt.tight_layout()
	plt.subplots_adjust(right=0.78)  # Make room for the legends

	# plt.show()
	plt.savefig(plot_output_fname)
	plt.close()


# 1. UMAP Visualization with both cluster and family information
def plot_umap_with_family_colors(motif_families, unique_families, unique_clusters, distance_matrix, labels, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain):
	# Reduce dimensions using UMAP
	umap_reducer = UMAP(metric='precomputed', random_state=42)
	umap_embedding = umap_reducer.fit_transform(distance_matrix)
	
	# Create figure with two subplots
	fig, ax = plt.subplots(figsize=(12, 10))

	# Plot colored by motif family
	for family in unique_families:
		for cluster in unique_clusters:
			# Find points that belong to both this family and this cluster
			mask = (motif_families == family) & (labels == cluster)
			if np.any(mask):  # Only plot if there are points in this combination
				ax.scatter(
					umap_embedding[mask, 0], 
					umap_embedding[mask, 1],
					marker=cluster_marker_map[cluster],
					color=family_color_map[family],
					s=100, 
					alpha=0.8,
					label=f"{family} (Cluster {cluster})"
				)

	family_handles = [plt.Line2D([0], [0], marker='o', color=family_color_map[family], 
								 linestyle='', markersize=10, label=family) 
					  for family in sorted(unique_families)]
	
	# Second, create custom handles for cluster markers (all in black for clarity)
	cluster_handles = [plt.Line2D([0], [0], marker=cluster_marker_map[cluster], 
								 color='black', linestyle='', markersize=10,
								 # label=f"Cluster {cluster}")
								 label=f'Cluster {cluster}' if cluster >= 0 else 'Noise') 
					  for cluster in unique_clusters]
	
	# Place the family legend on the right
	first_legend = ax.legend(handles=family_handles, title="Motif Family", 
						   loc='upper left', bbox_to_anchor=(1.05, 1))
	ax.add_artist(first_legend)
	
	# Place the cluster legend below the family legend
	ax.legend(handles=cluster_handles, title= clustering_method + " Cluster", 
			 loc='upper left', bbox_to_anchor=(1.05, 0.6))
	
	# plt.title('UMAP: RNA Motifs Colored by Family, Markers by Cluster (' + pdb_chain + ')')
	plt.title('Clustering motifs in ' + pdb_chain + ' with ' + clustering_method + ' and visualizing output with UMAP', fontsize=16)
	plt.tight_layout()
	plt.subplots_adjust(right=0.8)
	# return fig
	# plt.show()
	plt.savefig(plot_output_fname)
	plt.close()

# 2. t-SNE visualization with family colors
def plot_tsne_with_family_colors(motif_families, unique_families, unique_clusters, distance_matrix, labels, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain):
	n_samples = len(labels)
	# Reduce dimensions using t-SNE
	tsne = TSNE(n_components=2, init='random', metric='precomputed', random_state=42, perplexity=min(30, n_samples//5))
	tsne_embedding = tsne.fit_transform(distance_matrix)
	
	# Create figure
	fig, ax = plt.subplots(figsize=(12, 10))
	
	# Plot each point with color by family and marker by cluster
	for family in unique_families:
		for cluster in unique_clusters:
			# Find points that belong to both this family and this cluster
			mask = (motif_families == family) & (labels == cluster)
			if np.any(mask):  # Only plot if there are points in this combination
				ax.scatter(
					tsne_embedding[mask, 0], 
					tsne_embedding[mask, 1],
					marker=cluster_marker_map[cluster],
					color=family_color_map[family],
					s=100, 
					alpha=0.8,
					label=f"{family} (Cluster {cluster})"
				)
	
	# Create separate legends for families and clusters
	family_handles = [plt.Line2D([0], [0], marker='o', color=family_color_map[family], 
								 linestyle='', markersize=10, label=family)
					  for family in sorted(unique_families)]
	
	cluster_handles = [plt.Line2D([0], [0], marker=cluster_marker_map[cluster], 
								 color='black', linestyle='', markersize=10,
								 # label=f"Cluster {cluster}") 
								label=f'Cluster {cluster}' if cluster >= 0 else 'Noise')
					  for cluster in unique_clusters]
	
	first_legend = ax.legend(handles=family_handles, title="Motif Family", 
						   loc='upper left', bbox_to_anchor=(1.05, 1))
	ax.add_artist(first_legend)
	
	ax.legend(handles=cluster_handles, title=clustering_method + " Cluster", 
			 loc='upper left', bbox_to_anchor=(1.05, 0.6))
	
	# plt.title('t-SNE: RNA Motifs Colored by Family, Markers by Cluster (' + pdb_chain + ')')
	plt.title('Clustering motifs in ' + pdb_chain + ' with ' + clustering_method + ' and visualizing output with t-SNE', fontsize=16)
	plt.tight_layout()
	plt.subplots_adjust(right=0.8)
	# plt.show()
	plt.savefig(plot_output_fname)
	plt.close()


# 3. Network graph with family colors and cluster markers
# def plot_network_with_family_colors_cluster_markers(motif_families, unique_families, unique_clusters, distance_matrix, labels, cluster_marker_map, family_color_map, distance_threshold=0.3):
# 	# Create a graph
# 	G = nx.Graph()
	
# 	# Add nodes with family and cluster attributes
# 	for i in range(len(distance_matrix)):
# 		G.add_node(i, cluster=int(labels[i]), family=motif_families[i])
		
# 	# Add edges for distances below threshold
# 	for i in range(len(distance_matrix)):
# 		for j in range(i+1, len(distance_matrix)):
# 			if distance_matrix[i, j] < distance_threshold:
# 				G.add_edge(i, j, weight=1-distance_matrix[i, j])  # Convert distance to similarity
	
# 	# Create plot
# 	fig, ax = plt.subplots(figsize=(12, 12))
	
# 	# Set positions using force-directed layout
# 	pos = nx.spring_layout(G, weight='weight', seed=42)
	
# 	# Draw edges
# 	nx.draw_networkx_edges(
# 		G, pos,
# 		width=0.5,
# 		alpha=0.3,
# 		edge_color='gray',
# 		ax=ax
# 	)
	
# 	# Plot nodes with color by family and marker by cluster
# 	for family in unique_families:
# 		for cluster in unique_clusters:
# 			# Find nodes that belong to both this family and this cluster
# 			node_list = [n for n, attr in G.nodes(data=True) 
# 					   if attr['family'] == family and attr['cluster'] == cluster]
			
# 			if node_list:  # Only plot if there are nodes in this combination
# 				nx.draw_networkx_nodes(
# 					G, pos,
# 					nodelist=node_list,
# 					node_color=family_color_map[family],
# 					node_shape=cluster_marker_map[cluster],
# 					node_size=200,
# 					alpha=0.8,
# 					ax=ax
# 				)
	
# 	# Create legends for families and clusters
# 	family_handles = [plt.Line2D([0], [0], marker='o', color=family_color_map[family], 
# 								 linestyle='', markersize=10, label=family) 
# 					  for family in unique_families]
	
# 	cluster_handles = [plt.Line2D([0], [0], marker=cluster_marker_map[cluster], 
# 								 color='black', linestyle='', markersize=10,
# 								 label=f"Cluster {cluster}") 
# 					  for cluster in unique_clusters]
	
# 	# Place legends
# 	first_legend = ax.legend(handles=family_handles, title="Motif Family", 
# 						   loc='upper left', bbox_to_anchor=(1.05, 1))
# 	ax.add_artist(first_legend)
	
# 	ax.legend(handles=cluster_handles, title="HDBSCAN Cluster", 
# 			 loc='upper left', bbox_to_anchor=(1.05, 0.6))
	
# 	plt.title(f"RNA Motif Network: Family Colors and Cluster Markers (threshold: {distance_threshold})")
# 	plt.axis('off')
# 	plt.tight_layout()
# 	plt.show()

# 4. Circular plot to better visualize family-cluster relationships
def plot_circular_layout_with_family_colors_cluster_markers(motif_families, unique_families, unique_clusters, distance_matrix, labels, cluster_marker_map, family_color_map):
	# Create a graph
	G = nx.Graph()
	
	# Add nodes with family and cluster attributes
	for i in range(len(distance_matrix)):
		G.add_node(i, cluster=int(labels[i]), family=motif_families[i])
	
	# Add edges based on a threshold
	threshold = np.percentile(distance_matrix.ravel(), 10)  # Use lowest 10% of distances
	for i in range(len(distance_matrix)):
		for j in range(i+1, len(distance_matrix)):
			if distance_matrix[i, j] < threshold:
				G.add_edge(i, j, weight=1-distance_matrix[i, j])
	
	# Create plot
	fig, ax = plt.subplots(figsize=(12, 12))
	
	# Use circular layout
	pos = nx.circular_layout(G)
	
	# Draw edges
	nx.draw_networkx_edges(
		G, pos,
		width=0.5,
		alpha=0.2,
		edge_color='gray',
		ax=ax
	)
	
	# Plot nodes with color by family and marker by cluster
	for family in unique_families:
		for cluster in unique_clusters:
			# Find nodes that belong to both this family and this cluster
			node_list = [n for n, attr in G.nodes(data=True) 
					   if attr['family'] == family and attr['cluster'] == cluster]
			
			if node_list:  # Only plot if there are nodes in this combination
				nx.draw_networkx_nodes(
					G, pos,
					nodelist=node_list,
					node_color=family_color_map[family],
					node_shape=cluster_marker_map[cluster],
					node_size=200,
					alpha=0.8,
					ax=ax
				)
	
	# Create legends
	family_handles = [plt.Line2D([0], [0], marker='o', color=family_color_map[family], 
								 linestyle='', markersize=10, label=family) 
					  for family in unique_families]
	
	cluster_handles = [plt.Line2D([0], [0], marker=cluster_marker_map[cluster], 
								 color='black', linestyle='', markersize=10,
								 label=f"Cluster {cluster}") 
					  for cluster in unique_clusters]
	
	# Place legends
	first_legend = ax.legend(handles=family_handles, title="Motif Family", 
						   loc='upper left', bbox_to_anchor=(1.05, 1))
	ax.add_artist(first_legend)
	
	ax.legend(handles=cluster_handles, title="HDBSCAN Cluster", 
			 loc='upper left', bbox_to_anchor=(1.05, 0.6))
	
	plt.title("Circular Layout: RNA Motifs by Family Colors and Cluster Markers")
	plt.axis('off')
	plt.tight_layout()
	plt.show()

# 5. Cluster-Family Comparison Heatmap
def plot_cluster_family_distribution(motif_families, unique_families, unique_clusters, distance_matrix, labels, cluster_marker_map, family_color_map, plot_output_fname, clustering_method, pdb_chain):
	# Create a cross-tabulation of clusters vs families
	cluster_family_counts = np.zeros((len(unique_clusters), len(unique_families)))
	
	for i, cluster in enumerate(unique_clusters):
		for j, family in enumerate(unique_families):
			cluster_family_counts[i, j] = np.sum((labels == cluster) & (motif_families == family))
	
	# Create a heatmap
	fig, ax = plt.subplots(figsize=(12, 8))
	
	# Normalize by row (cluster) to show distribution within each cluster
	row_sums = cluster_family_counts.sum(axis=1, keepdims=True)
	normalized_counts = cluster_family_counts / row_sums
	
	# Plot heatmap
	im = ax.imshow(normalized_counts, cmap='YlOrRd')
	
	# Add labels
	ax.set_xticks(np.arange(len(unique_families)))
	ax.set_yticks(np.arange(len(unique_clusters)))
	ax.set_xticklabels(unique_families)
	ax.set_yticklabels([f"Cluster {c}" if c >= 0 else 'Noise' for c in unique_clusters])
	
	# Rotate x labels for better readability
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
	
	# Add colorbar
	cbar = ax.figure.colorbar(im, ax=ax)
	cbar.set_label('Fraction of cluster members')
	
	# Add text annotations showing the actual counts
	for i in range(len(unique_clusters)):
		for j in range(len(unique_families)):
			count = int(cluster_family_counts[i, j])
			if count > 0:
				text = ax.text(j, i, count, ha="center", va="center", 
							   color="black" if normalized_counts[i, j] < 0.7 else "white")
	
	ax.set_title("Distribution of Motif Families within " + clustering_method + " Clusters (" + pdb_chain + ")")
	fig.tight_layout()
	# plt.show()
	plt.savefig(plot_output_fname)
	plt.close()

# 6. Combined UMAP visualization with detailed annotations
def plot_detailed_umap(motif_families, unique_families, unique_clusters, distance_matrix, labels, cluster_marker_map, family_color_map, plot_output_fname):
	# Reduce dimensions using UMAP
	umap_reducer = UMAP(metric='precomputed', random_state=42)
	umap_embedding = umap_reducer.fit_transform(distance_matrix)
	
	# Create figure
	fig, ax = plt.subplots(figsize=(14, 12))
	
	# Plot each point with color by family and marker by cluster
	for family in unique_families:
		for cluster in unique_clusters:
			# Find points that belong to both this family and this cluster
			mask = (motif_families == family) & (labels == cluster)
			if np.any(mask):  # Only plot if there are points in this combination
				ax.scatter(
					umap_embedding[mask, 0], 
					umap_embedding[mask, 1],
					marker=cluster_marker_map[cluster],
					color=family_color_map[family],
					s=100, 
					alpha=0.8,
					edgecolors='black' if cluster != -1 else 'gray',  # Highlight non-noise points
					linewidth=1 if cluster != -1 else 0.5
				)
	
	# Add annotations for cluster centers
	for cluster in unique_clusters:
		if cluster != -1:  # Skip noise points
			cluster_mask = labels == cluster
			center_x = np.mean(umap_embedding[cluster_mask, 0])
			center_y = np.mean(umap_embedding[cluster_mask, 1])
			
			# Add a text label for the cluster
			ax.text(center_x, center_y, f"C{cluster}", 
					fontsize=16, fontweight='bold', 
					ha='center', va='center',
					bbox=dict(facecolor='white', alpha=0.7, edgecolor='black', boxstyle='round,pad=0.3'))
	
	# Create separate legends for families and clusters
	family_handles = [plt.Line2D([0], [0], marker='o', color=family_color_map[family], 
								 linestyle='', markersize=10, label=family) 
					  for family in unique_families]
	
	cluster_handles = [plt.Line2D([0], [0], marker=cluster_marker_map[cluster], 
								 color='black', linestyle='', markersize=10,
								 # label=f"Cluster {cluster}") 
								label=f'Cluster {cluster}' if cluster >= 0 else 'Noise')
					  for cluster in unique_clusters]
	
	# Place the legends
	first_legend = ax.legend(handles=family_handles, title="Motif Family", 
						   loc='upper left', bbox_to_anchor=(1.05, 1))
	ax.add_artist(first_legend)
	
	ax.legend(handles=cluster_handles, title="HDBSCAN Cluster", 
			 loc='upper left', bbox_to_anchor=(1.05, 0.6))
	
	# Add metrics for cluster-family correspondence
	from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
	
	# Filter out noise points for metrics
	mask = labels != -1
	if sum(mask) > 1:
		ari = adjusted_rand_score(motif_families[mask], labels[mask])
		ami = adjusted_mutual_info_score(motif_families[mask], labels[mask])
		
		# Add metrics to the plot
		plt.figtext(0.02, 0.02, f"Adjusted Rand Index: {ari:.3f}\nAdjusted Mutual Info: {ami:.3f}", 
				   fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
	
	plt.title('Detailed UMAP: RNA Motifs by Family (color) and Cluster (marker)')
	plt.tight_layout()
	plt.subplots_adjust(right=0.8)
	# plt.show()
	plt.savefig(plot_output_fname)
	plt.close()

def get_shortcode(family):
	if family.startswith('IL_') or family.startswith('HL_'):
		return get_known_motif_shortcode(family[3:])
	return family

# def get_shortcoded_list(family_list):
# 	shortcoded_list = []
# 	for fam in family_list:
# 		if fam.startswith('IL_') or fam.startswith('HL_'):
# 			shortcoded_list.append(get_known_motif_shortcode(fam[3:]))
# 		else:
# 			shortcoded_list.append(fam)
# 	return shortcoded_list