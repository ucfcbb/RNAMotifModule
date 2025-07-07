# import matplotlib.pyplot as plt
# import seaborn as sns
import os
import sys
import copy
import logging
import random
import numpy as np
np.random.seed(42)
from Bio.PDB import *
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
from graph_utils import *
from utils import *
from module import *
from logger import *
from tqdm import tqdm
import multiprocessing as mp

from collections import Counter
from scipy.stats import binom, pearsonr
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

import pandas as pd
# from sklearn.preprocessing import MinMaxScaler
pd.set_option('display.max_rows', None)  # Displays all rows
pd.set_option('display.max_columns', None)  # Displays all columns

# def plot_motif_heatmap(family_group_dict_with_ids):
#	 # Convert motif groups to a matrix
#	 motifs = sorted(set(fam_id for module_id, sorted_family_group_tuple in family_group_dict_with_ids.keys() for fam_id in sorted_family_group_tuple))
#	 matrix = [[motif_group_frequencies.get(frozenset([m1, m2]), 0) for m2 in motifs] for m1 in motifs]
	
#	 plt.figure(figsize=(12, 10))
#	 sns.heatmap(matrix, xticklabels=motifs, yticklabels=motifs, annot=True, cmap='YlOrRd')
#	 plt.title('Motif Co-occurrence Frequency')
#	 plt.tight_layout()
#	 plt.show()

def _generate_modules_for_random_data_worker(p):
	generate_modules_for_random_data(*p)

# def generate_modules_for_random_data(round_id, randomized_pdb_chainwise_loops, families, distance_threshold_to_be_nearest_residue, directories):
def generate_modules_for_random_data(round_id, pdb_chainwise_residue_data, pdb_chainwise_loops, families, family_group_dict_with_ids, distance_threshold_to_be_nearest_residue, directories):
	randomized_pdb_chainwise_loops = generate_randomized_dataset(pdb_chainwise_residue_data, pdb_chainwise_loops, directories)
	family_group_dict_for_random_data = generate_motif_module_info_from_pdb_chainwise_data(randomized_pdb_chainwise_loops, families, distance_threshold_to_be_nearest_residue)
	temp_files_dir = os.path.join(directories.temp_dir, 'modules_from_randomized_data')
	create_directory(temp_files_dir)
	family_group_dict_with_ids = dict(sorted(family_group_dict_with_ids.items(), key=lambda item: len(item[1]), reverse=True))
	
	counts = []
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		if sorted_family_group_tuple in family_group_dict_for_random_data:
			counts.append(str(len(family_group_dict_for_random_data[sorted_family_group_tuple])))
		else:
			counts.append('0')

	fp = open(os.path.join(temp_files_dir, 'modules_for_randomized_dataset_' + str(round_id) + '.txt'), 'w')
	fp.write('\n'.join(counts))
	fp.close()

def calculate_p_value_and_expected_freq(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, directories, number_of_multiprocess, randomized_round=100):
	family_group_dict_with_ids = dict(sorted(family_group_dict_with_ids.items(), key=lambda item: len(item[1]), reverse=True))
	p_value_dict = {}
	counter_dict = {}
	expected_freq_dict = {}

	# impl_type_to_test = 'multiprocessing'
	impl_type_to_test = 'serial'
	logger.info('Generating ' + str(randomized_round) + ' randomized dataset(s) to calculate expected frequency and p-value.')
	start_time = time.time()


	##### Multi-processing implementation (not feasible due to parallel writing):
	if impl_type_to_test == 'multiprocessing':
		parameter_list = []
		for i in range(randomized_round):
			# randomized_pdb_chainwise_loops = generate_randomized_dataset(pdb_chainwise_residue_data, pdb_chainwise_loops, directories)
			parameter_list.append((i, pdb_chainwise_residue_data, pdb_chainwise_loops, families, family_group_dict_with_ids, distance_threshold_to_be_nearest_residue, directories))
		pool = mp.Pool(number_of_multiprocess)
		pool.map(_generate_modules_for_random_data_worker, parameter_list)
		wait_for_certain_time_according_to_wait_factor(len(parameter_list))
		temp_files_dir = os.path.join(directories.temp_dir, 'modules_from_randomized_data')
		counts_by_rounds = []
		for i in range(randomized_round):
			fp = open(os.path.join(temp_files_dir, 'modules_for_randomized_dataset_' + str(i) + '.txt'))
			lines = fp.readlines()
			fp.close()
			counts = []
			for line in lines:
				counts.append(line.strip())
			counts_by_rounds.append(counts)

		j = 0
		for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
			if sorted_family_group_tuple not in counter_dict:
				counter_dict[sorted_family_group_tuple] = []
			for i in range(randomized_round):
				counter_dict[sorted_family_group_tuple].append(int(counts_by_rounds[i][j]))
			j += 1
	##### Multi-processing implementation ends

	else:
	##### Serial implementation:
		for i in range(randomized_round):
			print('round: ' + str(i+1))
			logger.info('Step1: Generate random dataset.')
			randomized_pdb_chainwise_loops = generate_randomized_dataset(pdb_chainwise_residue_data, pdb_chainwise_loops, directories)
			logger.info('Step2: Identify modules in the random dataset.')
			family_group_dict_for_random_data = generate_motif_module_info_from_pdb_chainwise_data(randomized_pdb_chainwise_loops, families, distance_threshold_to_be_nearest_residue)
			
			# print(family_group_dict_for_random_data.keys())

			for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
				if sorted_family_group_tuple in family_group_dict_for_random_data:
					if sorted_family_group_tuple not in counter_dict:
						counter_dict[sorted_family_group_tuple] = []
					counter_dict[sorted_family_group_tuple].append(len(family_group_dict_for_random_data[sorted_family_group_tuple]))
				else:
					if sorted_family_group_tuple not in counter_dict:
						counter_dict[sorted_family_group_tuple] = []
					counter_dict[sorted_family_group_tuple].append(0)
	##### Serial implementation ends

	logger.info('Calculating expected frequency and p-value.')
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		if sorted_family_group_tuple in counter_dict:
			p_value_dict[sorted_family_group_tuple] = sum(1 for val in counter_dict[sorted_family_group_tuple] if val >= len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])) / randomized_round
			expected_freq_dict[sorted_family_group_tuple] = sum(counter_dict[sorted_family_group_tuple]) / randomized_round

	logger.info('Done')
	logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

	return p_value_dict, expected_freq_dict, counter_dict

# HL/IL merged
# def generate_randomized_families(families):
# 	randomized_families = {}
# 	family_list = []
# 	loops_list = []
# 	family_loop_count = {}
# 	for family_id in families:
# 		family_list.append(family_id)
# 		randomized_families[family_id] = []
# 		loops = families[family_id]
# 		family_loop_count[family_id] = len(loops)
# 		loops_list += loops

# 	for loop in loops_list:
# 		while True:
# 			random_ind = random.randint(0, len(family_list)-1)
# 			random_family_id = family_list[random_ind]
# 			if len(randomized_families[random_family_id]) < family_loop_count[random_family_id]:
# 				randomized_families[random_family_id].append(loop)
# 				break

# 	return randomized_families

def generate_randomized_families(families):
	randomized_families = {}
	junctionwise_family_list = {}
	loops_list = []
	family_loop_count = {}
	for family_id in families:
		loops = families[family_id]
		if len(loops) == 0:
			continue
		junc_cnt = get_junction_count(loops[0])
		if junc_cnt not in junctionwise_family_list:
			junctionwise_family_list[junc_cnt] = []
		junctionwise_family_list[junc_cnt].append(family_id)

		# family_list.append(family_id)
		randomized_families[family_id] = []
		
		family_loop_count[family_id] = len(loops)
		loops_list += loops

	for loop in loops_list:
		junc_cnt = get_junction_count(loop)
		while True:
			random_ind = random.randint(0, len(junctionwise_family_list[junc_cnt])-1)
			random_family_id = junctionwise_family_list[junc_cnt][random_ind]
			if len(randomized_families[random_family_id]) < family_loop_count[random_family_id]:
				randomized_families[random_family_id].append(loop)
				break

	return randomized_families

def update_coord_data_in_pdb_chainwise_loops(pdb_chainwise_loops, pdb_chainwise_residue_data):
	new_pdb_chainwise_loops = {}
	for i, pdb_id in enumerate(pdb_chainwise_loops):
		new_pdb_chainwise_loops[pdb_id] = {}
		for j, chain_id in enumerate(pdb_chainwise_loops[pdb_id]):
			new_pdb_chainwise_loops[pdb_id][chain_id] = []
			loop_info_dict, min_coords, max_coords = pdb_chainwise_residue_data[pdb_id][chain_id]
			for loop in loop_info_dict:
				residuewise_motif_coords = loop_info_dict[loop]
				new_pdb_chainwise_loops[pdb_id][chain_id].append((loop, residuewise_motif_coords))

	return new_pdb_chainwise_loops

def calculate_p_value_and_expected_freq_v2a(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, directories, output_dir, write_to_file, number_of_multiprocess, randomized_round=100):
	family_group_dict_with_ids = dict(sorted(family_group_dict_with_ids.items(), key=lambda item: len(item[1]), reverse=True))
	p_value_dict = {}
	counter_dict = {}
	expected_freq_dict = {}

	# impl_type_to_test = 'multiprocessing'
	# impl_type_to_test = 'serial'
	print('')
	logger.info('Generating ' + str(randomized_round) + ' randomized dataset(s) to calculate expected frequency and p-value.')
	start_time = time.time()


	for i in range(randomized_round):
		print('')
		print('Round: ' + str(i+1))
		logger.info('Step1: Generate random family dataset.')
		randomized_families = generate_randomized_families(families)
		# print(randomized_families)
		# for family_id in randomized_families:
		# 	print(family_id + ',' + ','.join(randomized_families[family_id]))
		# sys.exit()
		logger.info('Step2: Identify modules in the random family dataset.')
		# print(pdb_chainwise_loops)
		# print(pdb_chainwise_residue_data)
		# sys.exit()
		pdb_chainwise_loops = update_coord_data_in_pdb_chainwise_loops(pdb_chainwise_loops, pdb_chainwise_residue_data)
		# print(pdb_chainwise_loops)
		# sys.exit()
		family_group_dict_for_random_data = generate_motif_module_info_from_pdb_chainwise_data(pdb_chainwise_loops, randomized_families, distance_threshold_to_be_nearest_residue)

		for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
			if sorted_family_group_tuple in family_group_dict_for_random_data:
				if sorted_family_group_tuple not in counter_dict:
					counter_dict[sorted_family_group_tuple] = []
				counter_dict[sorted_family_group_tuple].append(len(family_group_dict_for_random_data[sorted_family_group_tuple]))
			else:
				if sorted_family_group_tuple not in counter_dict:
					counter_dict[sorted_family_group_tuple] = []
				counter_dict[sorted_family_group_tuple].append(0)

	logger.info('Calculating expected frequency and p-value.')
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		if sorted_family_group_tuple in counter_dict:
			p_value_dict[sorted_family_group_tuple] = sum(1 for val in counter_dict[sorted_family_group_tuple] if val >= len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])) / randomized_round
			expected_freq_dict[sorted_family_group_tuple] = sum(counter_dict[sorted_family_group_tuple]) / randomized_round

	logger.info('Done')
	logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

	plots_dir = os.path.join(output_dir, 'plots')
	stats_dir = os.path.join(output_dir, 'stats')
	create_directory(plots_dir)
	create_directory(stats_dir)

	plot_output_fname = os.path.join(plots_dir, 'stats_on_rand_data_bubble_plot.png')
	plot_title = 'Motif Modules: Counts, p-values, and Observed/Expected Ratios based on randomized data'
	create_bubble_plot(family_group_dict_with_ids, expected_freq_dict, p_value_dict, plot_output_fname, plot_title)


	# Display and Write to file
	print('')
	print('Motif Module statistics based on randomized data')
	print('Module\tObservedFreq\tExpectedFreq\tp-value\tz-score')

	output_fname = os.path.join(stats_dir, 'motif_module_stats.txt')
	if write_to_file:
		fp = open(output_fname, 'w')
		fp.write('Module\tObservedFreq\tExpectedFreq\tp-value\tz-score\n')
		fp.close()

	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		observed_freq = len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])
		# expected_freq = expected_freq_dict[(module_id, sorted_family_group_tuple)]
		# chi2, p_value = chi_square_test(observed_freq, expected_freq)
		# chi2, p_value = chi_square_test(observed_freq, expected_freq, total_chains)
		# manual_chi2 = manual_chi_square([observed_freq], [expected_freq])
		# _, fisher_p_value = fisher_exact([[observed_freq, total_chains - observed_freq], [expected_freq, total_chains - expected_freq]])

		p_value_r = 0.0
		expected_freq_r = 0.0
		freq_list = counter_dict[sorted_family_group_tuple]

		mean_randomized = np.mean(freq_list)
		std_randomized = np.std(freq_list)
		z_score = 0.0
		if std_randomized > 0.0:
			z_score = (observed_freq - mean_randomized) / std_randomized
		if sorted_family_group_tuple in p_value_dict:
			p_value_r = p_value_dict[sorted_family_group_tuple]
			expected_freq_r = expected_freq_dict[sorted_family_group_tuple]
			# freq_list = counter_dict[sorted_family_group_tuple]

		# print(str(sorted_family_group_tuple) + '\t' + str(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])) + '\t' + str(expected_freq_dict[(module_id, sorted_family_group_tuple)]) + '\t' + str(chi2) + '\t' + str(p_value))
		# print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq:.2f}\t{observed_freq/expected_freq:.4f}\t{chi2:.2f}\t{manual_chi2:.2f}\t{p_value:.4f}\t{fisher_p_value:.4f}")
		# print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq:.2f}\t{observed_freq/expected_freq:.4f}\t{manual_chi2:.2f}\t{fisher_p_value:.4f}\t{p_value_r:.4f}\t{expected_freq_r:.4f}\t{z_score:.4f}\t{freq_list}")
		# print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq_r:.2f}\t{p_value_r:.4f}\t{z_score:.4f}\t{freq_list}")
		print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq_r:.2f}\t{p_value_r:.4f}\t{z_score:.4f}")
		if write_to_file:
			p_value_text = str(round(p_value_r, 3))
			if p_value_r == 0.0:
				p_value_text = '<' + str(1/randomized_round)
			fp = open(output_fname, 'a')
			fp.write(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq_r:.2f}\t{p_value_text}\t{z_score:.4f}\n")
			fp.close()

	# return p_value_dict, expected_freq_dict, counter_dict

def calculate_enrichment(families, family_group_dict_with_ids, output_dir, write_to_file):
	print('')
	logger.info('Analyzing the pairwise enrichment between motif families in modules.')

	all_families = list(families.keys())
	cooccurrence = np.zeros((len(all_families), len(all_families)))
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		for i, fam1 in enumerate(all_families):
			for j, fam2 in enumerate(all_families):
				if fam1 in sorted_family_group_tuple and fam2 in sorted_family_group_tuple:
					cooccurrence[i, j] += 1
	# print(cooccurrence)
	# sys.exit()
	odds_ratio_matrix = np.zeros((len(all_families), len(all_families)))
	p_value_matrix = np.ones((len(all_families), len(all_families)))

	# Create your odds ratio DataFrame
	odds_ratio_data = pd.DataFrame(index=all_families, columns=all_families)
	significance_matrix = pd.DataFrame(False, index=all_families, columns=all_families)


	enrichment_results = []
	for i, fam1 in enumerate(all_families):
		for j, fam2 in enumerate(all_families):
			if i >= j:  # Skip duplicates and self-comparisons
				continue
				
			# Create contingency table
			a = cooccurrence[i, j]  # Modules with both fam1 and fam2
			b = sum(1 for module_id, sorted_family_group_tuple in family_group_dict_with_ids if 
					any(f == fam1 for f in sorted_family_group_tuple) and 
					not any(f == fam2 for f in sorted_family_group_tuple))
			c = sum(1 for module_id, sorted_family_group_tuple in family_group_dict_with_ids if 
					not any(f == fam1 for f in sorted_family_group_tuple) and 
					any(f == fam2 for f in sorted_family_group_tuple))
			d = sum(1 for module_id, sorted_family_group_tuple in family_group_dict_with_ids if 
					not any(f == fam1 for f in sorted_family_group_tuple) and 
					not any(f == fam2 for f in sorted_family_group_tuple))
			
			# Fisher's exact test
			odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
			# odds_ratio_matrix[i, j] = odds_ratio_matrix[j, i] = odds_ratio
			# p_value_matrix[i, j] = p_value_matrix[j, i] = p_value

			odds_ratio_data.loc[fam1, fam2] = odds_ratio
			odds_ratio_data.loc[fam2, fam1] = odds_ratio
			if p_value < 0.05:
				significance_matrix.loc[fam1, fam2] = True
				significance_matrix.loc[fam2, fam1] = True
			
			enrichment_results.append({
				"Family1": fam1,
				"Family2": fam2,
				"Modules_with_both": a,
				"Modules_with_only_fam1": b,
				"Modules_with_only_fam2": c,
				"Modules_with_neither": d,
				"Odds_ratio": odds_ratio,
				"P_value": p_value
			})
	# print(enrichment_results)
	# sys.exit()

	p_values = [result["P_value"] for result in enrichment_results]
	# corrected_p = multipletests(p_values, method='fdr_bh')[1]  # Benjamini-Hochberg

	# print(odds_ratio_data)
	# print('odds_ratio_data.dtypes')
	# print(odds_ratio_data.dtypes)
	# odds_ratio_data = odds_ratio_data.astype(float, errors='coerce')
	# print(odds_ratio_data.dtypes)

	shortcode_dict = {}
	for fam in all_families:
		if fam.startswith('IL_') or fam.startswith('HL_'):
			shortcode_dict[fam] = get_known_motif_shortcode(fam[3:])

	odds_ratio_data = odds_ratio_data.rename(columns=shortcode_dict)
	odds_ratio_data = odds_ratio_data.rename(index=shortcode_dict)
	significance_matrix = significance_matrix.rename(columns=shortcode_dict)
	significance_matrix = significance_matrix.rename(index=shortcode_dict)


	plots_dir = os.path.join(output_dir, 'plots')
	stats_dir = os.path.join(output_dir, 'stats')
	create_directory(plots_dir)
	create_directory(stats_dir)
	plot_output_fname = os.path.join(plots_dir, 'pairwise_enrichment_between_motif_families_in_modules.png')

	create_detailed_enrichment_heatmap(
		odds_ratio_data,
		significance_matrix,
		plot_output_fname,
		title="Pairwise RNA Motif Module Enrichment",
		x_label="Odds Ratio",
		vmin=0.1,
		vmax=9.0,
		significance_symbol="*",
		triangle="upper"
	)

	
	# sys.exit()

	# Write in file:
	if write_to_file:
		output_fname = os.path.join(stats_dir, 'pairwise_enrichment_between_motif_families_in_modules.tsv')
		fp = open(output_fname, 'w')
		fp.write('Family1\tFamily2\tOdds_ratio\tP_value\n')
		for d in enrichment_results:
			# fp.write(d["Family1"] + "\t" + d["Family2"] + "\t" + str(d["Odds_ratio"]) + "\t" + str(d["P_value"]) + "\n")
			fp.write(f"{d['Family1']}\t{d['Family2']}\t{d['Odds_ratio']:.4f}\t{d['P_value']:.4f}\n")
		fp.close()


	# Display results
	print('')
	# results_df = pd.DataFrame(enrichment_results)
	# print("Family1, Family2, Odds_ratio, P_value, Corrected_P_value, Significant")
	print('Pairwise Enrichment between Motif Families in Modules')
	print("Family1, Family2, Odds_ratio, P_value")
	for d in enrichment_results:
		# print(d["Family1"] + ", " + d["Family2"] + ", " + str(d["Odds_ratio"]) + ", " + str(d["P_value"]) + ", " + str(d["Corrected_P_value"]) + ", " + str(d["Significant"]))
		print(d["Family1"] + ", " + d["Family2"] + ", " + str(d["Odds_ratio"]) + ", " + str(d["P_value"]))
	# print(results_df[["Family1", "Family2", "Odds_ratio", "P_value", "Corrected_P_value", "Significant"]])
	# # sys.exit()

	# mask = np.zeros_like(odds_ratio_data, dtype=bool)
	# mask[np.triu_indices_from(mask)] = True

	print('')
	logger.info('Analyzing specific module compositions (e.g., modules with exactly 2 KT, 1 SR, 1 CL)')

	compositions = []
	composition_counts = Counter()
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		comp = Counter()
		for family in sorted_family_group_tuple:
			comp[family] += 1
		compositions.append(dict(comp))
		composition_counts[tuple(sorted(comp.items()))] = len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])

	# print(compositions)
	# sys.exit()

	# composition_counts = Counter(tuple(sorted(comp.items())) for comp in compositions)
	# print(composition_counts)
	# sys.exit()

	all_motifs = Counter()
	# for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
	for family in families:
		all_motifs[family] += len(families[family])
	total_motifs = sum(all_motifs.values())
	family_probs = {f: c/total_motifs for f, c in all_motifs.items()}

	# target_comp = Counter({"KT": 2, "SR": 1, "CL": 1})
	# observed_count = sum(1 for comp in compositions if Counter(comp) == target_comp)
	total_modules = sum(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)]) for module_id, sorted_family_group_tuple in family_group_dict_with_ids)
	# print(total_motifs)
	# sys.exit()
	expected_prob_dict = {}
	expected_count_dict = {}
	p_value_dict = {}
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		observed_count = len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])
		expected_prob = 1
		for family in sorted_family_group_tuple:
			expected_prob *= family_probs[family]

		expected_prob_dict[sorted_family_group_tuple] = expected_prob
		expected_count_dict[sorted_family_group_tuple] = round(expected_prob * len(family_group_dict_with_ids), 0)

		p_value = 1 - binom.cdf(observed_count-1, total_modules, expected_prob)
		p_value_dict[sorted_family_group_tuple] = p_value

		# print(module_id, sorted_family_group_tuple, observed_count, p_value)

	
	# expected_prob = 0.05  # This would be calculated based on family probabilities
	# p_value = 1 - binom.cdf(observed_count-1, total_modules, expected_prob)

	# print(f"Composition {dict(target_comp)}: Observed={observed_count}, p-value={p_value:.5f}")
	# sys.exit()
	plots_dir = os.path.join(output_dir, 'plots')
	stats_dir = os.path.join(output_dir, 'stats')
	create_directory(plots_dir)
	create_directory(stats_dir)
	plot_output_fname = os.path.join(plots_dir, 'specific_module_composition_analysis.png')
	plot_title = 'Motif Modules: Counts, p-values, and Observed/Expected Ratios based on enrichment analysis'
	create_bubble_plot(family_group_dict_with_ids, expected_count_dict, p_value_dict, plot_output_fname, plot_title)

	# Write in file:
	if write_to_file:
		output_fname = os.path.join(stats_dir, 'specific_module_composition_analysis.tsv')
		fp = open(output_fname, 'w')
		fp.write('Module\tObservedFreq\tExpectedProb\tExpectedFreq\tp-value\n')
		for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
			observed_freq = len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])
			expected_prob = expected_prob_dict[sorted_family_group_tuple]
			expected_count = expected_count_dict[sorted_family_group_tuple]
			p_value = p_value_dict[sorted_family_group_tuple]
			fp.write(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_prob:.2f}\t{expected_count:.0f}\t{p_value:.4f}\n")
		fp.close()

	# display results:
	print('Specific Module Composition Analysis')
	print('Module, ObservedFreq, ExpectedProb, ExpectedFreq, p-value')
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		observed_freq = len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])
		expected_prob = expected_prob_dict[sorted_family_group_tuple]
		expected_count = expected_count_dict[sorted_family_group_tuple]
		p_value = p_value_dict[sorted_family_group_tuple]
		print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_prob:.2f}\t{expected_count:.0f}\t{p_value:.4f}")

def family_count_correlation(families, family_group_dict_with_ids, output_dir, write_to_file):
	print('')
	logger.info('Analyzing if certain families tend to have multiple instances together.')

	all_families = list(families.keys())
	total_modules = sum(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)]) for module_id, sorted_family_group_tuple in family_group_dict_with_ids)
	count_matrix = np.zeros((total_modules, len(all_families)))
	print(total_modules)

	i = 0
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		for k in range(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])):
			for j, family in enumerate(all_families):
				count_matrix[i, j] = sum(1 for fam in sorted_family_group_tuple if fam == family)

			i += 1

	# print(count_matrix[0])
	# sys.exit()
	# Create count matrix (rows=modules, columns=motif families)
	# count_matrix = np.zeros((len(modules), len(all_families)))
	# for i, module in enumerate(modules):
	# 	for j, family in enumerate(all_families):
	# 		count_matrix[i, j] = sum(count for fam, count in module["motifs"] if fam == family)

	# Calculate count correlations
	# corr_matrix = np.zeros((len(all_families), len(all_families)))
	corr_matrix = pd.DataFrame(index=all_families, columns=all_families)
	# pval_matrix = np.zeros((len(all_families), len(all_families)))
	pval_matrix = pd.DataFrame(index=all_families, columns=all_families)
	# sign_matrix = np.zeros((len(all_families), len(all_families)))
	sign_matrix = pd.DataFrame(False, index=all_families, columns=all_families)

	correlation_results = []
	for i in range(len(all_families)):
		fam1 = all_families[i]
		for j in range(i+1, len(all_families)):
			fam2 = all_families[j]
			corr, pval = pearsonr(count_matrix[:, i], count_matrix[:, j])
			# corr_matrix[i, j] = corr_matrix[j, i] = corr
			corr_matrix.loc[fam1, fam2] = corr_matrix.loc[fam2, fam1] = corr
			# pval_matrix[i, j] = pval_matrix[j, i] = pval
			pval_matrix.loc[fam1, fam2] = pval_matrix.loc[fam2, fam1] = pval
			if pval < 0.05:
				sign_matrix.loc[fam1, fam2] = sign_matrix.loc[fam2, fam1] = True

			correlation_results.append({
				"Family1": fam1,
				"Family2": fam2,
				"Pearson_corr_coeff": corr,
				"P_value": pval
			})

	# Display correlation results
	# corr_df = pd.DataFrame(corr_matrix, index=all_families, columns=all_families)
	# print("Count correlations between motif families:")
	# print(corr_df)

	# pval_df = pd.DataFrame(pval_matrix, index=all_families, columns=all_families)
	# print("Count p-values between motif families:")
	# print(pval_df)
	# sys.exit()

	shortcode_dict = {}
	for fam in all_families:
		if fam.startswith('IL_') or fam.startswith('HL_'):
			shortcode_dict[fam] = get_known_motif_shortcode(fam[3:])

	corr_matrix = corr_matrix.rename(columns=shortcode_dict)
	corr_matrix = corr_matrix.rename(index=shortcode_dict)
	sign_matrix = sign_matrix.rename(columns=shortcode_dict)
	sign_matrix = sign_matrix.rename(index=shortcode_dict)

	plots_dir = os.path.join(output_dir, 'plots')
	stats_dir = os.path.join(output_dir, 'stats')
	create_directory(plots_dir)
	create_directory(stats_dir)
	plot_output_fname = os.path.join(plots_dir, 'family_count_correlation_analysis.png')

	create_detailed_enrichment_heatmap(
		corr_matrix,
		sign_matrix,
		plot_output_fname,
		title="Family Count Correlation Analysis",
		x_label="Pearson correlation coefficient",
		cmap='RdBu_r',
		vmin=-1.0,
		vmax=1.0,
		center=0.0,
		significance_symbol="*",
		triangle="upper"
	)

	# Write in file:
	if write_to_file:
		output_fname = os.path.join(stats_dir, 'family_count_correlation_analysis.tsv')
		fp = open(output_fname, 'w')
		fp.write('Family1\tFamily2\tPearson_corr_coeff\tP_value\n')
		for d in correlation_results:
			# fp.write(d["Family1"] + "\t" + d["Family2"] + "\t" + str(d["Pearson_corr_coeff"]) + "\t" + str(d["P_value"]) + "\n")
			fp.write(f"{d['Family1']}\t{d['Family2']}\t{d['Pearson_corr_coeff']:.4f}\t{d['P_value']:.4f}\n")
		fp.close()

	print('')
	print('Family Count Correlation Analysis')
	print("Family1, Family2, tPearson_corr_coeff, P_value")
	for d in correlation_results:
		# print(d["Family1"] + ", " + d["Family2"] + ", " + str(d["Odds_ratio"]) + ", " + str(d["P_value"]) + ", " + str(d["Corrected_P_value"]) + ", " + str(d["Significant"]))
		print(f"{d['Family1']}\t{d['Family2']}\t{d['Pearson_corr_coeff']:.4f}\t{d['P_value']:.4f}")

def calculate_p_value_and_expected_freq_v2b(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, directories, output_dir, write_to_file, number_of_multiprocess, randomized_round=100):
	family_group_dict_with_ids = dict(sorted(family_group_dict_with_ids.items(), key=lambda item: len(item[1]), reverse=True))

	# print(list(families.keys()))
	# sys.exit()

	# Background-Corrected Enrichment Analysis

	# Specific Motif Count Enrichment
	#1. analyze the pairwise enrichment between motif families in modules - MAYBE NOT KEEP IT
	#2. For analyzing specific combinations (e.g., modules with exactly 2 KT, 1 SR, 1 CL)- MAYBE KEEP IT
	calculate_enrichment(families, family_group_dict_with_ids, output_dir, write_to_file)
	# Family Count Correlation Analysis
	# To analyze if certain families tend to have multiple instances together - MAYBE NOT KEEP IT - THINK
	family_count_correlation(families, family_group_dict_with_ids, output_dir, write_to_file)

def get_shortcoded_tuple(family_group_tuple):
	family_group = list(family_group_tuple)
	shortcoded_list = []
	for fam in family_group:
		if fam.startswith('IL_') or fam.startswith('HL_'):
			shortcoded_list.append(get_known_motif_shortcode(fam[3:]))
		else:
			shortcoded_list.append(fam)
	return tuple(shortcoded_list)

def get_shorter_tuple_with_counter(family_group_tuple):
	family_group = list(family_group_tuple)
	shortcoded_numbered_list = []
	counter = Counter()
	for fam in family_group:
		counter[fam] += 1
	counter = dict(counter)
	for fam in sorted(counter):
		shortcoded_numbered_list.append(str(counter[fam]) + fam)
	return tuple(shortcoded_numbered_list)

def create_bubble_plot(family_group_dict_with_ids, expected_freq_dict_r, p_value_dict, plot_output_fname, plot_title):
	module_types = []
	counts = []
	p_values = []
	obs_exp_ratios = []

	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		observed_freq = len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])
		p_value_r = 0.0
		expected_freq_r = 0.0
		if sorted_family_group_tuple in p_value_dict:
			p_value_r = p_value_dict[sorted_family_group_tuple]
			expected_freq_r = expected_freq_dict_r[sorted_family_group_tuple]

		if observed_freq > 1:
			module_types.append(str(get_shorter_tuple_with_counter(get_shortcoded_tuple(sorted_family_group_tuple))))
			counts.append(observed_freq)
			p_values.append(p_value_r)
			if expected_freq_r > 0.0:
				obs_exp_ratios.append(observed_freq/expected_freq_r)
			else:
				obs_exp_ratios.append(np.nan)

	# max_obs_exp_ratio = max(obs_exp_ratios)
	max_obs_exp_ratio = max(r for r in obs_exp_ratios if not np.isnan(r))

	if max_obs_exp_ratio > 3.0:
		max_obs_exp_ratio = 3.0
	# print('Nan items')
	# print(np.sum(np.isnan(obs_exp_ratios)))
	obs_exp_ratios = list(map(lambda x: (max_obs_exp_ratio + 1.5) if np.isnan(x) else 3.0 if x > 3.0 else x, obs_exp_ratios))

	# Create the bubble plot
	fig, ax = plt.subplots(figsize=(12, 6))


	# for i, (mod, count, p, ratio) in enumerate(zip(module_types, counts, p_values, obs_exp_ratios)):
	# 	if not np.isnan(ratio):
	# 		# Regular bubble - size based on obs/exp ratio
	# 		scatter = ax.scatter([i], [count], 
	# 						   s=[ratio*100],  # Size based on ratio, scaled up for visibility
	# 						   c=[p],  # Color based on p-value
	# 						   cmap='viridis_r',  # Reverse colormap so darker = more significant
	# 						   alpha=0.7)
	# 	else:
	# 		# Special case bubble - slightly larger with red border
	# 		special_size = (max_obs_exp_ratio + 1) * 100  # Make them slightly larger than the max
	# 		ax.scatter([i], [count], s=[special_size], 
	# 				 c=[p],  # Fill color
	# 				 cmap='viridis_r',
	# 				 edgecolor='red',  # Border color
	# 				 linewidth=2,  # Border width
	# 				 alpha=0.7)



	# Scatter plot with varying bubble size
	scatter = ax.scatter(module_types, counts, 
						s=[r*100 for r in obs_exp_ratios],  # Size based on ratio, scaled up for visibility
						c=p_values,  # Color based on p-value
						cmap='viridis_r',  # Reverse colormap so darker = more significant
						alpha=0.7)


	# # Now add hatched patches for special cases
	# for i, (count, ratio) in enumerate(zip(counts, obs_exp_ratios)):
	# 	if ratio == max_obs_exp_ratio + 1.5:  # This is a special case
	# 		special_size = np.sqrt((max_obs_exp_ratio + 1.5) * 0.25 / np.pi)  # Convert to radius (circle patch uses radius, not area)
			
	# 		# Create a Circle patch with hatching
	# 		circle = plt.Circle(
	# 			(i, count),				# (x, y) center coordinates
	# 			radius=special_size,	   # Circle radius 
	# 			# facecolor='lightblue',	 # Fill color
	# 			fill=False,
	# 			edgecolor='red',		   # Border color
	# 			linewidth=2,			   # Border width
	# 			alpha=0.7,				 # Transparency
	# 			hatch='///'			   # Hatching pattern
	# 		)
			
	# 		# Add the patch to the plot
	# 		ax.add_patch(circle)

	# Now add stars for special cases (expected=0)
	special_points = [(i, count) for i, (count, ratio) in 
					 enumerate(zip(counts, obs_exp_ratios)) if ratio == max_obs_exp_ratio + 1.5]

	if special_points:
		special_x, special_y = zip(*special_points)
		ax.scatter(special_x, special_y,
				 marker='*',			 # Star marker
				 s=100,				  # Size of the star
				 color='lightgray',			# Color of the star
				 edgecolor='black',	  # Border color
				 linewidth=0.5)			# Border width



	# Customize the plot
	ax.set_xlabel('Motif Modules')
	ax.set_ylabel('Count')
	ax.set_title(plot_title)

	# Add colorbar for p-values
	cbar = plt.colorbar(scatter)
	cbar.set_label('p-value')

	# # Add significance markers
	# for i, p in enumerate(p_values):
	# 	if p <= 0.001:
	# 		sig = '***'
	# 	elif p <= 0.01:
	# 		sig = '**'
	# 	elif p <= 0.05:
	# 		sig = '*'
	# 	else:
	# 		sig = ''
		
	# 	if sig:  # Only add markers for significant values
	# 		ax.annotate(sig, (i, counts[i]), 
	# 				   xytext=(0, 5), textcoords='offset points',
	# 				   ha='center')

	# Add legend for bubble size
	sizes = [0.5, 1.0, 1.5, 2.0, 3.0]#, max_obs_exp_ratio + 1.5]
	labels = ['0.5', '1.0', '1.5', '2.0', '>=3.0']#, 'Expected=0']
	legend_bubbles = []
	for size in sizes:
		legend_bubbles.append(plt.scatter([], [], s=size*100, color='gray', alpha=0.7))

	# # Add special case to legend
	# star_marker = plt.scatter([], [], marker='*', s=(max_obs_exp_ratio + 1.5)*100, color='lightgray', edgecolor='black', linewidth=0.5, alpha=0.7)
	# legend_bubbles.append(star_marker)
	# labels.append('Expected = 0')

	# Add single star marker for special case
	star_marker = plt.scatter([], [], marker='*', s=300, color='lightgray', edgecolor='black', linewidth=0.5)
	legend_bubbles.append(star_marker)
	labels.append('Expected=0')

	# # Create a special entry for the legend with both bubble and star
	# # First create a bubble
	# ax.scatter([0], [0], s=(max_obs_exp_ratio + 1.5)*100, color='gray', alpha=0.5)  # Off-screen placeholder
	# # Then create a star on top (this will be our legend entry)
	# special_marker = ax.scatter([0], [0], s=(max_obs_exp_ratio + 1.5)*100, color='gray', alpha=0.5, label='_nolegend_')
	# star_marker = ax.scatter([0], [0], marker='*', s=100, color='lightgray', edgecolor='black', linewidth=0.5)
	
	# # For the legend, we need to handle this separately
	# from matplotlib.lines import Line2D
	# special_entry = Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
	# 					markersize=np.sqrt((max_obs_exp_ratio + 1.5)*2), alpha=0.5, label='_nolegend_')
	# star_entry = Line2D([0], [0], marker='*', color='w', markerfacecolor='lightgray',
	# 				  markersize=12, markeredgecolor='black', markeredgewidth=0.5)
	
	# # Combine both elements into one legend entry
	# import matplotlib.container as container
	# combined_entry = container.Container([special_entry, star_entry], label='Expected=0')
	# legend_bubbles.append(combined_entry)

	# # Add special case to legend
	# special_bubble = plt.scatter([], [], s=(max_obs_exp_ratio+1)*100, 
	#							facecolor='gray', edgecolor='red', 
	#							linewidth=2, alpha=0.7)
	# legend_bubbles.append(special_bubble)
	# labels.append('Expected = 0')

	# ax.legend(legend_bubbles, labels, scatterpoints=1, labelspacing=1.5, borderpad=1.5,
	# 		  title='Observed/Expected Ratio', loc='upper right')
	ax.legend(legend_bubbles, labels, scatterpoints=1, labelspacing=1.25, borderpad=1.15,
			  title='Observed/Expected Ratio', loc='upper right')

	# plt.tight_layout()
	# Rotate labels and adjust alignment
	ax.set_xticklabels(module_types, rotation=90, ha='center', va='top')

	# Add some padding at the bottom to ensure labels are visible
	plt.tight_layout(pad=3.0)
	# plt.show()

	# plots_dir = os.path.join(output_dir, 'plots')
	# create_directory(plots_dir)
	# plot_output_fname = os.path.join(plots_dir, 'stats_on_rand_data_bubble_plot.png')
	plt.savefig(plot_output_fname)
	# sys.exit()

def prepare_stats(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, total_chains, directories, output_dir, number_of_multiprocess, randomized_round, write_to_file=True):
	# output_fname = os.path.join(output_dir, 'motif_module_stats.txt')

	print('')
	logger.info('Calculating statistical significance ...')
	# print('stats start')
	# p_value_dict, expected_freq_dict_r, counter_dict = calculate_p_value_and_expected_freq(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, directories, number_of_multiprocess, randomized_round)

	# 
	# p_value_dict, expected_freq_dict_r, counter_dict = calculate_p_value_and_expected_freq_v2a(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, directories, output_dir, number_of_multiprocess, randomized_round)
	calculate_p_value_and_expected_freq_v2a(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, directories, output_dir, write_to_file, number_of_multiprocess, randomized_round)

	#
	calculate_p_value_and_expected_freq_v2b(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, directories, output_dir, write_to_file, number_of_multiprocess, randomized_round)

	# sys.exit()
	# expected_freq_dict = calculate_expected_frequency(family_group_dict_with_ids, families, total_chains)

	# print(chi2_contingency())

	# print('Module\tObservedFreq\tExpectedFreq\tRatio\tchi2\tmanual_chi2\tp_value\tfisher_p_value\tp_value_r\texpected_freq_r')
	# print('Module\tObservedFreq\tExpectedFreq\tRatio\tmanual_chi2\tfisher_p_value\tp_value_r\texpected_freq_r\tz_score')
	# print('Module\tObservedFreq\tExpectedFreq\tp-value\tz-score')
	# if write_to_file:
	# 	fp = open(output_fname, 'w')
	# 	fp.write('Module\tObservedFreq\tExpectedFreq\tp-value\tz-score\n')
	# 	fp.close()

	# for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
	# 	observed_freq = len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])
	# 	# expected_freq = expected_freq_dict[(module_id, sorted_family_group_tuple)]
	# 	# chi2, p_value = chi_square_test(observed_freq, expected_freq)
	# 	# chi2, p_value = chi_square_test(observed_freq, expected_freq, total_chains)
	# 	# manual_chi2 = manual_chi_square([observed_freq], [expected_freq])
	# 	# _, fisher_p_value = fisher_exact([[observed_freq, total_chains - observed_freq], [expected_freq, total_chains - expected_freq]])

	# 	p_value_r = 0.0
	# 	expected_freq_r = 0.0
	# 	freq_list = counter_dict[sorted_family_group_tuple]

	# 	mean_randomized = np.mean(freq_list)
	# 	std_randomized = np.std(freq_list)
	# 	z_score = 0.0
	# 	if std_randomized > 0.0:
	# 		z_score = (observed_freq - mean_randomized) / std_randomized
	# 	if sorted_family_group_tuple in p_value_dict:
	# 		p_value_r = p_value_dict[sorted_family_group_tuple]
	# 		expected_freq_r = expected_freq_dict_r[sorted_family_group_tuple]
	# 		# freq_list = counter_dict[sorted_family_group_tuple]

	# 	# print(str(sorted_family_group_tuple) + '\t' + str(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])) + '\t' + str(expected_freq_dict[(module_id, sorted_family_group_tuple)]) + '\t' + str(chi2) + '\t' + str(p_value))
	# 	# print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq:.2f}\t{observed_freq/expected_freq:.4f}\t{chi2:.2f}\t{manual_chi2:.2f}\t{p_value:.4f}\t{fisher_p_value:.4f}")
	# 	# print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq:.2f}\t{observed_freq/expected_freq:.4f}\t{manual_chi2:.2f}\t{fisher_p_value:.4f}\t{p_value_r:.4f}\t{expected_freq_r:.4f}\t{z_score:.4f}\t{freq_list}")
	# 	# print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq_r:.2f}\t{p_value_r:.4f}\t{z_score:.4f}\t{freq_list}")
	# 	print(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq_r:.2f}\t{p_value_r:.4f}\t{z_score:.4f}")
	# 	if write_to_file:
	# 		p_value_text = str(round(p_value_r, 3))
	# 		if p_value_r == 0.0:
	# 			p_value_text = '<' + str(1/randomized_round)
	# 		fp = open(output_fname, 'a')
	# 		fp.write(f"{sorted_family_group_tuple}\t{observed_freq}\t{expected_freq_r:.2f}\t{p_value_text}\t{z_score:.4f}\n")
	# 		fp.close()

# def calculate_expected_frequency(family_group_dict_with_ids, families, total_chains):
# 	expected_freq_dict = {}
# 	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
# 		prob_all = 1.0
# 		for family_id in sorted_family_group_tuple:
# 			prob_fam = len(families[family_id]) / total_chains
# 			prob_all *= prob_fam
# 		expected_freq_dict[(module_id, sorted_family_group_tuple)] = prob_all * total_chains

# 	return expected_freq_dict

# def chi_square_test(observed_freq, expected_freq, total_chains):
# 	# print(observed_freq, expected_freq, manual_chi_square([observed_freq], [expected_freq]))

# 	# chi2, p_value, dof, expected = chi2_contingency([observed_freq, expected_freq])
# 	# return chi2, p_value
# 	# return manual_chi_square([observed_freq], [expected_freq])

# 	observed = [observed_freq, total_chains - observed_freq]
# 	expected = [expected_freq, total_chains - expected_freq]
# 	chi2, p_value, _, _ = chi2_contingency([observed, expected])

def manual_chi_square(observed, expected):
	return sum((o - e) ** 2 / e for o, e in zip(observed, expected) if e != 0)

def generate_randomized_dataset(pdb_chainwise_residue_data, pdb_chainwise_loops, directories):
	randomized_pdb_chainwise_loops = {}
	all_coords = []
	# min_x = min_y = min_z = 99999
	# max_x = max_y = max_z = -99999
	pdb_count = len(pdb_chainwise_loops)
	for i, pdb_id in enumerate(pdb_chainwise_loops):
		print('Processing ' + str(pdb_id) + ' (' + str(i+1) + '/' + str(pdb_count) + ')')
		chain_cnt = len(pdb_chainwise_loops[pdb_id])

		randomized_pdb_chainwise_loops[pdb_id] = {}
		# pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
		# parser = FastMMCIFParser(QUIET=True)
		# pdb_structure = parser.get_structure('struct', pdb_fname)
		for j, chain_id in enumerate(pdb_chainwise_loops[pdb_id]):
			print('processing chain ' + str(chain_id) + ' (' + str(j+1) + '/' + str(chain_cnt) + ')')

			randomized_pdb_chainwise_loops[pdb_id][chain_id] = []
			# chain_structure = pdb_structure[0][chain_id]
			# residues = chain_structure.get_residues()
			# residues = pdb_chainwise_residue_data
			loop_info_dict, min_coords, max_coords = pdb_chainwise_residue_data[pdb_id][chain_id]
			# for res in residues:
			# 	for atom in res:
			# 		point = res[atom.name].get_coord()
			# 		all_coords.append(point)
			# 		# min_x, min_y, min_z = min(point[0], min_x), min(point[1], min_y), min(point[2], min_z)
			# 		# max_x, max_y, max_z = max(point[0], max_x), max(point[1], max_y), max(point[2], max_z)

			# # print(pdb_id, chain_id)
			# # print(all_coords)
			# # print(pdb_id, chain_id)
			# min_coords = np.min(all_coords, axis=0)
			# max_coords = np.max(all_coords, axis=0)
			# # print(min_coords)
			# # print(max_coords)
			# # print(min_x, min_y, min_z)
			# # print(max_x, max_y, max_z)
			# # sys.exit()

			# for loop in pdb_chainwise_loops[pdb_id][chain_id]:
			for loop in tqdm(pdb_chainwise_loops[pdb_id][chain_id], ncols=75):
				# motif_residues = get_residue_list_of_a_loop(loop, directories, chain_structure)
				## motif_residues = loop_info_dict[loop]
				# print(motif_residues)
				# sys.exit()
				## motif_coords = extract_motif_coordinates(motif_residues)
				residuewise_motif_coords = loop_info_dict[loop]
				motif_coords = []
				for coords in residuewise_motif_coords:
					motif_coords += coords
				# motif_coords = loop_info_dict[loop]
				# print('original')
				# print(len(motif_coords))
				# print(motif_coords)
				# sys.exit()
				all_randomized_coords = randomize_motif(motif_coords, min_coords, max_coords)

				residuewise_randomized_coords = []
				x = 0
				for coords in residuewise_motif_coords:
					randomized_coords = []
					for coord in coords:
						randomized_coords.append(all_randomized_coords[x])
						x += 1
					residuewise_randomized_coords.append(randomized_coords)

				# print('original')
				# print(motif_coords)
				# sys.exit()
				# print('randomized')
				# print(len(all_randomized_coords))
				# print(all_randomized_coords)
				# sys.exit()
				## randomized_motif_residues = update_motif_residues_with_randomized_coords(motif_residues, randomized_coords)
				## randomized_pdb_chainwise_loops[pdb_id][chain_id].append((loop, randomized_motif_residues))
				randomized_pdb_chainwise_loops[pdb_id][chain_id].append((loop, residuewise_randomized_coords))

	return randomized_pdb_chainwise_loops

def extract_motif_coordinates(motif_residues):
	"""Extracts the coordinates of a motif from the PDB structure."""
	residuewise_motif_coords = []
	for residue in motif_residues:
		coords = []
		for atom in residue:
			coords.append(residue[atom.name].get_coord().tolist())
		# print(coords)
		# sys.exit()
		residuewise_motif_coords.append(coords)
	# print(residuewise_motif_coords)
	# sys.exit()
	return residuewise_motif_coords

def randomize_motif(motif_coords, min_coords, max_coords):
	"""Applies random rotation and translation to motif coordinates."""
	motif_coords = np.array(motif_coords)
	
	# Generate random rotation matrix
	theta_x, theta_y, theta_z = np.random.uniform(0, 2 * np.pi, 3)
	# theta_x, theta_y, theta_z = (0, 0, 0)
	rotation_x = np.array([[1, 0, 0],
						   [0, np.cos(theta_x), -np.sin(theta_x)],
						   [0, np.sin(theta_x), np.cos(theta_x)]])
	rotation_y = np.array([[np.cos(theta_y), 0, np.sin(theta_y)],
						   [0, 1, 0],
						   [-np.sin(theta_y), 0, np.cos(theta_y)]])
	rotation_z = np.array([[np.cos(theta_z), -np.sin(theta_z), 0],
						   [np.sin(theta_z), np.cos(theta_z), 0],
						   [0, 0, 1]])
	rotation_matrix = rotation_x @ rotation_y @ rotation_z

	# Apply rotation
	rotated_coords = motif_coords @ rotation_matrix.T

	# Apply random translation
	# translation_vector = np.random.uniform(-10, 10, 3)  # Adjust range as needed
	translation_vector = np.random.uniform(min_coords, max_coords, 3)
	randomized_coords = rotated_coords + translation_vector

	return randomized_coords.tolist()

def validate_structure(motif_coords, all_coords, threshold=3.5):
	"""Validates the randomized structure for steric clashes."""
	clashes = 0
	for coord in motif_coords:
		dists = distance.cdist([coord], all_coords, 'euclidean')
		clashes += np.sum(dists < threshold)  # Threshold for clash detection
	return clashes == 0  # True if no clashes

# all_coords = [atom.get_coord() for atom in structure.get_atoms()]
# randomized_coords = randomize_motif(motif_coords)
# is_valid = validate_structure(randomized_coords, all_coords)
# print("Structure is valid:", is_valid)

def update_structure_with_motif(structure, motif_residues, randomized_coords):
	"""Updates the structure with new motif coordinates."""
	i = 0
	for residue in structure.get_residues():
		if residue.get_id()[1] in motif_residues:
			for atom in residue:
				atom.set_coord(randomized_coords[i])
				i += 1

def update_motif_residues_with_randomized_coords(motif_residues, randomized_coords):
	randomized_motif_residues = copy.deepcopy(motif_residues)
	i = 0
	for residue in randomized_motif_residues:
		for atom in residue:
			residue[atom.name].set_coord(randomized_coords[i])
			i += 1
	return randomized_motif_residues

# update_structure_with_motif(structure, motif_residues, randomized_coords)

# Save the updated structure
# io = PDBIO()
# io.set_structure(structure)
# io.save("randomized_rna.pdb")



def prepare_nearest_residue_data_stat(filtered_nearest_residue_data_dict):
	# in each module
	modulewise_proteins_dict = {}
	all_protein_list = ['ARG', 'LYS', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'ASN', 'GLN', 'HIS', 'GLY', 'GLU']
	for module_id, family_group_tuple in filtered_nearest_residue_data_dict:
		modulewise_proteins_dict[module_id, family_group_tuple] = []

		# in each loop group
		for filtered_nearest_residue_data in filtered_nearest_residue_data_dict[module_id, family_group_tuple]:
			cnt_dict = {}
			# chainwise items
			for chain_id in sorted(filtered_nearest_residue_data):
				# for r, res in filtered_nearest_residue_data[chain_id]:
				for r in filtered_nearest_residue_data[chain_id]:
					res_name = str(r).strip().split('(')[1].rstrip(')')
					# if len(res_name) > 1:
						# all_protein_list.append(res_name)
					if res_name not in cnt_dict:
						cnt_dict[res_name] = 0
					cnt_dict[res_name] += 1

			modulewise_proteins_dict[module_id, family_group_tuple].append(cnt_dict)

	# all_protein_list = list(set(all_protein_list))
	modulewise_proteins_stat_dict = {}
	for module_id, family_group_tuple in filtered_nearest_residue_data_dict:
		modulewise_proteins_stat_dict[module_id, family_group_tuple] = {}
		for prt in all_protein_list:
			numbers = []
			for cnt_dict in modulewise_proteins_dict[module_id, family_group_tuple]:
				if prt in cnt_dict:
					numbers.append(cnt_dict[prt])
				else:
					numbers.append(0)
			max_val, min_val, avg_val = max(numbers), min(numbers), sum(numbers)/len(numbers)
			modulewise_proteins_stat_dict[module_id, family_group_tuple][prt] = (max_val, min_val, avg_val)

	return modulewise_proteins_stat_dict


def find_loops_participating_in_modules(family_group_dict_with_ids):
	participating_loops = []
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		loop_list_of_lists = family_group_dict_with_ids[(module_id, sorted_family_group_tuple)]
		# print(loop_list_of_lists)
		# sys.exit()
		for id, loop_list in loop_list_of_lists:
			# print(loop_list)
			# sys.exit()
			for loop in loop_list:
				participating_loops.append(loop)
	# print(participating_loops)
	return list(set(participating_loops))

def find_loops_not_participating_in_modules(families, family_group_dict_with_ids, output_dir, write_to_file=True):
	output_fname = os.path.join(output_dir, 'nonmodule_loops_info.txt')

	familywise_loops_non_module = {}
	loops_participating_in_modules = find_loops_participating_in_modules(family_group_dict_with_ids)
	loops_non_module = []
	# print('Loops not participating in modules:')
	if write_to_file:
		fp = open(output_fname, 'w')
		fp.write('FamilyName\tNon-moduleLoopCount\tLoops\n')
		fp.close()

	for family_id in families:
		familywise_loops_non_module[family_id] = []
		loops = families[family_id]
		# print(family_id + ':')
		for loop in loops:
			if loop not in loops_participating_in_modules:
				familywise_loops_non_module[family_id].append(loop)
				loops_non_module.append(loop)
				# print(loop)

	# for family_id in familywise_loops_non_module:
	# 	print(family_id + ': ' + str(len(familywise_loops_non_module[family_id])))

	if write_to_file:
		fp = open(output_fname, 'a')
		for family_id in familywise_loops_non_module:
			fp.write(family_id + '\t' + str(len(familywise_loops_non_module[family_id])) + '\t' + ','.join(familywise_loops_non_module[family_id]) + '\n')
		fp.close()

	# sys.exit()

	return familywise_loops_non_module, loops_non_module



# Function to create a more detailed enrichment heatmap with customization options
def create_detailed_enrichment_heatmap(
	odds_ratio_data, 
	significance_matrix, 
	plot_output_fname,
	title,
	x_label,
	cmap=None,
	vmin=0.25,
	vmax=4.0,
	center=1.0,
	significance_symbol="*",
	show_diagonal=False,
	triangle="lower"
):
	"""
	Create a customized enrichment heatmap showing odds ratios with significance indicators.
	
	Parameters:
	- odds_ratio_data: DataFrame with odds ratios
	- significance_matrix: DataFrame indicating statistical significance (same shape as odds_ratio_data)
	- title: Title for the plot
	- cmap: Custom colormap (defaults to diverging red-white-blue)
	- vmin, vmax: Bounds for color scale
	- significance_symbol: Symbol to use for indicating significance
	- show_diagonal: Whether to show diagonal values
	- triangle: Which triangle to show - "lower", "upper", or "both"
	
	Returns:
	- Figure object
	"""

	odds_ratio_data = odds_ratio_data.astype(float)
	# nan_mask = np.isnan(odds_ratio_data)
	nan_mask = odds_ratio_data.isna()
	# Create masks based on requested triangle
	mask = np.zeros_like(odds_ratio_data, dtype=bool)
	
	if triangle == "lower":
		mask[np.triu_indices_from(mask, k=1 if show_diagonal else 0)] = True
	elif triangle == "upper":
		mask[np.tril_indices_from(mask, k=-1 if show_diagonal else 0)] = True
	elif triangle == "both":
		if not show_diagonal:
			mask[np.diag_indices_from(mask)] = True

	# Combine masks - this will mask both upper triangle and NaN values
	combined_mask = np.logical_or(mask, nan_mask)
	# odds_ratio_data = odds_ratio_data.fillna(0)
	plot_data = odds_ratio_data.copy()
	plot_data = plot_data.fillna(0)  # Temporary fill just for plotting

	# scaler = MinMaxScaler()
	# columns = plot_data.columns.tolist()
	# plot_data[columns] = scaler.fit_transform(plot_data[columns])
	
	# Default colormap if none provided
	if cmap is None:
		cmap = LinearSegmentedColormap.from_list(
			"enrichment_cmap", 
			[(0, "blue"), (0.5, "white"), (1, "red")],
			N=256
		)
	
	# Set up the figure
	fig, ax = plt.subplots(figsize=(10, 8))
	
	# # Create the heatmap
	# sns.heatmap(
	#	 odds_ratio_data,
	#	 mask=mask,
	#	 cmap=cmap,
	#	 vmin=vmin,
	#	 vmax=vmax,
	#	 center=1.0,
	#	 annot=True,
	#	 fmt=".1f",
	#	 linewidths=0.5,
	#	 ax=ax,
	#	 cbar_kws={"label": "Odds Ratio", "shrink": 0.8}
	# )
	# Create the heatmap
	sns.heatmap(
		plot_data,
		mask=combined_mask,
		cmap=cmap,
		vmin=vmin,
		vmax=vmax,
		center=center,
		annot=True,
		fmt=".1f",
		linewidths=0.5,
		ax=ax,
		cbar_kws={"label": x_label, "shrink": 0.8}
	)

	# For any visible cells that had NaN values originally, you can add special annotation
	# for i, j in zip(*np.where(np.logical_and(nan_mask, ~mask))):
	# 	ax.text(j + 0.5, i + 0.5, "N/A", ha='center', va='center', color='black')

	
	# Add hatching to cells with NaN but not in upper triangle
	for i, j in zip(*np.where(np.logical_and(nan_mask, ~mask))):
		ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='lightgray', hatch='///'))
		ax.text(j + 0.5, i + 0.5, "N/A", ha='center', va='center', color='black')

	# # Find cells with NaN and add hatching
	# for i, j in zip(*np.where(np.isnan(plot_data.values))):
	#	 ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='gray', hatch='///'))
	
	# Add significance indicators
	if significance_symbol:
		for i, row in enumerate(odds_ratio_data.index):
			for j, col in enumerate(odds_ratio_data.columns):
				# Only add symbols where we have values and significance
				if not mask[i, j] and not pd.isna(odds_ratio_data.iloc[i, j]) and significance_matrix.iloc[i, j]:
					# Find the text object for this cell and modify it
					for text_obj in ax.texts:
						if text_obj.get_position() == (j + 0.5, i + 0.5):
							current_text = text_obj.get_text()
							text_obj.set_text(f"{current_text}{significance_symbol}")
							break
	
	# Customize appearance
	ax.set_title(title, fontsize=16, pad=20)
	ax.set_xlabel("Motif Family", fontsize=12)
	ax.set_ylabel("Motif Family", fontsize=12)
	
	# Add legend for significance
	if significance_symbol:
		from matplotlib.lines import Line2D
		legend_elements = [
			Line2D([0], [0], marker='o', color='w', markerfacecolor='w', 
				  markersize=0, label=f"{significance_symbol} p < 0.05")
		]
		ax.legend(handles=legend_elements, loc='lower left')
	
	plt.tight_layout()
	# plt.show()

	plt.savefig(plot_output_fname)

	# return fig