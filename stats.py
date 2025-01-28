# import matplotlib.pyplot as plt
# import seaborn as sns
import os
import sys
import copy
import logging
import numpy as np
from Bio.PDB import *
from scipy.stats import chi2_contingency, fisher_exact
from graph_utils import *
from utils import *
from module import *
from logger import *
from tqdm import tqdm
import multiprocessing as mp

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
			print('round: ' + str(i))
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

def prepare_stats(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, total_chains, directories, output_dir, number_of_multiprocess, randomized_round, write_to_file=True):
	output_fname = os.path.join(output_dir, 'motif_module_stats.txt')

	logger.info('Generating statistical significance.')
	# print('stats start')
	p_value_dict, expected_freq_dict_r, counter_dict = calculate_p_value_and_expected_freq(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, directories, number_of_multiprocess, randomized_round)

	# expected_freq_dict = calculate_expected_frequency(family_group_dict_with_ids, families, total_chains)

	# print(chi2_contingency())

	# print('Module\tObservedFreq\tExpectedFreq\tRatio\tchi2\tmanual_chi2\tp_value\tfisher_p_value\tp_value_r\texpected_freq_r')
	# print('Module\tObservedFreq\tExpectedFreq\tRatio\tmanual_chi2\tfisher_p_value\tp_value_r\texpected_freq_r\tz_score')
	print('Module\tObservedFreq\tExpectedFreq\tp-value\tz-score')
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
			expected_freq_r = expected_freq_dict_r[sorted_family_group_tuple]
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
