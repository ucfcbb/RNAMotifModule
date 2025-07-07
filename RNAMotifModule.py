from collections import Counter
import matplotlib.pyplot as plt
import pickle
from tqdm import tqdm


import sys
import os
import time
import argparse
import logging
import collections
import multiprocessing as mp

# import networkx as nx
# from sklearn.cluster import AgglomerativeClustering

# from graphkernels.kernels import WeisfeilerLehman
# from sklearn.manifold import SpectralEmbedding
# from sklearn.cluster import OPTICS

# from grakel import Graph, Kernel
# from grakel.utils import graph_from_networkx
# from grakel.kernels import WeisfeilerLehman, VertexHistogram
# from sklearn.manifold import SpectralEmbedding
# from sklearn.cluster import OPTICS

# from karateclub import Graph2Vec
# from sklearn.cluster import KMeans


from classes import *
from config import *
from logger import *
from prepare import *
from converter import *
from partial_pdb_generator import *
# from profile import *
from utils import *
from graph_utils import *
# from input_utils import *
# from interaction_utils import *
# from scoring_utils import *
# from search import *

# from benchmark_helper import *
from stats import *
from module import *

from clustering import *

def main():

	parser = argparse.ArgumentParser(description='Identify structural motif modules from corresponding RNA structures.')

	parser.add_argument('-i', '--input', nargs='+', default=['input/known_families_IL.csv', 'input/known_families_HL.csv'], help='Space separated list of input files.')
	parser.add_argument('-d', '--distance', nargs='?', default='5.0', help='Distance threshold to be identified as spatially close.')
	parser.add_argument('-o', '--outputsubdir', nargs='?', default='', const='', help='Subdirectory inside the "output" directory to save the results.')
	parser.add_argument('-s', '--stat', nargs='?', default=False, const=True, help='Continue to generate statistical significance of the motif modules.')
	parser.add_argument('-c', '--cluster', nargs='?', default=False, const=True, help='Continue to cluster chainwise motifs to idenitfy motif modules from clustering point of view.')
	parser.add_argument('-x', '--randomizedrounds', nargs='?', default='100', const='100', help='Randomized dataset to generate to calculate expected frequency of the motif modules.')
	parser.add_argument('-r', '--residue', nargs='?', default=False, const=True, help='Continue to generate nearest residue data.')
	parser.add_argument('-t', '--residuestat', nargs='?', default=False, const=True, help='Continue to generate nearest residue stat.')
	parser.add_argument('-n', '--nonmodule', nargs='?', default=False, const=True, help='Include non-module information.')
	parser.add_argument('-p', '--partialpdb', nargs='?', default=False, const=True, help='Generate partial PDB files to visualize in PyMOL or similar tools.')
	parser.add_argument('-e', '--extension', nargs='?', default='0', const='0', help='Residues to extend beyond loop boundary to generate the partial PDB (loop.cif) files.')
	parser.add_argument('-k', '--pickle', nargs='?', default=False, const=True, help='Utilize Python pickle (recommended) to re-use previously generated data. Need to be deleted manually for a fresh run.')
	parser.add_argument('-a', '--atomset', nargs='?', default='0', const='0', help='Atom set to use in calculating distance between residues (0: all, 1: backbone and sugar, 2: backbone).')

	try:
		args = parser.parse_args()
	except Exception as e:
		parser.print_help()
		sys.exit()

	input_files = args.input
	distance_threshold_to_be_nearest_residue = float(args.distance)
	output_subdir_name = args.outputsubdir
	generate_stat = args.stat # True if args.stat.lower() == 'true' else False
	generate_clusters = args.cluster
	randomized_round = int(args.randomizedrounds)
	include_non_module = args.nonmodule # True if args.nonmodule.lower() == 'true' else False
	generate_nearest_res_data = args.residue # True if args.residue.lower() == 'true' else False
	generate_nearest_res_stat = args.residuestat # True if args.residuestat.lower() == 'true' else False
	generate_partial_pdb = args.partialpdb # True if args.partialpdb.lower() == 'true' else False
	loop_cif_extension = int(args.extension)
	utilize_pickle = args.pickle # True if args.pickle.lower() == 'true' else False
	atom_set_choice = int(args.atomset)

	if generate_nearest_res_stat or generate_partial_pdb:
		generate_nearest_res_data = True

	# print(input_files)
	# print(distance_threshold_to_be_nearest_residue)
	# print(generate_stat)
	# sys.exit()
	# print(include_non_module)
	# print(generate_partial_pdb)
	# print(utilize_pickle)

	# sys.exit()

	# distance_threshold_to_be_nearest_residue = 5

	directories = get_base_dir_names()
	directories.partial_pdb_dir = directories.partial_pdb_dir.replace('*', str(loop_cif_extension))
	create_required_directories(directories)

	input_index_type, annotation_source, items_per_chunk, mp_number_of_process, content_download_attempts = get_misc_params()

	# input_fname = 'input/known_families_IL.csv'
	# # input_fname = 'input/Subcluster_output_IL_PDB.csv.in'
	# # input_fname = 'input/cluster_output_PDB.csv'
	# # output_subdir_name = ''
	# # loop_cif_extension = 0
	# # directories.partial_pdb_dir = directories.partial_pdb_dir.replace('*', str(loop_cif_extension))

	# # families = {}
	# fp_input = open(input_fname)
	# loop_list = csv_to_list(fp_input.readlines())
	# fp_input.close()

	# input_fname = 'input/known_families_HL.csv'
	# # input_fname = 'input/Subcluster_output_HL_PDB.csv.in'
	# # input_fname = 'input/cluster_output_PDB.csv'

	# # families = {}
	# fp_input = open(input_fname)
	# # print(len(loop_list))
	# loop_list += csv_to_list(fp_input.readlines())
	# # print(len(loop_list))
	# fp_input.close()

	# fp_input = open('input/cluster_annotation_IL_HL.csv')
	# cl_ann_list = csv_to_list(fp_input.readlines())
	# fp_input.close()

	families = {}
	loop_list = []
	for input_fname in input_files:
		fp_input = open(input_fname)
		loop_list += csv_to_list(fp_input.readlines())
		fp_input.close()

	# print('\n'.join(list(map(lambda x: str(x), loop_list))))
	# sys.exit()

	# cl_annotation_dict = {}
	# for item in cl_ann_list:
	# 	if len(item) == 2:
	# 		loop_type = item[1].strip().split('_')[0]
	# 		cl_annotation_dict[loop_type + '_subcluster_' + str(item[0])] = item[1]

	
	for item in loop_list:
		if len(item) > 1:
			# families[item[0]] = map(lambda x: str(strToNode(x)), item[1:]) # item[1:]
			way_counts = get_way_counts(item[1:])
			loop_type = ''
			if len(way_counts) == 1:
				if way_counts[0] == 1:
					loop_type = 'HL'
				elif way_counts[0] == 2:
					loop_type = 'IL'

			# if loop_type + '_' + item[0] in cl_annotation_dict:
			# 	if cl_annotation_dict[loop_type + '_' + item[0]] in families:
			# 		families[cl_annotation_dict[loop_type + '_' + item[0]]] += item[1:]
			# 	else:
			# 		families[cl_annotation_dict[loop_type + '_' + item[0]]] = item[1:]
			# else:
			families[loop_type + '_' + item[0]] = item[1:]

	# print(families.keys())
	# sys.exit()
	prepare_data(families, directories, annotation_source, content_download_attempts, mp_number_of_process)
	if input_index_type == 'pdb':
		families = convert_a_cluster_from_PDB_to_FASTA(families, directories)
		for family_id in families:
			families[family_id] = list(map(lambda x: str(strToNode(x)), families[family_id]))

	loop_count = 0
	loop_node_list_str = []
	for family_id in families:
		loops = families[family_id]
		loop_count += len(loops)
		for loop in loops:
			loop = str(strToNode(loop))
			loop_node_list_str.append(loop)
	
	duplicates = [item for item, count in collections.Counter(loop_node_list_str).items() if count > 1]
	if len(duplicates) > 0:
		print('duplicates:')
		print(duplicates)
		# sys.exit()

	loop_node_list_str = sorted(list(set(loop_node_list_str)))
	# loop_node_list_str = sorted(loop_node_list_str)

	logger.info(str(loop_count) + ' loops (' + str(len(loop_node_list_str)) + ' unique) found in ' + str(len(families)) + ' famil' + ('ies' if len(families) > 1 else 'y') + '.')
	print('')
	# sys.exit()

	prepare_loop_files(loop_node_list_str, directories, annotation_source, mp_number_of_process, get_env())	#chkd
	prepare_partial_pdbs(directories, loop_node_list_str, loop_cif_extension, mp_number_of_process)

	output_dir = directories.output_dir
	if len(output_subdir_name) > 0:
		output_dir = os.path.join(output_dir, output_subdir_name)
	create_directory(output_dir)

	

	pdb_chainwise_loops = {}
	# c = 0
	for family_id in families:
		# if family_id not in ['E-loop', 'Kink-turn']:
		# 	continue
		for loop in families[family_id]:
			pdb_chain, regions = loop.strip().split(':')

			# if pdb_chain not in ['4V9F_0', '5J7L_DA', '7RQB_1A', '4V5O_BA', '7A0S_X', '3J7A_A', '4V88_A6', '4LFB_A']:
			# 	continue
			
			pdb_id, chain_id = pdb_chain.strip().split('_')

			if pdb_id not in pdb_chainwise_loops:
				pdb_chainwise_loops[pdb_id] = {}

			if chain_id not in pdb_chainwise_loops[pdb_id]:
				pdb_chainwise_loops[pdb_id][chain_id] = []

			pdb_chainwise_loops[pdb_id][chain_id].append(loop)
			# c += 1

	pdb_count = len(pdb_chainwise_loops)
	chain_count = sum([len(pdb_chainwise_loops[pdb_id]) for pdb_id in pdb_chainwise_loops])
	loop_count = sum([len(pdb_chainwise_loops[pdb_id][chain_id]) for pdb_id in pdb_chainwise_loops for chain_id in pdb_chainwise_loops[pdb_id]])

	print('Total PDBs: ' + str(pdb_count))
	print('Total Chains: ' + str(chain_count))
	print('Total loops: ' + str(loop_count))


	pdb_chainwise_residue_data = load_pdb_chainwise_residue_data_for_all_loops(pdb_chainwise_loops, directories, utilize_pickle)
	# print(pdb_chainwise_residue_data)
	# sys.exit()


	# randomized_pdb_chainwise_loops = generate_randomized_dataset(pdb_chainwise_loops, directories)
	# print(c)

	# for pdb_id in pdb_chainwise_loops:
	# 	for chain_id in pdb_chainwise_loops[pdb_id]:
	# 		print(pdb_id + '_' + chain_id + ': ' + str(len(pdb_chainwise_loops[pdb_id][chain_id])))

	# total_chains = 0
	# for pdb_id in pdb_chainwise_loops:
	# 	total_chains += len(pdb_chainwise_loops[pdb_id])

	# print(total_chains)
	# print(sum([len(pdb_chainwise_loops[pdb_id]) for pdb_id in pdb_chainwise_loops]))
	# sys.exit()

	spatial_proximity_data = load_spatial_proximity_data(pdb_chainwise_loops, directories, atom_set_choice, utilize_pickle)
	# print(spatial_proximity_data)
	# sys.exit()


	###### Experimental #######
	# graphs = []
	# for i, pdb_id in enumerate(pdb_chainwise_loops):
	# 	for j, chain_id in enumerate(pdb_chainwise_loops[pdb_id]):
	# 		loops = pdb_chainwise_loops[pdb_id][chain_id]
	# 		g = create_chain_graph(loops, families, pdb_id, chain_id, spatial_proximity_data, distance_threshold_to_be_nearest_residue)
	# 		graphs.append(g)

	# # print(len(graphs))
	# group_counts, total_groupwise_loops = analyze_motif_groups(graphs)
	# # print(group_counts)
	# frequent_groups = get_frequent_motif_groups(group_counts, min_support=0.0, min_size=2)
	# # print(frequent_groups)

	# print("Frequent Motif Groups:")
	# for group, count in frequent_groups:
	# 	print(f"{', '.join(sorted(group))}: {count} {total_groupwise_loops[group]}")

	###### Experimental #######


	# print_overall_nearest_motif_data(spatial_proximity_data, families)
	# print('\nFiltered or strict:')
	# print_overall_nearest_motif_data(filtered_spatial_proximity_data, families)

	# for pdb_id in pdb_chainwise_loops:
	# 	for chain_id in pdb_chainwise_loops[pdb_id]:
	# 		loops = pdb_chainwise_loops[pdb_id][chain_id]
	# 		graphs = []
	# 		g = create_chain_graph(loops, families, pdb_id, chain_id, filtered_spatial_proximity_data)
	# 		graphs.append(g)
	# 		group_counts = analyze_motif_groups(graphs)
	# 		frequent_groups = get_frequent_motif_groups(group_counts, min_support=0.1, min_size=2)

	# 		print("Frequent Motif Groups:")
	# 		for group, count in frequent_groups:
	# 			print(f"{', '.join(sorted(group))}: {count}")

	# visualize_motif_groups(frequent_groups)

	nearest_loop_data = {}
	for i, pdb_id in enumerate(pdb_chainwise_loops):
		nearest_loop_data[pdb_id] = {}
		for j, chain_id in enumerate(pdb_chainwise_loops[pdb_id]):
			nearest_loop_data[pdb_id][chain_id] = {}
			loops = pdb_chainwise_loops[pdb_id][chain_id]
			for loop in loops:
				if loop in nearest_loop_data[pdb_id][chain_id]:
					nearest_loop_data[pdb_id][chain_id][loop] += get_nearest_loop_list(pdb_id, chain_id, pdb_chainwise_loops, loop, spatial_proximity_data, distance_threshold_to_be_nearest_residue)
				else:
					nearest_loop_data[pdb_id][chain_id][loop] = get_nearest_loop_list(pdb_id, chain_id, pdb_chainwise_loops, loop, spatial_proximity_data, distance_threshold_to_be_nearest_residue)

	family_group_dict = get_motif_module_frequencies(nearest_loop_data, families)
	family_group_dict_with_ids = update_dict_with_ids(family_group_dict)

	print_overall_nearest_motif_data(family_group_dict_with_ids, directories, output_dir)

	# print('Testing clustering:')
	# clustering_analysis(families, pdb_chainwise_loops, pdb_chainwise_residue_data, spatial_proximity_data, distance_threshold_to_be_nearest_residue)

	# this part done, commenting temporarily
	# print("preparing stats")
	if generate_stat:
		# print('why here, then')
		# print(type(generate_stat))
		# print(generate_stat == False)
		# print(generate_stat == 'False')
		# sys.exit()
		prepare_stats(pdb_chainwise_residue_data, pdb_chainwise_loops, family_group_dict_with_ids, families, distance_threshold_to_be_nearest_residue, loop_count, directories, output_dir, mp_number_of_process, randomized_round)

	if generate_clusters:
		clustering_analysis(families, pdb_chainwise_loops, pdb_chainwise_residue_data, spatial_proximity_data, distance_threshold_to_be_nearest_residue, output_dir)
		# sys.exit()


	if generate_nearest_res_data:
		logger.info('Preparing nearest residue data.')
		# pickle_fname = 'filtered_nearest_residue_data-all_HL_IL.pkl3'
		pickle_fname = 'filtered_nearest_residue_data.pkl3'
		pickle_fname = os.path.join(directories.pickles_dir, pickle_fname)
		write_to_file = True
		output_fname = os.path.join(output_dir, 'nearest_protein_stat.txt')

		filtered_nearest_residue_data_dict = {}

		# if utilize_pickle and os.path.isfile(pickle_fname):
		# 	pf = open(pickle_fname, 'rb')
		# 	filtered_nearest_residue_data_dict = pickle.load(pf)
		# 	pf.close()
		# 	# print(spatial_proximity_data)
		# 	# sys.exit()
		# 	# generate_nearest_residue_data_write_only(filtered_nearest_residue_data_dict)
		# else:
		# 	filtered_nearest_residue_data_dict = generate_nearest_residue_data(family_group_dict_with_ids, directories, distance_threshold_to_be_nearest_residue, output_dir)
		# 	pf = open(pickle_fname, 'wb')
		# 	pickle.dump(filtered_nearest_residue_data_dict, pf)

		filtered_nearest_residue_data_dict = generate_nearest_residue_data(family_group_dict_with_ids, directories, distance_threshold_to_be_nearest_residue, output_dir)
		pf = open(pickle_fname, 'wb')
		pickle.dump(filtered_nearest_residue_data_dict, pf)

		if generate_nearest_res_stat:

			modulewise_proteins_stat_dict = prepare_nearest_residue_data_stat(filtered_nearest_residue_data_dict)

			if write_to_file:
				fp = open(output_fname, 'w')
				fp.write('Module\tAminoAcid1: (max, min, avg)\tAminoAcid2: (max, min, avg)\tAminoAcid3: (max, min, avg)\n')
				fp.close()

			for module_id, family_group_tuple in modulewise_proteins_stat_dict:
				group_names = ','.join(list(family_group_tuple))
				sorted_stat_dict = dict(sorted(modulewise_proteins_stat_dict[module_id, family_group_tuple].items(), key=lambda x: (x[1][0], x[1][2]), reverse=True))
				# print(str(family_group_tuple) + ': ' + str(modulewise_proteins_stat_dict[module_id, family_group_tuple]))
				sorted_stat_dict = {key : (sorted_stat_dict[key][0], sorted_stat_dict[key][1], round(sorted_stat_dict[key][2], 2)) for key in sorted_stat_dict}
				# print(group_names + ': ' + str(dict(list(sorted_stat_dict.items())[:3])))
				if write_to_file:
					item_line = ''
					# list(v)
					# list(map(lambda x: str(x), list(v)))
					flag = False
					temp_dict = dict(list(sorted_stat_dict.items())[:3])
					for k in temp_dict:
						if flag:
							item_line += '\t'
						item_line += str(k) + ': (' + ','.join(list(map(lambda x: str(x), list(temp_dict[k])))) + ')'
						flag = True

					fp = open(output_fname, 'a')
					fp.write(group_names + '\t' + item_line + '\n')
					fp.close()
			# sys.exit()

		if generate_partial_pdb:
			generate_partial_pdb_files_for_modules(directories, output_dir, family_group_dict_with_ids, loop_cif_extension, mp_number_of_process)

		if include_non_module:
			logger.info('Generating nearest residue data for non-module loops.')
			familywise_loops_non_module, all_non_module_loops = find_loops_not_participating_in_modules(families, family_group_dict_with_ids, output_dir)
			find_nearest_residue_data_for_non_module_loops(familywise_loops_non_module, all_non_module_loops, pdb_chainwise_loops, directories, distance_threshold_to_be_nearest_residue, output_dir)

def find_nearest_residue_data_for_non_module_loops(familywise_loops_non_module, all_non_module_loops, pdb_chainwise_loops, directories, distance_threshold_to_be_nearest_residue, output_dir, write_to_file=True):
	# logger.info('Finding nearest residue data for non-module loops.')
	# generate_filtered_nearest_protein_info = True
	# distance_threshold_to_be_nearest_residue = 5.0
	nearest_residue_data_dir = os.path.join(output_dir, 'nearest_residue_data_non_module_loops')
	create_directory(nearest_residue_data_dir)

	# nearest_residue_data_dict = {}
	# filtered_nearest_residue_data_dict = {}

	pdb_organism_details = load_pdb_organism_data(directories)
	# module_id = 1


	# pickle_fname = 'filtered_nearest_residue_data_for_non_module_loops-all_HL_IL.pkl3'
	pickle_fname = 'filtered_nearest_residue_data_for_non_module_loops.pkl3'
	pickle_fname = os.path.join(directories.pickles_dir, pickle_fname)

	nearest_residue_data_dict = {}
	filtered_nearest_residue_data_dict = {}
	if os.path.isfile(pickle_fname):
		pf = open(pickle_fname, 'rb')
		filtered_nearest_residue_data_dict = pickle.load(pf)
		pf.close()
		# print(spatial_proximity_data)
		# sys.exit()
	else:
		# filtered_nearest_residue_data_dict = generate_nearest_residue_data(family_group_dict_with_ids, directories)

		prev_pdb_id = None
		prev_chain_id = None
		pdb_structure = None
		for pdb_id in pdb_chainwise_loops:
			pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
			parser = FastMMCIFParser(QUIET=True)
			pdb_structure = parser.get_structure('struct', pdb_fname)
			for chain_id in pdb_chainwise_loops[pdb_id]:
				chain_structure = pdb_structure[0][chain_id]
				
				for loop in pdb_chainwise_loops[pdb_id][chain_id]:
					if loop in all_non_module_loops:
						# print([loop])
						# sys.exit()
						box_coord, loop_all_res_list = get_loops_bounding_box_coord([loop], directories, chain_structure)
						nearest_residue_data = get_nearest_residue_data_from_bounding_box(directories, pdb_id, chain_id, box_coord, pdb_structure, distance_threshold_to_be_nearest_residue)

						if len(nearest_residue_data) > 0:
							nearest_residue_data_dict[loop] = nearest_residue_data
							filtered_nearest_residue_data_dict[loop] = {}
							# if generate_filtered_nearest_protein_info == True:
							# 2nd step: get nearest info from loop atom coordinates
							filtered_nearest_residue_data = get_filtered_nearest_residue_data(loop_all_res_list, nearest_residue_data, distance_threshold_to_be_nearest_residue)
							filtered_nearest_residue_data_dict[loop] = filtered_nearest_residue_data

		pf = open(pickle_fname, 'wb')
		pickle.dump(filtered_nearest_residue_data_dict, pf)

	print('writing up found data.')
	for family_id in familywise_loops_non_module:
		filtered_nearest_residue_fname = os.path.join(nearest_residue_data_dir, 'nearest_residues_' + str(family_id) +'_filtered.txt') #'nearest_residue_data_non_module_loops/nearest_residues_' + str(family_id) +'_filtered.txt'
		# if generate_filtered_nearest_protein_info == True:
		fp_output = open(filtered_nearest_residue_fname, 'w')
		fp_output.write('Filtered nearest residue data for ' + str(family_id) + ' loops:\n\n')
		fp_output.close()


		loops = familywise_loops_non_module[family_id]
		for loop in loops:
			if loop in filtered_nearest_residue_data_dict:

				fp_output = open(filtered_nearest_residue_fname, 'a')
				write_nearest_residue_data_for_a_loop(directories, fp_output, loop, filtered_nearest_residue_data_dict[loop], pdb_organism_details, True)
				fp_output.close()

			else:
				fp_output = open(filtered_nearest_residue_fname, 'a')
				write_nearest_residue_data_for_a_loop(directories, fp_output, loop, {}, pdb_organism_details, True)
				# fp_output.write('Nothing found!\n\n')
				fp_output.close()


	print('preparing stat.')
	all_protein_list = ['ARG', 'LYS', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'ASN', 'GLN', 'HIS', 'GLY', 'GLU']
	familywise_proteins_dict = {}
	for family_id in familywise_loops_non_module:
		familywise_proteins_dict[family_id] = []
		for loop in familywise_loops_non_module[family_id]:
			cnt_dict = {}
			# chainwise items
			if loop in filtered_nearest_residue_data_dict:
				# print(filtered_nearest_residue_data_dict[loop])
				# sys.exit()
				for chain_id in sorted(filtered_nearest_residue_data_dict[loop]):
					# for r, res in filtered_nearest_residue_data[chain_id]:
					for r in filtered_nearest_residue_data_dict[loop][chain_id]:
						res_name = str(r).strip().split('(')[1].rstrip(')')
						# if len(res_name) > 1:
							# all_protein_list.append(res_name)
						if res_name not in cnt_dict:
							cnt_dict[res_name] = 0
						cnt_dict[res_name] += 1

			familywise_proteins_dict[family_id].append(cnt_dict)

	familywise_proteins_stat_dict = {}
	for family_id in familywise_proteins_dict:
		familywise_proteins_stat_dict[family_id] = {}
		for prt in all_protein_list:
			numbers = []
			for cnt_dict in familywise_proteins_dict[family_id]:
				if prt in cnt_dict:
					numbers.append(cnt_dict[prt])
				else:
					numbers.append(0)
			# print(numbers)
			# sys.exit()
			max_val = 0
			min_val = 999
			avg_val = 0.0
			if len(numbers) > 0:
				max_val, min_val, avg_val = max(numbers), min(numbers), sum(numbers)/len(numbers)
			familywise_proteins_stat_dict[family_id][prt] = (max_val, min_val, avg_val)


	write_to_file = True
	output_fname = os.path.join(output_dir, 'nearest_protein_stat_nonmodule.txt')

	if write_to_file:
		fp = open(output_fname, 'w')
		fp.write('Module\tAminoAcid1: (max, min, avg)\tAminoAcid2: (max, min, avg)\tAminoAcid3: (max, min, avg)\n')
		fp.close()

	for family_id in familywise_proteins_stat_dict:
		# group_names = ','.join(list(family_group_tuple))
		sorted_stat_dict = dict(sorted(familywise_proteins_stat_dict[family_id].items(), key=lambda x: (x[1][0], x[1][2]), reverse=True))
		# print(str(family_group_tuple) + ': ' + str(modulewise_proteins_stat_dict[module_id, family_group_tuple]))
		sorted_stat_dict = {key : (sorted_stat_dict[key][0], sorted_stat_dict[key][1], round(sorted_stat_dict[key][2], 2)) for key in sorted_stat_dict}
		# print(family_id + ': ' + str(dict(list(sorted_stat_dict.items())[:3])))

		if write_to_file:
			item_line = ''
			# list(v)
			# list(map(lambda x: str(x), list(v)))
			flag = False
			temp_dict = dict(list(sorted_stat_dict.items())[:3])
			for k in temp_dict:
				if flag:
					item_line += '\t'
				item_line += str(k) + ': (' + ','.join(list(map(lambda x: str(x), list(temp_dict[k])))) + ')'
				flag = True

			fp = open(output_fname, 'a')
			fp.write(family_id + '\t' + item_line + '\n')
			fp.close()


		# filtered_nearest_residue_data_dict[module_id, family_group_tuple] = []
		# commenting temporarily
		# nearest_residue_fname = 'nearest_residue_data/nearest_residues_module_' + str(module_id) +'.txt'
		# fp_output = open(nearest_residue_fname, 'w')
		# fp_output.write('Nearest residue data for "Module' + str(module_id) + '" loops ' + str(family_group_tuple) + ':\n\n')
		# fp_output.close()

		# commenting temporarily
		# filtered_nearest_residue_fname = 'nearest_residue_data/nearest_residues_module_' + str(module_id) +'_filtered.txt'
		# if generate_filtered_nearest_protein_info == True:
			# fp_output = open(filtered_nearest_residue_fname, 'w')
			# fp_output.write('Filtered nearest residue data for "Module' + str(module_id) + '" loops ' + str(family_group_tuple) + ':\n\n')
			# fp_output.close()

		# loop_list_of_lists = family_group_dict_with_ids[(module_id, family_group_tuple)]
	# 	loops = familywise_loops_non_module[family_id]

	# 	prev_pdb_id = None
	# 	prev_chain_id = None
	# 	pdb_structure = None
	# 	for loop_group_id, loop_list in loop_list_of_lists:
	# 		pdb_chain, regions = loop_list[0].strip().split(':')
	# 		pdb_id, chain_id = pdb_chain.strip().split('_')
			
	# 		if pdb_id != prev_pdb_id:
	# 			pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
	# 			parser = FastMMCIFParser(QUIET=True)
	# 			pdb_structure = parser.get_structure('struct', pdb_fname)

	# 		nearest_residue_data, filtered_nearest_residue_data = find_nearest_residues_from_a_loop_list(loop_list, pdb_structure, directories)
	# 		filtered_nearest_residue_data_dict[module_id, family_group_tuple].append(filtered_nearest_residue_data)

	# 		# commenting temporarily
	# 		# fp_output = open(nearest_residue_fname, 'a')
	# 		# write_nearest_residue_data(directories, pdb_structure, pdb_id, fp_output, module_id, loop_group_id, loop_list, nearest_residue_data, pdb_organism_details)
	# 		# fp_output.close()

	# 		# commenting temporarily
	# 		# fp_output = open(filtered_nearest_residue_fname, 'a')
	# 		# write_nearest_residue_data(directories, pdb_structure, pdb_id, fp_output, module_id, loop_group_id, loop_list, filtered_nearest_residue_data, pdb_organism_details, True)
	# 		# fp_output.close()

	# 	# module_id += 1
	# return filtered_nearest_residue_data_dict

# def update_dict_with_ids(family_group_dict):
# 	family_group_dict_with_ids = {}

# 	family_group_dict = dict(sorted(family_group_dict.items(), key=lambda item: len(item[1]), reverse=True))
# 	for module_id, family_group_tuple in enumerate(family_group_dict):
# 		family_group_dict_with_ids[(module_id+1, family_group_tuple)] = []
# 		for loop_group_id, loop_list in enumerate(family_group_dict[family_group_tuple]):
# 			family_group_dict_with_ids[(module_id+1, family_group_tuple)].append((loop_group_id+1, loop_list))

# 	return family_group_dict_with_ids

def generate_partial_pdb_files_for_modules(directories, output_dir, family_group_dict_with_ids, loop_cif_extension, number_of_multiprocess):
	logger.info('Generating partial PDB files for motif modules. (Loop extension being used: ' + str(loop_cif_extension) + ')')
	pdbx_dir = directories.pdbx_dir
	partial_pdbx_dir = os.path.join(output_dir, 'partial_pdb_for_modules')
	create_directory(partial_pdbx_dir)

	include_adjacent_info = True

	parameter_list = []

	for module_id, family_group_tuple in family_group_dict_with_ids:
		loop_list_of_lists = family_group_dict_with_ids[(module_id, family_group_tuple)]
		prev_pdb_id = None
		prev_chain_id = None
		pdb_structure = None
		for loop_group_id, loop_list in loop_list_of_lists:
			pdb_chain, regions = loop_list[0].strip().split(':')
			pdb_id, chain_id = pdb_chain.strip().split('_')
			
			if pdb_id != prev_pdb_id:
				pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
				parser = FastMMCIFParser(QUIET=True)
				pdb_structure = parser.get_structure('struct', pdb_fname)

			pdb_chain_dict, pdb_chain_adj_dict = generate_pdb_chain_dict_for_modules(directories, module_id, loop_group_id, loop_list, loop_cif_extension, include_adjacent_info)

			# print('pdb_chain_dict')
			# print(pdb_chain_dict)
			# print('pdb_chain_adj_dict')
			# print(pdb_chain_adj_dict)

			# if pdb_id in pdb_chain_adj_dict:
			# 	parameter_list.append((pdbx_dir, partial_pdbx_dir, pdb_id, module_id, loop_group_id, pdb_chain_dict[pdb_id], pdb_chain_adj_dict[pdb_id]))
			# else:
			# 	parameter_list.append((pdbx_dir, partial_pdbx_dir, pdb_id, module_id, loop_group_id, pdb_chain_dict[pdb_id], []))

			generate_partial_pdb_files_for_single_module(pdbx_dir, partial_pdbx_dir, pdb_id, module_id, loop_group_id, pdb_chain_dict[pdb_id], pdb_chain_adj_dict[pdb_id])

	# pool = mp.Pool(number_of_multiprocess)
	# pool.map(_generate_partial_pdb_files_worker, parameter_list)

# def _generate_partial_pdb_files_worker(p):
#	 generate_partial_pdb_files_for_single_pdb(*p)

# def generate_nearest_residue_data_write_only(filtered_nearest_residue_data_dict, output_dir):
# 	nearest_residue_data_dir = os.path.join(output_dir, 'nearest_residue_data')

# 	for module_id, family_group_tuple in filtered_nearest_residue_data_dict:
# 		filtered_nearest_residue_fname = os.path.join(nearest_residue_data_dir, 'nearest_residues_module_' + str(module_id) +'_filtered.txt')
# 		fp_output = open(filtered_nearest_residue_fname, 'w')
# 		fp_output.write('Filtered nearest residue data for "Module' + str(module_id) + '" loops ' + str(family_group_tuple) + ':\n\n')
# 		fp_output.close()

# 		loop_list_of_lists = family_group_dict_with_ids[(module_id, family_group_tuple)]

# 		for loop_group_id, loop_list in loop_list_of_lists:

# 		for filtered_nearest_residue_data in filtered_nearest_residue_data_dict[(module_id, family_group_tuple)]:
# 			fp_output = open(filtered_nearest_residue_fname, 'a')
# 			fp_output.write('Loop group: ' + str(loop_group_id) + '\n')
# 			for chain_id in filtered_nearest_residue_data:
# 				for r in filtered_nearest_residue_data[chain_id]:
# 			write_nearest_residue_data(directories, pdb_structure, pdb_id, fp_output, module_id, loop_group_id, loop_list, filtered_nearest_residue_data, pdb_organism_details, True)
# 			fp_output.close()



def generate_nearest_residue_data(family_group_dict_with_ids, directories, distance_threshold_to_be_nearest_residue, output_dir):
	# generate_filtered_nearest_protein_info = True
	family_group_dict_with_ids = dict(sorted(family_group_dict_with_ids.items(), key=lambda item: len(item[1]), reverse=True))
	# family_group_dict_with_adj_prtn = {}

	# nearest_residue_fname = '.'.join(input_fname.strip().split('.')[:-1]) + '_nearest_residues.txt'
	nearest_residue_data_dir = os.path.join(output_dir, 'nearest_residue_data')
	# create_directory('nearest_residue_data')
	create_directory(nearest_residue_data_dir)

	nearest_residue_data_dict = {}
	filtered_nearest_residue_data_dict = {}

	pdb_organism_details = load_pdb_organism_data(directories)
	# module_id = 1

	for module_id, family_group_tuple in family_group_dict_with_ids:
		filtered_nearest_residue_data_dict[module_id, family_group_tuple] = []
		# commenting temporarily
		nearest_residue_fname = os.path.join(nearest_residue_data_dir, 'nearest_residues_module_' + str(module_id) +'.txt')
		fp_output = open(nearest_residue_fname, 'w')
		fp_output.write('Nearest residue data for "Module' + str(module_id) + '" loops ' + str(family_group_tuple) + ':\n\n')
		fp_output.close()

		# commenting temporarily
		filtered_nearest_residue_fname = os.path.join(nearest_residue_data_dir, 'nearest_residues_module_' + str(module_id) +'_filtered.txt')
		# if generate_filtered_nearest_protein_info == True:
		fp_output = open(filtered_nearest_residue_fname, 'w')
		fp_output.write('Filtered nearest residue data for "Module' + str(module_id) + '" loops ' + str(family_group_tuple) + ':\n\n')
		fp_output.close()

		loop_list_of_lists = family_group_dict_with_ids[(module_id, family_group_tuple)]

		prev_pdb_id = None
		prev_chain_id = None
		pdb_structure = None
		for loop_group_id, loop_list in loop_list_of_lists:
			pdb_chain, regions = loop_list[0].strip().split(':')
			pdb_id, chain_id = pdb_chain.strip().split('_')
			
			if pdb_id != prev_pdb_id:
				pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
				parser = FastMMCIFParser(QUIET=True)
				pdb_structure = parser.get_structure('struct', pdb_fname)

			nearest_residue_data, filtered_nearest_residue_data = find_nearest_residues_from_a_loop_list(loop_list, pdb_structure, directories, distance_threshold_to_be_nearest_residue)
			filtered_nearest_residue_data_dict[module_id, family_group_tuple].append(filtered_nearest_residue_data)

			# commenting temporarily
			fp_output = open(nearest_residue_fname, 'a')
			is_filtered_residue_data = False
			write_nearest_residue_data(directories, pdb_structure, pdb_id, fp_output, module_id, loop_group_id, loop_list, nearest_residue_data, is_filtered_residue_data, pdb_organism_details)
			fp_output.close()

			# commenting temporarily
			fp_output = open(filtered_nearest_residue_fname, 'a')
			is_filtered_residue_data = True
			write_nearest_residue_data(directories, pdb_structure, pdb_id, fp_output, module_id, loop_group_id, loop_list, filtered_nearest_residue_data, is_filtered_residue_data, pdb_organism_details, True)
			fp_output.close()

		# module_id += 1
	return filtered_nearest_residue_data_dict

def find_nearest_residues_from_a_loop_list(loop_list, pdb_structure, directories, distance_threshold_to_be_nearest_residue):
	# generate_filtered_nearest_protein_info = True
	# distance_threshold_to_be_nearest_residue = 5.0

	pdb_chain, regions = loop_list[0].strip().split(':')
	pdb_id, chain_id = pdb_chain.strip().split('_')
	chain_structure = pdb_structure[0][chain_id]
	box_coord, loops_all_res_list = get_loops_bounding_box_coord(loop_list, directories, chain_structure)

	nearest_residue_data = get_nearest_residue_data_from_bounding_box(directories, pdb_id, chain_id, box_coord, pdb_structure, distance_threshold_to_be_nearest_residue)
	filtered_nearest_residue_data = {}

	if len(nearest_residue_data) > 0:
			# nearest_residue_data_dict[loop] = nearest_residue_data

			# if generate_filtered_nearest_protein_info == True:
			# 2nd step: get nearest info from loop atom coordinates
			filtered_nearest_residue_data = get_filtered_nearest_residue_data(loops_all_res_list, nearest_residue_data, distance_threshold_to_be_nearest_residue)
			
			# if len(filtered_nearest_residue_data) > 0:
				# filtered_nearest_residue_data_dict[loop] = filtered_nearest_residue_data
	return nearest_residue_data, filtered_nearest_residue_data

def get_nearest_residue_data_from_bounding_box(directories, pdb_id, loop_chain_id, box_coord, pdb_structure, distance_threshold, selected_chain_ids=None):
	nearest_residue_data = {}

	if pdb_structure == None:
		pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
		parser = FastMMCIFParser(QUIET=True)
		pdb_structure = parser.get_structure('struct', pdb_fname)

	for chain in pdb_structure[0]:
		chain_id = chain.get_id()
		if chain_id == loop_chain_id:
			continue
		if selected_chain_ids != None and chain_id not in selected_chain_ids:
			continue
		residues = chain.get_residues()
		for res in residues:
			for atom in res:
				# x, y, z = res[atom.name].get_vector()
				if is_atom_coord_in_box(res[atom.name].get_vector(), box_coord, distance_threshold):
					if chain_id not in nearest_residue_data:
						nearest_residue_data[chain_id] = []

					_, seqnum, icode = res.get_id()
					icode = icode.strip()
					r = Residue(Chainindex(chain_id, seqnum, icode), res.get_resname())
					nearest_residue_data[chain_id].append((r, res))
					# nearest_residue_data[chain_id].append(r)
					break
					# resname = res.get_resname()
	return nearest_residue_data

def get_filtered_nearest_residue_data(loop_all_res_list, nearest_residue_data, distance_threshold):
	filtered_nearest_residue_data = {}
	for chain_id in nearest_residue_data:
		for _, res in nearest_residue_data[chain_id]:
			for atom in res:
				if is_atom_close_to_any_loop_atom(res[atom.name].get_vector(), loop_all_res_list, distance_threshold) == True:
					if chain_id not in filtered_nearest_residue_data:
						filtered_nearest_residue_data[chain_id] = []

					_, seqnum, icode = res.get_id()
					icode = icode.strip()
					r = Residue(Chainindex(chain_id, seqnum, icode), res.get_resname())
					# filtered_nearest_residue_data[chain_id].append((r, res))
					filtered_nearest_residue_data[chain_id].append(r)
					break

	return filtered_nearest_residue_data

def get_loops_bounding_box_coord(loops, directories, chain_structure):
	# print(loops)
	# sys.exit()
	min_x = min_y = min_z = 99999
	max_x = max_y = max_z = -99999

	loops_all_res_list = []

	# loop_cif_extension = 5
	# partial_pdbx_dir = os.path.join(directories.data_dir, 'pdbx_extracted_ext' + str(loop_cif_extension))
	# print(directories.partial_pdb_dir)
	# directories.partial_pdb_dir = directories.partial_pdb_dir.replace('*', str(loop_cif_extension))
	# print(directories.partial_pdb_dir)
	# sys.exit()

	pdb_chain = loops[0].strip().split(':')[0]
	pdb_id, chain_id = pdb_chain.strip().split('_')

	index_list_for_loops = []
	for loop in loops:
		# pdb_fn = os.path.join(directories.partial_pdb_dir, loop.replace(':', '_') + '.cif')
		index_list = get_pdb_index_list(directories, loop)
		# print(index_list)
		index_list = list(map(lambda x: (x[0], int(x[1]), x[2]), index_list))

		index_list_for_loops += index_list

	# parser = FastMMCIFParser(QUIET=True)
	# structure = parser.get_structure('struct', pdb_fn)

	# chain = structure[0][chain_id]
	residues = chain_structure.get_residues()

	# print(index_list)
	# residue_coord_dict = {}
	for r in residues:
		# print(r.get_id())
		hetflag, resseq, icode = r.get_id()
		icode = icode.strip()
		ind = (chain_id, resseq, icode)
		# print('Searching ')
		# print(ind)
		if ind in index_list_for_loops:
			# if ind not in residue_coord_dict:
			#	 residue_coord_dict[ind] = []
			loops_all_res_list.append(r)

			for atom in r:
				point = r[atom.name].get_vector()

				min_x, min_y, min_z = min(point[0], min_x), min(point[1], min_y), min(point[2], min_z)
				max_x, max_y, max_z = max(point[0], max_x), max(point[1], max_y), max(point[2], max_z)
				# residue_coord_dict[ind].append(r[atom.name].get_vector())

	return (min_x, min_y, min_z, max_x, max_y, max_z), loops_all_res_list

# def get_nearest_loop_list(pdb_id, chain_id, pdb_chainwise_loops, loop1, spatial_proximity_data, distance_threshold):
# 	nearest_loop_list = []
# 	for loop2 in pdb_chainwise_loops[pdb_id][chain_id]:
# 		if loop1 == loop2:
# 			continue
# 		if spatial_proximity_data[pdb_id][chain_id][loop1][loop2] <= distance_threshold:
# 			nearest_loop_list.append(loop2)

# 	return nearest_loop_list

def create_chain_graph(motifs, families, pdb_id, chain_id, spatial_proximity_data, distance_threshold):
	G = nx.Graph()
	for i, motif in enumerate(motifs):
		G.add_node(i, family=get_family_id(motif, families), motif=str(motif))
	
	for i in range(len(motifs)):
		for j in range(len(motifs)):
			if i == j:
				continue
			if spatial_proximity_data[pdb_id][chain_id][motifs[i]][motifs[j]] <= distance_threshold:
				G.add_edge(i, j)

		# if motifs[i] not in filtered_spatial_proximity_data[pdb_id][chain_id]:
		# 	continue

		# filtered_nearest_loop_list = filtered_spatial_proximity_data[pdb_id][chain_id][motifs[i]]
		# loop_only_list = [loop for loop, res_list in filtered_nearest_loop_list]
		# for j in range(len(motifs)):
		# 	if motifs[j] in loop_only_list:
		# 		G.add_edge(i, j)
	
	return G

def find_motif_groups(graph):
	components = list(nx.connected_components(graph))
	# print(components)
	# sys.exit()

	frozenset_list = []
	frozensetwise_loops = {}
	for component in components:
		family_list = []
		loop_list = []
		for n in component:
			# print(graph.nodes[n])
			family_list.append(graph.nodes[n]['family'])
			loop_list.append(graph.nodes[n]['motif'])
		# print('')
		frozenset_list.append(frozenset(family_list))
		if frozenset(family_list) not in frozensetwise_loops:
			frozensetwise_loops[frozenset(family_list)] = []
		frozensetwise_loops[frozenset(family_list)].append(loop_list)

	# print(frozenset_list)
	# print([frozenset(graph.nodes[n]['family'] for n in component) for component in components])
	# sys.exit()

	# return [frozenset(graph.nodes[n]['family'] for n in component) for component in components]
	return frozenset_list, frozensetwise_loops

def analyze_motif_groups(chain_graphs):
	group_counts = Counter()
	total_groupwise_loops = {}
	
	for graph in chain_graphs:
		groups, groupwise_loops = find_motif_groups(graph)
		group_counts.update(groups)
		for group in groups:
			if group not in total_groupwise_loops:
				total_groupwise_loops[group] = []
			total_groupwise_loops[group] += groupwise_loops[group]

	
	return group_counts, total_groupwise_loops

def get_frequent_motif_groups(group_counts, min_support, min_size=2):
	total_chains = sum(group_counts.values())
	frequent_groups = [
		(group, count) for group, count in group_counts.items()
		if len(group) >= min_size and count / total_chains >= min_support
	]
	return sorted(frequent_groups, key=lambda x: x[1], reverse=True)

# def visualize_motif_groups(frequent_groups, top_n=20):
# 	groups, counts = zip(*frequent_groups[:top_n])
# 	group_labels = [', '.join(sorted(group)) for group in groups]
	
# 	plt.figure(figsize=(12, 8))
# 	plt.barh(range(len(counts)), counts, align='center')
# 	plt.yticks(range(len(counts)), group_labels)
# 	plt.xlabel('Frequency')
# 	plt.ylabel('Motif Groups')
# 	plt.title(f'Top {top_n} Frequent Motif Groups')
# 	plt.gca().invert_yaxis()  # To display the most frequent at the top
# 	plt.tight_layout()
# 	plt.show()

# def get_motif_module_frequencies(spatial_proximity_data, families):
# 	family_group_dict = {}
# 	for pdb_id in spatial_proximity_data:
# 		for chain_id in spatial_proximity_data[pdb_id]:
# 			for loop1 in spatial_proximity_data[pdb_id][chain_id]:
# 				if len(spatial_proximity_data[pdb_id][chain_id][loop1]) == 0:
# 					continue

# 				family_group = []
# 				loop_group = []
				
# 				family_group.append(get_family_id(loop1, families))
# 				loop_group.append(loop1)

# 				for loop2 in spatial_proximity_data[pdb_id][chain_id][loop1]:
# 					family_group.append(get_family_id(loop2, families))
# 					loop_group.append(loop2)

# 				sorted_family_group = sorted(family_group)
# 				sorted_loop_group = sorted(loop_group)

# 				sorted_family_group_tuple = tuple(sorted_family_group)

# 				if sorted_family_group_tuple not in family_group_dict:
# 					family_group_dict[sorted_family_group_tuple] = []
# 					family_group_dict[sorted_family_group_tuple].append(sorted_loop_group)
# 				else:
# 					if sorted_loop_group not in family_group_dict[sorted_family_group_tuple]:
# 						family_group_dict[sorted_family_group_tuple].append(sorted_loop_group)

# 	family_group_dict = dict(sorted(family_group_dict.items(), key=lambda item: len(item[1]), reverse=True))

# 	return family_group_dict

def print_overall_nearest_motif_data(family_group_dict_with_ids, directories, output_dir, write_to_file=True):
	output_fname = os.path.join(output_dir, 'identified_motif_modules.txt')
	
	if not write_to_file:
		print('Overall:')
	else:
		fp = open(output_fname, 'w')
		fp.write('Modules\tObservedFreq\tInstances\n')
		fp.close()

	# family_group_dict = {}
	# for pdb_id in spatial_proximity_data:
	# 	for chain_id in spatial_proximity_data[pdb_id]:
	# 		for loop1 in spatial_proximity_data[pdb_id][chain_id]:
	# 			if len(spatial_proximity_data[pdb_id][chain_id][loop1]) == 0:
	# 				continue

	# 			family_group = []
	# 			loop_group = []
				
	# 			family_group.append(get_family_id(loop1, families))
	# 			loop_group.append(loop1)

	# 			for loop2 in spatial_proximity_data[pdb_id][chain_id][loop1]:
	# 				family_group.append(get_family_id(loop2, families))
	# 				loop_group.append(loop2)

	# 			sorted_family_group = sorted(family_group)
	# 			sorted_loop_group = sorted(loop_group)

	# 			sorted_family_group_tuple = tuple(sorted_family_group)

	# 			if sorted_family_group_tuple not in family_group_dict:
	# 				family_group_dict[sorted_family_group_tuple] = []
	# 				family_group_dict[sorted_family_group_tuple].append(sorted_loop_group)
	# 			else:
	# 				if sorted_loop_group not in family_group_dict[sorted_family_group_tuple]:
	# 					family_group_dict[sorted_family_group_tuple].append(sorted_loop_group)

	# print(list(family_group_dict_with_ids.items())[0])
	# sys.exit()
	family_group_dict_with_ids = dict(sorted(family_group_dict_with_ids.items(), key=lambda item: len(item[1]), reverse=True))
	for module_id, sorted_family_group_tuple in family_group_dict_with_ids:
		# if len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)]) < 5:
		# 	print(str(sorted_family_group_tuple) + ' (' + str(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])) + ')')
		# else:
		# 	print(str(sorted_family_group_tuple) + ' (' + str(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])) + ')' + str(list(map(lambda x: [convert_a_loop_from_FASTA_to_PDB(item, directories) for item in x[1]], family_group_dict_with_ids[(module_id, sorted_family_group_tuple)]))))

		# printing FASTA ind
		if not write_to_file:
			print(str(sorted_family_group_tuple) + ' (' + str(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])) + ')' + str(list(map(lambda x: [item for item in x[1]], family_group_dict_with_ids[(module_id, sorted_family_group_tuple)]))))

		else:
			# '(' + ','.join(list(sorted_family_group_tuple)) + ')'
			# str(tuple([str(convert_a_loop_from_FASTA_to_PDB(item, directories)) for item in x[1]]))
			# '(' + ','.join([str(convert_a_loop_from_FASTA_to_PDB(item, directories)) for item in x[1]]) + ')'
			fp = open(os.path.join(output_dir, 'identified_motif_modules.txt'), 'a')
			fp.write('(' + ','.join(list(sorted_family_group_tuple)) + ')' + '\t' + str(len(family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])) + '\t' + str(','.join(list(map(lambda x: '(' + ','.join([str(convert_a_loop_from_FASTA_to_PDB(item, directories)) for item in x[1]]) + ')', family_group_dict_with_ids[(module_id, sorted_family_group_tuple)])))) + '\n')
			fp.close()

def print_nearest_loop_info(pdb_id, chain_id, spatial_proximity_data_in_a_chain, families):
	family_group_dict = {}

	print('PDB: ' + str(pdb_id) + ', chain: ' + str(chain_id))

	for loop1 in spatial_proximity_data_in_a_chain:
		if len(spatial_proximity_data_in_a_chain[loop1]) == 0:
			continue

		family_group = []
		loop_group = []
		
		family_group.append(get_family_id(loop1, families))
		loop_group.append(loop1)

		print(loop1 + '(' + get_family_id(loop1, families) + ') [' + str(len(spatial_proximity_data_in_a_chain[loop1])) + ']: ')
		for loop2, nearest_res_list in spatial_proximity_data_in_a_chain[loop1]:

			family_group.append(get_family_id(loop2, families))
			loop_group.append(loop2)

			print(loop2 + '(' + get_family_id(loop2, families) + ')', end='\t')
		print('')

		# print(family_group)
		# print(sorted(family_group))
		# sys.exit()
		sorted_family_group = sorted(family_group)
		sorted_loop_group = sorted(loop_group)

		sorted_family_group_tuple = tuple(sorted_family_group)

		if sorted_family_group_tuple not in family_group_dict:
			family_group_dict[sorted_family_group_tuple] = []
			family_group_dict[sorted_family_group_tuple].append(sorted_loop_group)
		else:
			if sorted_loop_group not in family_group_dict[sorted_family_group_tuple]:
				family_group_dict[sorted_family_group_tuple].append(sorted_loop_group)

	# print(family_group_dict)
	# sys.exit()
	print('\nSummary:')
	# print(family_group_dict.items())
	# sys.exit()
	family_group_dict = dict(sorted(family_group_dict.items(), key=lambda item: len(item[1]), reverse=True))
	for sorted_family_group_tuple in family_group_dict:
		print(str(sorted_family_group_tuple) + ' (' + str(len(family_group_dict[sorted_family_group_tuple])) + ')')
	print('\n')

# def get_family_id(loop, families):
# 	for family_id in families:
# 		if loop in families[family_id]:
# 			return family_id
# 	return 'N/A'

def get_loop_bounding_box_coord(loop, directories, chain_structure):
	min_x = min_y = min_z = 99999
	max_x = max_y = max_z = -99999

	loop_all_res_list = []

	# loop_cif_extension = 5
	# partial_pdbx_dir = os.path.join(directories.data_dir, 'pdbx_extracted_ext' + str(loop_cif_extension))
	# print(directories.partial_pdb_dir)
	# directories.partial_pdb_dir = directories.partial_pdb_dir.replace('*', str(loop_cif_extension))
	# print(directories.partial_pdb_dir)
	# sys.exit()

	pdb_chain = loop.strip().split(':')[0]
	pdb_id, chain_id = pdb_chain.strip().split('_')

	# pdb_fn = os.path.join(directories.partial_pdb_dir, loop.replace(':', '_') + '.cif')
	index_list = get_pdb_index_list(directories, loop)
	# print(index_list)
	index_list = list(map(lambda x: (x[0], int(x[1]), x[2]), index_list))

	# parser = FastMMCIFParser(QUIET=True)
	# structure = parser.get_structure('struct', pdb_fn)

	# chain = structure[0][chain_id]
	residues = chain_structure.get_residues()

	# print(index_list)
	# residue_coord_dict = {}
	for r in residues:
		# print(r.get_id())
		hetflag, resseq, icode = r.get_id()
		icode = icode.strip()
		ind = (chain_id, resseq, icode)
		# print('Searching ')
		# print(ind)
		if ind in index_list:
			# if ind not in residue_coord_dict:
			#	 residue_coord_dict[ind] = []
			loop_all_res_list.append(r)

			for atom in r:
				point = r[atom.name].get_vector()

				min_x, min_y, min_z = min(point[0], min_x), min(point[1], min_y), min(point[2], min_z)
				max_x, max_y, max_z = max(point[0], max_x), max(point[1], max_y), max(point[2], max_z)
				# residue_coord_dict[ind].append(r[atom.name].get_vector())

	return (min_x, min_y, min_z, max_x, max_y, max_z), loop_all_res_list

# def get_residue_list_of_a_loop(loop, directories, chain_structure):
# 	min_x = min_y = min_z = 99999
# 	max_x = max_y = max_z = -99999

# 	loop_all_res_list = []

# 	# loop_cif_extension = 5
# 	# partial_pdbx_dir = os.path.join(directories.data_dir, 'pdbx_extracted_ext' + str(loop_cif_extension))
# 	# print(directories.partial_pdb_dir)
# 	# directories.partial_pdb_dir = directories.partial_pdb_dir.replace('*', str(loop_cif_extension))
# 	# print(directories.partial_pdb_dir)
# 	# sys.exit()

# 	pdb_chain = loop.strip().split(':')[0]
# 	pdb_id, chain_id = pdb_chain.strip().split('_')

# 	# pdb_fn = os.path.join(directories.partial_pdb_dir, loop.replace(':', '_') + '.cif')
# 	index_list = get_pdb_index_list(directories, loop)
# 	# print(index_list)
# 	index_list = list(map(lambda x: (x[0], int(x[1]), x[2]), index_list))

# 	# parser = FastMMCIFParser(QUIET=True)
# 	# structure = parser.get_structure('struct', pdb_fn)

# 	# chain = structure[0][chain_id]
# 	residues = chain_structure.get_residues()

# 	# print(index_list)
# 	# residue_coord_dict = {}
# 	for r in residues:
# 		# print(r.get_id())
# 		hetflag, resseq, icode = r.get_id()
# 		icode = icode.strip()
# 		ind = (chain_id, resseq, icode)
# 		# print('Searching ')
# 		# print(ind)
# 		if ind in index_list:
# 			# if ind not in residue_coord_dict:
# 			#	 residue_coord_dict[ind] = []
# 			loop_all_res_list.append(r)

# 			# for atom in r:
# 			# 	point = r[atom.name].get_vector()

# 			# 	min_x, min_y, min_z = min(point[0], min_x), min(point[1], min_y), min(point[2], min_z)
# 			# 	max_x, max_y, max_z = max(point[0], max_x), max(point[1], max_y), max(point[2], max_z)
# 			# 	# residue_coord_dict[ind].append(r[atom.name].get_vector())

# 	return loop_all_res_list

# def get_spatial_proximity_data_of_a_loop(loop1, loop_info_dict, loops):
# 	spatial_proximity_data_of_a_loop = {}

# 	residue_list1 = loop_info_dict[loop1]

# 	for loop2 in loops:
# 		min_distance = 0
# 		if loop1 != loop2:
# 			residue_list2 = loop_info_dict[loop2]
# 			min_distance = find_minimum_distance_between_two_residue_lists(residue_list1, residue_list2)

# 		spatial_proximity_data_of_a_loop[loop2] = min_distance

# 	return spatial_proximity_data_of_a_loop

# def find_minimum_distance_between_two_residue_lists(residue_list1, residue_list2):
# 	min_distance = 999999

# 	for res1 in residue_list1:
# 		for res2 in residue_list2:
# 			distance = find_minimum_distance_between_two_residues(res1, res2)
# 			min_distance = min(min_distance, distance)

# 	return min_distance

# def find_minimum_distance_between_two_residues(res1, res2):
# 	min_distance = 999999

# 	for atom1 in res1:
# 		point1 = res1[atom1.name].get_vector()
# 		for atom2 in res2:
# 			point2 = res2[atom2.name].get_vector()
# 			distance = distance_between_points(point1, point2)
# 			min_distance = min(min_distance, distance)

# 	return min_distance

def get_filtered_nearest_loops_data(loop, loop_info_dict, loops, nearest_loop_list, distance_threshold):
# def get_filtered_nearest_residue_data(loop_all_res_list, nearest_residue_data, distance_threshold):
	filtered_nearest_loop_list = []
	box_coord1, residues1 = loop_info_dict[loop]

	# print(loop)
	for loop2 in loops:
		# print(loop2)
		if loop == loop2:
			continue

		box_coord2, residues2 = loop_info_dict[loop2]
		nearest_res_list = []
		for res in residues2:
			for atom in res:
				# x, y, z = res[atom.name].get_vector()
				if is_atom_close_to_any_loop_atom(res[atom.name].get_vector(), residues1, distance_threshold) == True:
				# if is_atom_coord_in_box(res[atom.name].get_vector(), box_coord1, distance_threshold):
					nearest_res_list.append(res)
					break

		if len(nearest_res_list) > 0:
			filtered_nearest_loop_list.append((loop2, nearest_res_list))

	return filtered_nearest_loop_list

def is_atom_close_to_any_loop_atom(point, loop_all_res_list, distance_threshold):
	# for res in loop_all_res_list:
	#	 min_x = min_y = min_z = 99999
	#	 max_x = max_y = max_z = -99999
	#	 for atom in res:
	#		 p = res[atom.name].get_vector()
	#		 min_x, min_y, min_z = min(p[0], min_x), min(p[1], min_y), min(p[2], min_z)
	#		 max_x, max_y, max_z = max(p[0], max_x), max(p[1], max_y), max(p[2], max_z)

	#	 if is_atom_coord_in_box(point, (min_x, min_y, min_z, max_x, max_y, max_z)) == True:
	#		 return True

	for res in loop_all_res_list:
		for atom in res:
			p = res[atom.name].get_vector()
			if distance_between_points(point, p) <= distance_threshold:
				return True

	return False

# def distance_between_points(p1, p2):
# 	x1, y1, z1 = p1
# 	x2, y2, z2 = p2
# 	return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def get_nearest_loops_from_bounding_box(loop, loop_info_dict, loops, distance_threshold):
	nearest_loop_list = []
	box_coord1, residues1 = loop_info_dict[loop]

	for loop2 in loops:
		if loop == loop2:
			continue

		box_coord2, residues2 = loop_info_dict[loop2]
		nearest_res_list = []
		for res in residues2:
			for atom in res:
				# x, y, z = res[atom.name].get_vector()
				if is_atom_coord_in_box(res[atom.name].get_vector(), box_coord1, distance_threshold):
					nearest_res_list.append(res)
					break

		if len(nearest_res_list) > 0:
			nearest_loop_list.append((loop2, nearest_res_list))

	return nearest_loop_list

def is_atom_coord_in_box(point, box_coord, distance_threshold):
	(min_x, min_y, min_z, max_x, max_y, max_z) = box_coord
	x, y, z = point
	if (min_x - distance_threshold <= x and x <= max_x + distance_threshold) and \
	(min_y - distance_threshold <= y and y <= max_y + distance_threshold) and \
	(min_z - distance_threshold <= z and z <= max_z + distance_threshold):
		return True

	return False

def write_nearest_residue_data(directories, pdb_structure, pdb_id, fp, module_id, loop_group_id, loops, nearest_data, is_filtered_residue_data, pdb_organism_details, write_adj_info_file=False):
	extended_distance_threshold_to_be_nearest_rna = 10
	extended_distance_threshold_to_be_nearest_protein = 20

	# lib_dir = 'my_lib'
	amino_acid_list = amino_acid_collection(os.path.join(directories.lib_dir, 'aminoacidlist.dat'))

	adjacent_info_dir = os.path.join(directories.data_dir, 'adjacent_info')
	create_directory(adjacent_info_dir)

	fp.write('Loop group: ' + str(loop_group_id) + '\n')
	pdb_loops = []
	for loop in loops:
		fp.write(loop + ' (FASTA)\n')
		pdb_loop = convert_a_loop_from_FASTA_to_PDB(loop, directories)
		pdb_loops.append(pdb_loop)
		fp.write(pdb_loop + ' (PDB)\n\n')
	
	pdb_chain, _ = loops[0].strip().split(':')
	pdb_id, chain_id = pdb_chain.strip().split('_')
	
	org_type, RNA_Type = get_organism_info(pdb_chain, pdb_organism_details)
	# fp.write(org_type.ljust(10) + ', ' + RNA_Type.ljust(30) + '\n')
	fp.write(chain_id.ljust(2) + ' (' + org_type.ljust(10) + ', ' + RNA_Type.ljust(30) + ')\n')

	adj_chain = None
	if len(nearest_data) > 0:
		adj_chain = list(nearest_data.keys())[0]
	adj_res_count = 0
	for chain_id in sorted(nearest_data):
		pdb_chain = str(pdb_id) + '_' + str(chain_id)
		
		org_type, RNA_Type = get_organism_info(pdb_chain, pdb_organism_details)

		# print(nearest_data[chain_id])

		if is_filtered_residue_data:
			# str_data = list(map(lambda x: str(x[0]), sorted(nearest_data[chain_id], key=lambda x: x[0].index.seqnum)))
			str_data = list(map(lambda x: str(x), sorted(nearest_data[chain_id], key=lambda x: x.index.seqnum)))
			# protein_cnt = reduce(lambda count, item: count + (item[0].symbol in amino_acid_list), nearest_data[chain_id], 0)
			protein_cnt = reduce(lambda count, item: count + (item.symbol in amino_acid_list), nearest_data[chain_id], 0)
		else:
			str_data = list(map(lambda x: str(x[0]), sorted(nearest_data[chain_id], key=lambda x: x[0].index.seqnum)))
			# str_data = list(map(lambda x: str(x), sorted(nearest_data[chain_id], key=lambda x: x.index.seqnum)))
			protein_cnt = reduce(lambda count, item: count + (item[0].symbol in amino_acid_list), nearest_data[chain_id], 0)
			# protein_cnt = reduce(lambda count, item: count + (item.symbol in amino_acid_list), nearest_data[chain_id], 0)

		if protein_cnt > adj_res_count:
			adj_chain = chain_id
			adj_res_count = protein_cnt

		fp.write(chain_id.ljust(2) + ' (' + org_type.ljust(10) + ', ' + RNA_Type.ljust(30) + ') [Total: ' + str(len(nearest_data[chain_id])) + ', Protein: ' + str(protein_cnt) + '] : ' + ', '.join(str_data) + '\n')

	# writing adjacent residue info for each loop in separate file
	if adj_chain == None and write_adj_info_file == True:
		adjacent_info_fname = os.path.join(adjacent_info_dir, 'module' + str(module_id) + '_loopgroup' + str(loop_group_id) + '.adj_info')
		fp_adj_info = open(adjacent_info_fname, 'w')
		# fp_adj_info.write(pdb_loop + ',' + str(loop_extension) + '\n')
		fp_adj_info.write(str(pdb_loops) + '\n')
		fp_adj_info.write('\n\n')
		fp_adj_info.write('1,0.5')
		fp_adj_info.close()

	if adj_chain != None and write_adj_info_file == True:
		adjacent_info_fname = os.path.join(adjacent_info_dir, 'module' + str(module_id) + '_loopgroup' + str(loop_group_id) + '.adj_info')
		fp_adj_info = open(adjacent_info_fname, 'w')
		# fp_adj_info.write(pdb_loop + ',' + str(loop_extension) + '\n')
		fp_adj_info.write(str(pdb_loops) + '\n')

		if is_filtered_residue_data:
			start_ind = min(list(map(lambda x: x.index.seqnum, nearest_data[adj_chain])))
			max_ind = max(list(map(lambda x: x.index.seqnum, nearest_data[adj_chain])))
		else:
			start_ind = min(list(map(lambda x: x[0].index.seqnum, nearest_data[adj_chain])))
			max_ind = max(list(map(lambda x: x[0].index.seqnum, nearest_data[adj_chain])))

		end_ind = min(start_ind + 5, max_ind)
		fp_adj_info.write(adj_chain + ':' + str(start_ind) + '-' + str(end_ind) + '\n')

		nearest_rna_data = {}
		nearest_protein_data = {}
		for chain_id in nearest_data:
			if len(nearest_data[chain_id]) > 0:

				if is_filtered_residue_data:
					if nearest_data[chain_id][0].symbol in amino_acid_list:
						nearest_protein_data[chain_id] = nearest_data[chain_id]
					else:
						nearest_rna_data[chain_id] = nearest_data[chain_id]
				else:
					if nearest_data[chain_id][0][0].symbol in amino_acid_list:
						nearest_protein_data[chain_id] = nearest_data[chain_id]
					else:
						nearest_rna_data[chain_id] = nearest_data[chain_id]

		# print('nearest data')
		# lengths = []
		# for chain_id in nearest_data:
		#	 lengths.append(str(len(nearest_data[chain_id])))
		# print(','.join(sorted(lengths)))

		# print('nearest rna data')
		# lengths = []
		# for chain_id in nearest_rna_data:
		#	 lengths.append(str(len(nearest_rna_data[chain_id])))
		# print(','.join(sorted(lengths)))

		# print('nearest pretein data')
		# lengths = []
		# for chain_id in nearest_protein_data:
		#	 lengths.append(str(len(nearest_protein_data[chain_id])))
		# print(','.join(sorted(lengths)))

		loop_chain_id = loops[0].strip().split(':')[0].strip().split('_')[1]
		chain_structure = pdb_structure[0][loop_chain_id]
		box_coord, loops_all_res_list = get_loops_bounding_box_coord(loops, directories, chain_structure)
		 

		load_regions = []

		if len(nearest_protein_data) > 0:
			nearest_protein_data = get_nearest_residue_data_from_bounding_box(directories, pdb_id, loop_chain_id, box_coord, None, extended_distance_threshold_to_be_nearest_protein, nearest_protein_data.keys())
			extended_nearest_protein_data = get_filtered_nearest_residue_data(loops_all_res_list, nearest_protein_data, extended_distance_threshold_to_be_nearest_protein)

			# print('extended nearest protein data')
			# lengths = []
			# for chain_id in extended_nearest_protein_data:
			#	 lengths.append(str(len(extended_nearest_protein_data[chain_id])))
			# print(','.join(sorted(lengths)))

			for chain_id in extended_nearest_protein_data:
				if len(extended_nearest_protein_data[chain_id]) > 0:
					if is_filtered_residue_data:
						sorted_data = sorted(extended_nearest_protein_data[chain_id], key=lambda x: x.index.seqnum)
					else:
						sorted_data = sorted(extended_nearest_protein_data[chain_id], key=lambda x: x[0].index.seqnum)
					# load_regions.append(chain_id + ':' + str(sorted_data[0][0].index.seqnum) + '-' + str(sorted_data[-1][0].index.seqnum))
					region_list = get_regions_list(sorted_data, is_filtered_residue_data)
					load_regions.append(chain_id + ':' + '_'.join(region_list))

		# print_a_dict_sorted(nearest_data)
		# print_a_dict_sorted(nearest_rna_data)
		# print_a_dict_sorted(nearest_protein_data)
		# print_a_dict_sorted(extended_nearest_rna_data)
		# print_a_dict_sorted(extended_nearest_protein_data)

		# fp_adj_info.write(','.join(nearest_data.keys()) + '\n')
		fp_adj_info.write(','.join(load_regions) + '\n')
		fp_adj_info.write('1,0.5')
		fp_adj_info.close()

	fp.write('\n\n')

def write_nearest_residue_data_for_a_loop(directories, fp, loop, nearest_data, pdb_organism_details, write_adj_info_file=False):
# def write_nearest_residue_data_for_a_loop(directories, pdb_structure, pdb_id, fp, module_id, loop_group_id, loops, nearest_data, pdb_organism_details, write_adj_info_file=False):
	extended_distance_threshold_to_be_nearest_rna = 10
	extended_distance_threshold_to_be_nearest_protein = 20

	# lib_dir = 'my_lib'
	amino_acid_list = amino_acid_collection(os.path.join(directories.lib_dir, 'aminoacidlist.dat'))

	# adjacent_info_dir = 'adjacent_info'
	# create_directory(adjacent_info_dir)

	# fp.write('Loop group: ' + str(loop_group_id) + '\n')
	# pdb_loops = []
	# for loop in loops:
	fp.write(loop + ' (FASTA)\n')
	pdb_loop = convert_a_loop_from_FASTA_to_PDB(loop, directories)
	# pdb_loops.append(pdb_loop)
	fp.write(pdb_loop + ' (PDB)\n\n')
	
	pdb_chain, _ = loop.strip().split(':')
	pdb_id, chain_id = pdb_chain.strip().split('_')
	
	org_type, RNA_Type = get_organism_info(pdb_chain, pdb_organism_details)
	# fp.write(org_type.ljust(10) + ', ' + RNA_Type.ljust(30) + '\n')
	fp.write(chain_id.ljust(2) + ' (' + org_type.ljust(10) + ', ' + RNA_Type.ljust(30) + ')\n')

	adj_chain = None
	if len(nearest_data) > 0:
		adj_chain = list(nearest_data.keys())[0]
	adj_res_count = 0
	for chain_id in sorted(nearest_data):
		pdb_chain = str(pdb_id) + '_' + str(chain_id)
		
		org_type, RNA_Type = get_organism_info(pdb_chain, pdb_organism_details)

		# str_data = list(map(lambda x: str(x[0]), sorted(nearest_data[chain_id], key=lambda x: x[0].index.seqnum)))
		str_data = list(map(lambda x: str(x), sorted(nearest_data[chain_id], key=lambda x: x.index.seqnum)))
		# protein_cnt = reduce(lambda count, item: count + (item[0].symbol in amino_acid_list), nearest_data[chain_id], 0)
		protein_cnt = reduce(lambda count, item: count + (item.symbol in amino_acid_list), nearest_data[chain_id], 0)

		if protein_cnt > adj_res_count:
			adj_chain = chain_id
			adj_res_count = protein_cnt

		fp.write(chain_id.ljust(2) + ' (' + org_type.ljust(10) + ', ' + RNA_Type.ljust(30) + ') [Total: ' + str(len(nearest_data[chain_id])) + ', Protein: ' + str(protein_cnt) + '] : ' + ', '.join(str_data) + '\n')

	# # writing adjacent residue info for each loop in separate file
	# if adj_chain == None and write_adj_info_file == True:
	# 	 adjacent_info_fname = os.path.join(adjacent_info_dir, 'module' + str(module_id) + '_loopgroup' + str(loop_group_id) + '.adj_info')
	# 	 fp_adj_info = open(adjacent_info_fname, 'w')
	# 	 # fp_adj_info.write(pdb_loop + ',' + str(loop_extension) + '\n')
	# 	 fp_adj_info.write(str(pdb_loops) + '\n')
	# 	 fp_adj_info.write('\n\n')
	# 	 fp_adj_info.write('1,0.5')
	# 	 fp_adj_info.close()

	# if adj_chain != None and write_adj_info_file == True:
	# 	 adjacent_info_fname = os.path.join(adjacent_info_dir, 'module' + str(module_id) + '_loopgroup' + str(loop_group_id) + '.adj_info')
	# 	 fp_adj_info = open(adjacent_info_fname, 'w')
	# 	 # fp_adj_info.write(pdb_loop + ',' + str(loop_extension) + '\n')
	# 	 fp_adj_info.write(str(pdb_loops) + '\n')

	# 	 start_ind = min(list(map(lambda x: x[0].index.seqnum, nearest_data[adj_chain])))
	# 	 max_ind = max(list(map(lambda x: x[0].index.seqnum, nearest_data[adj_chain])))
	# 	 end_ind = min(start_ind + 5, max_ind)
	# 	 fp_adj_info.write(adj_chain + ':' + str(start_ind) + '-' + str(end_ind) + '\n')

	# 	 nearest_rna_data = {}
	# 	 nearest_protein_data = {}
	# 	 for chain_id in nearest_data:
	# 		 if len(nearest_data[chain_id]) > 0:
	# 			 if nearest_data[chain_id][0][0].symbol in amino_acid_list:
	# 				 nearest_protein_data[chain_id] = nearest_data[chain_id]
	# 			 else:
	# 				 nearest_rna_data[chain_id] = nearest_data[chain_id]

	# 	 # print('nearest data')
	# 	 # lengths = []
	# 	 # for chain_id in nearest_data:
	# 	 #	 lengths.append(str(len(nearest_data[chain_id])))
	# 	 # print(','.join(sorted(lengths)))

	# 	 # print('nearest rna data')
	# 	 # lengths = []
	# 	 # for chain_id in nearest_rna_data:
	# 	 #	 lengths.append(str(len(nearest_rna_data[chain_id])))
	# 	 # print(','.join(sorted(lengths)))

	# 	 # print('nearest pretein data')
	# 	 # lengths = []
	# 	 # for chain_id in nearest_protein_data:
	# 	 #	 lengths.append(str(len(nearest_protein_data[chain_id])))
	# 	 # print(','.join(sorted(lengths)))

	# 	 loop_chain_id = loops[0].strip().split(':')[0].strip().split('_')[1]
	# 	 chain_structure = pdb_structure[0][loop_chain_id]
	# 	 box_coord, loops_all_res_list = get_loops_bounding_box_coord(loops, directories, chain_structure)
		 

	# 	 load_regions = []
	# 	 ####### Temporarily skipping RNA chains #######
	# 	 # if len(nearest_rna_data) > 0:
	# 	 #	 nearest_rna_data, pdb_structure = get_nearest_residue_data_from_bounding_box(pdb_id, loop_chain_id, box_coord, None, extended_distance_threshold_to_be_nearest_rna, nearest_rna_data.keys())
	# 	 #	 extended_nearest_rna_data = get_filtered_nearest_residue_data(loop_all_res_list, nearest_rna_data, extended_distance_threshold_to_be_nearest_rna)

	# 	 #	 for chain_id in extended_nearest_rna_data:
	# 	 #		 if len(extended_nearest_rna_data[chain_id]) > 0:
	# 	 #			 sorted_data = sorted(extended_nearest_rna_data[chain_id], key=lambda x: x[0].index.seqnum)
	# 	 #			 region_list = get_regions_list(sorted_data)
	# 	 #			 load_regions.append(chain_id + ':' + '_'.join(region_list))

	# 	 if len(nearest_protein_data) > 0:
	# 		 nearest_protein_data = get_nearest_residue_data_from_bounding_box(directories, pdb_id, loop_chain_id, box_coord, None, extended_distance_threshold_to_be_nearest_protein, nearest_protein_data.keys())
	# 		 extended_nearest_protein_data = get_filtered_nearest_residue_data(loops_all_res_list, nearest_protein_data, extended_distance_threshold_to_be_nearest_protein)

	# 		 # print('extended nearest protein data')
	# 		 # lengths = []
	# 		 # for chain_id in extended_nearest_protein_data:
	# 		 #	 lengths.append(str(len(extended_nearest_protein_data[chain_id])))
	# 		 # print(','.join(sorted(lengths)))

	# 		 for chain_id in extended_nearest_protein_data:
	# 			 if len(extended_nearest_protein_data[chain_id]) > 0:
	# 				 sorted_data = sorted(extended_nearest_protein_data[chain_id], key=lambda x: x[0].index.seqnum)
	# 				 # load_regions.append(chain_id + ':' + str(sorted_data[0][0].index.seqnum) + '-' + str(sorted_data[-1][0].index.seqnum))
	# 				 region_list = get_regions_list(sorted_data)
	# 				 load_regions.append(chain_id + ':' + '_'.join(region_list))

	# 	 # print_a_dict_sorted(nearest_data)
	# 	 # print_a_dict_sorted(nearest_rna_data)
	# 	 # print_a_dict_sorted(nearest_protein_data)
	# 	 # print_a_dict_sorted(extended_nearest_rna_data)
	# 	 # print_a_dict_sorted(extended_nearest_protein_data)

	# 	 # fp_adj_info.write(','.join(nearest_data.keys()) + '\n')
	# 	 fp_adj_info.write(','.join(load_regions) + '\n')
	# 	 fp_adj_info.write('1,0.5')
	# 	 fp_adj_info.close()

	fp.write('\n\n')

def load_pdb_organism_data(directories):
	# lib_dir = 'my_lib'
	pdb_organism_details = read_pdb_chain_organism_details(os.path.join(directories.lib_dir, 'PDB_Chain_Organism_Details.tsv'))
	pdb_organism_details_scrapped = read_pdb_chain_organism_details(os.path.join(directories.lib_dir, 'PDB_Chain_Organism_Details_scrapped.tsv'))

	for pdb_chain in pdb_organism_details_scrapped:
		if pdb_chain not in pdb_organism_details:
			pdb_organism_details[pdb_chain] = pdb_organism_details_scrapped[pdb_chain]

	return pdb_organism_details

def get_organism_info(pdb_chain, pdb_organism_details):
	org_type = ''
	RNA_Type = ''
	if pdb_chain in pdb_organism_details:
		RNA_Type, organism, org_class, org_type, pdb_source = pdb_organism_details[pdb_chain]

	RNA_Type_lower = RNA_Type.lower()
	if 'ribosomal subunit protein' in RNA_Type_lower:
		ind = RNA_Type_lower.index('ribosomal subunit protein')
		cnt = len('ribosomal subunit protein')
		RNA_Type = RNA_Type[:ind] + 'Ribo-Prot' + RNA_Type[ind+cnt:]

	if 'ribosomal protein' in RNA_Type_lower:
		ind = RNA_Type_lower.index('ribosomal protein')
		cnt = len('ribosomal protein')
		RNA_Type = RNA_Type[:ind] + 'Ribo-Prot' + RNA_Type[ind+cnt:]

	org_type = org_type if len(org_type) <= 10 else org_type[:11]
	RNA_Type = RNA_Type if len(RNA_Type) <= 30 else RNA_Type[:31]

	return org_type, RNA_Type

def read_pdb_chain_organism_details(fname):
	fp = open(fname)
	first_line = True
	pdb_organism_details = {}
	for line in fp.readlines():
		if first_line:
			first_line = False
			continue
		pieces = line.strip('\n').strip('\r').split('\t')
		# For each chain, store RNA Types, Organism, Class, Type (Manually Defined), Source
		if len(pieces) > 0:
			pdb_organism_details[pieces[0].strip()] = pieces[1:]
	fp.close()
	return pdb_organism_details

def get_regions_list(sorted_data, is_filtered_residue_data):
	# print(sorted_data)
	region_list = []
	if is_filtered_residue_data:
		s = sorted_data[0].index.seqnum
	else:
		s = sorted_data[0][0].index.seqnum
	prev_i_seqnum = s
	for i in range(1, len(sorted_data)):
		if is_filtered_residue_data:
			i_seqnum = sorted_data[i].index.seqnum
		else:
			i_seqnum = sorted_data[i][0].index.seqnum

		if i_seqnum - prev_i_seqnum > 5:
			region_list.append((s, prev_i_seqnum))
			s = i_seqnum
		prev_i_seqnum = i_seqnum
	region_list.append((s, prev_i_seqnum))
	region_list = list(map(lambda x: str(x[0]) + '-' + str(x[1]), region_list))

	return region_list

def load_chain_data(chain_structure):
	residues = list(chain_structure.get_residues())
	all_coords = []
	for res in residues:
		for atom in res:
			point = res[atom.name].get_coord()
			all_coords.append(point)
	min_coords = np.min(all_coords, axis=0)
	max_coords = np.max(all_coords, axis=0)

	return residues, min_coords, max_coords

def get_loop_specific_residues(loop, residues, chain_id, directories):
	loop_specific_residue_list = []

	index_list = get_pdb_index_list(directories, loop)
	index_list = list(map(lambda x: (x[0], int(x[1]), x[2]), index_list))

	# print(loop)
	# print('index_list')
	# print(index_list)
	# print(residues)
	# sys.exit()

	for residue in residues:
		# print(r.get_id())
		hetflag, resseq, icode = residue.get_id()
		icode = icode.strip()
		ind = (chain_id, resseq, icode)
		# print('Searching ')
		# print(ind)
		if ind in index_list:
			# if ind not in residue_coord_dict:
			#	residue_coord_dict[ind] = []
			loop_specific_residue_list.append(residue)

			# for atom in r:
			#   point = r[atom.name].get_vector()

			#   min_x, min_y, min_z = min(point[0], min_x), min(point[1], min_y), min(point[2], min_z)
			#   max_x, max_y, max_z = max(point[0], max_x), max(point[1], max_y), max(point[2], max_z)
			#   # residue_coord_dict[ind].append(r[atom.name].get_vector())

	return loop_specific_residue_list

def load_pdb_chainwise_residue_data_for_all_loops(pdb_chainwise_loops, directories, utilize_pickle, pickle_fname='pdb_chainwise_residue_data.pkl3'):
	pickle_fname = os.path.join(directories.pickles_dir, pickle_fname)
	# pickle_fname = 'load_pdb_chainwise_residue_data_for_all_loops-all_HL_IL.pkl3'
	pdb_chainwise_residue_data = {}

	print('')
	if utilize_pickle and os.path.isfile(pickle_fname):
		logger.info('Loading PDB-chain-wise residue data from a previously generated pickle file.')
		pf = open(pickle_fname, 'rb')
		pdb_chainwise_residue_data = pickle.load(pf)
		pf.close()
	
	else:
		# pdb_chainwise_residue_data = {}

		logger.info('Generating PDB-chain-wise residue data.')
		start_time = time.time()

		pdb_count = len(pdb_chainwise_loops)
		for i, pdb_id in enumerate(pdb_chainwise_loops):
			print('Processing ' + str(pdb_id) + ' (' + str(i+1) + '/' + str(pdb_count) + ')')
			pdb_chainwise_residue_data[pdb_id] = {}

			chain_cnt = len(pdb_chainwise_loops[pdb_id])
			print(str(chain_cnt) + ' chain(s) to process')

			pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
			parser = FastMMCIFParser(QUIET=True)
			pdb_structure = parser.get_structure('struct', pdb_fname)

			for j, chain_id in enumerate(pdb_chainwise_loops[pdb_id]):
				print('processing chain ' + str(chain_id) + ' (' + str(j+1) + '/' + str(chain_cnt) + ')')
				pdb_chainwise_residue_data[pdb_id][chain_id] = {}

				loops = pdb_chainwise_loops[pdb_id][chain_id]
				print(str(len(loops)) + ' loops')
				chain_structure = pdb_structure[0][chain_id]

				residues, min_coords, max_coords = load_chain_data(chain_structure)

				loop_info_dict = {}
				for loop in tqdm(loops, ncols=75):
					loop_residues = get_loop_specific_residues(loop, residues, chain_id, directories)
					loop_coords = extract_motif_coordinates(loop_residues)
					loop_info_dict[loop] = loop_coords
				# print(loop_info_dict)
				# sys.exit()
				pdb_chainwise_residue_data[pdb_id][chain_id] = (loop_info_dict, min_coords, max_coords)

				# loop_info_dict = {}
				# for loop in loops:
				# 	# box_coord, loop_all_res_list = get_loop_bounding_box_coord(loop, directories, chain_structure)
				# 	loop_info_dict[loop] = get_residue_list_of_a_loop(loop, directories, chain_structure)
				# 	# loop_info_dict[loop] = (box_coord, loop_all_res_list)

				# for loop in tqdm(loops, ncols=75):

		logger.info('Done')
		logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

		pf = open(pickle_fname, 'wb')
		pickle.dump(pdb_chainwise_residue_data, pf)

	return pdb_chainwise_residue_data

def load_spatial_proximity_data(pdb_chainwise_loops, directories, atom_set_choice, utilize_pickle, pickle_fname='filtered_spatial_proximity_data.pkl3'):
	# pickle_fname = 'filtered_spatial_proximity_data_all_pair.pkl3'
	pickle_fname = os.path.join(directories.pickles_dir, pickle_fname)
	# pickle_fname = 'filtered_spatial_proximity_data_all_pair-all_HL_IL.pkl3'
	spatial_proximity_data = {}

	if utilize_pickle and os.path.isfile(pickle_fname):
		logger.info('Loading spatial proximity data from a previously generated pickle file.')

		pf = open(pickle_fname, 'rb')
		spatial_proximity_data = pickle.load(pf)
		pf.close()
		# print(spatial_proximity_data)
		# sys.exit()
	else:
		# prev_pdb_id = None
		# pdb_structure = None
		# spatial_proximity_data = {}
		# filtered_spatial_proximity_data = {}
		# graphs = []
		logger.info('Generating spatial proximity data.')
		start_time = time.time()

		pdb_count = len(pdb_chainwise_loops)
		for i, pdb_id in enumerate(pdb_chainwise_loops):
			# pdb_id = '8BUU'

			print('Processing ' + str(pdb_id) + ' (' + str(i+1) + '/' + str(pdb_count) + ')')
			chain_cnt = len(pdb_chainwise_loops[pdb_id])
			print(str(chain_cnt) + ' chain(s) to process')

			spatial_proximity_data[pdb_id] = {}
			# filtered_spatial_proximity_data[pdb_id] = {}

			pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
			parser = FastMMCIFParser(QUIET=True)
			pdb_structure = parser.get_structure('struct', pdb_fname)

			for j, chain_id in enumerate(pdb_chainwise_loops[pdb_id]):
				# chain_id = 'A'
				print('processing chain ' + str(chain_id) + ' (' + str(j+1) + '/' + str(chain_cnt) + ')')

				spatial_proximity_data[pdb_id][chain_id] = {}
				# filtered_spatial_proximity_data[pdb_id][chain_id] = {}

				loops = pdb_chainwise_loops[pdb_id][chain_id]
				print(str(len(loops)) + ' loops')

				chain_structure = pdb_structure[0][chain_id]

				loop_info_dict = {}
				for loop in loops:
					# box_coord, loop_all_res_list = get_loop_bounding_box_coord(loop, directories, chain_structure)
					loop_info_dict[loop] = get_residue_list_of_a_loop(loop, directories, chain_structure)
					# loop_info_dict[loop] = (box_coord, loop_all_res_list)

				for loop in tqdm(loops, ncols=75):

					spatial_proximity_data[pdb_id][chain_id][loop] = get_spatial_proximity_data_of_a_loop(loop, loop_info_dict, loops, atom_set_choice)

					# # 1st step: narrowing the search space
					# nearest_loop_list = get_nearest_loops_from_bounding_box(loop, loop_info_dict, loops, distance_threshold_to_be_nearest_residue)

					# if len(nearest_loop_list) > 0:
					# 	spatial_proximity_data[pdb_id][chain_id][loop] = nearest_loop_list

					# 	filtered_nearest_loop_list = get_filtered_nearest_loops_data(loop, loop_info_dict, loops, nearest_loop_list, distance_threshold_to_be_nearest_residue)

					# 	if len(filtered_nearest_loop_list) > 0:
					# 		filtered_spatial_proximity_data[pdb_id][chain_id][loop] = filtered_nearest_loop_list


				# if len(spatial_proximity_data[pdb_id][chain_id]) > 0:			
				# 	print_nearest_loop_info(pdb_id, chain_id, spatial_proximity_data[pdb_id][chain_id], families)
				# 	# nearest_residue_data, pdb_structure = get_nearest_residue_data_from_bounding_box(pdb_id, chain_id, box_coord, pdb_structure, distance_threshold_to_be_nearest_residue)

				# 	if len(filtered_spatial_proximity_data[pdb_id][chain_id]) > 0:
				# 		print('Precise:')
				# 		print_nearest_loop_info(pdb_id, chain_id, filtered_spatial_proximity_data[pdb_id][chain_id], families)

						
				# g = create_chain_graph(loops, families, pdb_id, chain_id, filtered_spatial_proximity_data)
				# graphs.append(g)
						
		logger.info('Done')
		logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

		pf = open(pickle_fname, 'wb')
		pickle.dump(spatial_proximity_data, pf)
		# sys.exit()

	return spatial_proximity_data

if __name__ == '__main__':
	main()
