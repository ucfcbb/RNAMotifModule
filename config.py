import sys
import os
import logging

from classes import *
from logger import *

load_search_input_from_chain = True

def get_env():
	return 'local'

def get_urls():
	pdbx_url = 'https://files.rcsb.org/download/'
	# fasta_url = 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?compressionType=uncompressed&structureIdList='
	fasta_url = 'https://www.rcsb.org/fasta/entry/'
	fr3d_url = "http://rna.bgsu.edu/rna3dhub/pdb/XXXX/interactions/fr3d/all/csv"
	dssr_url = 'http://skmatic.x3dna.org/pdb/XXXX/XXXX.out'

	return pdbx_url, fasta_url, fr3d_url, dssr_url

def get_misc_params():
	annotation_source = 'merged'
	input_index_type = 'pdb'

	items_per_chunk = 5
	mp_number_of_process = 16
	content_download_attempts = 7

	return input_index_type, annotation_source, items_per_chunk, mp_number_of_process, content_download_attempts

def get_scores(penalties, bonuses):

	hbond_match_base = 2.0
	weight_isosteric = 2.0
	weight_nonadjacent_stacking = 1.0
	weight_adjacent_stacking = 0.2
	weight_sequence = 0.1
	
	asym_nuc = -0.2
	asym_loop = -2.0
	cons_stacking = 0.5
	stack_to_bulge = -1.5
	stack_to_internal_asym = -1.0
	stack_to_internal_sym = -0.5

	scores = Scores(penalties, bonuses, hbond_match_base, weight_isosteric, weight_nonadjacent_stacking, weight_adjacent_stacking, weight_sequence, asym_nuc, asym_loop, cons_stacking, stack_to_bulge, stack_to_internal_asym, stack_to_internal_sym)

	return scores

def get_mat_file_names():
	isosteric_file_name = 'data/mat/iso.mat'
	stacking_file_name = 'data/mat/stk.mat'
	nucleotide_file_name = 'data/mat/nuc.mat'

	return isosteric_file_name, stacking_file_name, nucleotide_file_name

def get_base_dir_names():
	root_dir = os.getcwd()
	base_dir = root_dir
	data_dir = os.path.join(root_dir, 'data/')

	pdbx_dir = os.path.join(data_dir, 'pdbx/')
	fasta_dir = os.path.join(data_dir, 'fasta/')
	pdb_fasta_mapping_dir = os.path.join(data_dir, 'pdb_fasta/')
	annotation_dir = os.path.join(data_dir, 'annotation/')
	loop_dir = os.path.join(data_dir, 'loops/')
	partial_pdb_dir = os.path.join(data_dir, 'pdbx_extracted_ext*/')
	output_dir = os.path.join(root_dir, 'output/')
	pickles_dir = os.path.join(data_dir, 'pickles/')
	temp_dir = os.path.join(root_dir, 'temp/')

	lib_dir = os.path.join(root_dir, 'my_lib/')
	dssr_dir = os.path.join(lib_dir, 'DSSR/')

	return Directories(base_dir, data_dir, pdbx_dir, fasta_dir, pdb_fasta_mapping_dir, annotation_dir, loop_dir, partial_pdb_dir, output_dir, temp_dir, lib_dir, dssr_dir, pickles_dir)

def get_scanx_dir():
	directories = get_base_dir_names()
	# return os.path.join(directories.lib_dir, 'RNAMotifScanX-release')
	return os.path.join(directories.lib_dir, 'RNAMotifScanX-release')

# def get_include_stackings_flag_for_profile_generation():
# 	return False

# def get_include_stackings_flag_for_searching():
# 	return True

def get_include_stackings_flag(mode):
	if mode.lower() == 'search':
		return False
	else:
		return False

def validate_params(args):

	mode = str(args.m)
	if len(mode) > 0:
		if mode.lower() != 'search':
			# print('Please provide "search" as the value of -m parameter to search motifs.')
			logger.error('Please provide \'search\' as the value of -m parameter to search motifs. Please disregard this parameter to use this tool in profile generation mode.')
			sys.exit()
		else:
			profile_fname = str(args.p)
			if len(profile_fname) == 0:
				# print('Please provide a profile (.pfl) file when using in "search" mode.')
				logger.error('Please provide a valid profile (.pfl) file with -p when using in \'search\' mode.')
				sys.exit()
			elif os.path.basename(profile_fname).strip().split('.')[-1] != 'pfl':
				logger.error('Please provide a profile (.pfl) file with -p when using in \'search\' mode.')
				sys.exit()
			elif os.path.isfile(profile_fname) == False:
				logger.error('File not exists! Please provide a valid profile (.pfl) file with -p when using in \'search\' mode.')
				sys.exit()

			# include_scanx_score = args.x
			# scanx_structure_file = str(args.xs)
			# if include_scanx_score == True:
			# 	if len(scanx_structure_file) == 0:
			# 		logger.error('Please provide a valid scanx structure (.struct) file with -xs to include scanx alignment score in \'search\' mode.')
			# 		sys.exit()
			# 	elif os.path.basename(scanx_structure_file).strip().split('.')[-1] != 'struct':
			# 		logger.error('Please provide a scanx structure (.struct) file with -xs to include scanx alignment score in \'search\' mode.')
			# 		sys.exit()
			# 	elif os.path.isfile(scanx_structure_file) == False:
			# 		logger.error('File not exists! Please provide a valid scanx structure (.struct) file with -xs to include scanx alignment score in \'search\' mode.')
			# 		sys.exit()

			pdb_chain_to_search = args.c
			if load_search_input_from_chain == True:
				if len(pdb_chain_to_search) == 0:
					logger.error('Please provide a valid chain in <PDB_ID>_<Chain_ID> format with -c when using in \'search\' mode.')
					sys.exit()
				elif '_' not in pdb_chain_to_search:
					logger.error('Please provide a valid chain in <PDB_ID>_<Chain_ID> format with -c when using in \'search\' mode.')
					sys.exit()