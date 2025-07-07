import os
import sys

from graph_utils import *
from functools import reduce
from tqdm import tqdm
from utils import *

def get_pdb_index_list(directories, lp):
	pdb_chain, regions = lp.split(':')
	# segments = regions.strip().split('_')
	# index_list = []
	r = list(map(lambda x: x.split('-'), regions.split('_')))
	index_list = reduce(lambda y, z: y+z, list(map(lambda x: list(range(int(x[0]), int(x[1])+1)), r)))
	pdb_res_map = load_pdb_res_map(directories, pdb_chain)
	pdb_index_list = pdb_pos_map(pdb_res_map, index_list)

	return pdb_index_list

def get_residue_list_of_a_loop(loop, directories, chain_structure):
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
            #    residue_coord_dict[ind] = []
            loop_all_res_list.append(r)

            # for atom in r:
            #   point = r[atom.name].get_vector()

            #   min_x, min_y, min_z = min(point[0], min_x), min(point[1], min_y), min(point[2], min_z)
            #   max_x, max_y, max_z = max(point[0], max_x), max(point[1], max_y), max(point[2], max_z)
            #   # residue_coord_dict[ind].append(r[atom.name].get_vector())

    return loop_all_res_list


def get_spatial_proximity_data_of_a_loop_from_coords(loop1, loop_info_dict, loops):
	spatial_proximity_data_of_a_loop = {}

	coord_list1 = loop_info_dict[loop1]

	for loop2 in loops:
		min_distance = 0
		if loop1 != loop2:
			coord_list2 = loop_info_dict[loop2]
			min_distance = find_minimum_distance_between_two_coord_lists(coord_list1, coord_list2)

		spatial_proximity_data_of_a_loop[loop2] = min_distance

	return spatial_proximity_data_of_a_loop

def find_minimum_distance_between_two_coord_lists(coord_list1, coord_list2):
	min_distance = 999999

	for coord1 in coord_list1:
		for coord2 in coord_list2:
			distance = find_minimum_distance_between_two_coords(coord1, coord2)
			min_distance = min(min_distance, distance)

	return min_distance

def find_minimum_distance_between_two_coords(coord1, coord2):
	min_distance = 999999

	for point1 in coord1:
		# point1 = res1[atom1.name].get_vector()
		for point2 in coord2:
			# point2 = res2[atom2.name].get_vector()
			distance = distance_between_points(point1, point2)
			min_distance = min(min_distance, distance)

	return min_distance




def get_spatial_proximity_data_of_a_loop(loop1, loop_info_dict, loops, atom_set_choice):
	spatial_proximity_data_of_a_loop = {}

	residue_list1 = loop_info_dict[loop1]

	for loop2 in loops:
		min_distance = 0
		if loop1 != loop2:
			residue_list2 = loop_info_dict[loop2]
			min_distance = find_minimum_distance_between_two_residue_lists(residue_list1, residue_list2, atom_set_choice)

		spatial_proximity_data_of_a_loop[loop2] = min_distance

	return spatial_proximity_data_of_a_loop

def find_minimum_distance_between_two_residue_lists(residue_list1, residue_list2, atom_set_choice):
	min_distance = 999999

	for res1 in residue_list1:
		for res2 in residue_list2:
			distance = find_minimum_distance_between_two_residues(res1, res2, atom_set_choice)
			# distance0 = find_minimum_distance_between_two_residues(res1, res2, 0)
			# distance1 = find_minimum_distance_between_two_residues(res1, res2, 1)
			# distance2 = find_minimum_distance_between_two_residues(res1, res2, 2)
			# print(distance, distance0, distance1, distance2)
			min_distance = min(min_distance, distance)

	return min_distance

def find_minimum_distance_between_two_residues(res1, res2, atom_set_choice=0):
	min_distance = 999999

	backbone_atoms, sugar_atoms = get_backbone_and_sugar_atoms()

	for atom1 in res1:
		if atom_set_choice == 1 and not (atom1.name in backbone_atoms or atom1.name in sugar_atoms):
			continue
		elif atom_set_choice == 2 and atom1.name not in backbone_atoms:
			continue

		point1 = res1[atom1.name].get_vector()
		for atom2 in res2:
			if atom_set_choice == 1 and not (atom2.name in backbone_atoms or atom2.name in sugar_atoms):
				continue
			elif atom_set_choice == 2 and atom2.name not in backbone_atoms:
				continue

			point2 = res2[atom2.name].get_vector()
			distance = distance_between_points(point1, point2)
			min_distance = min(min_distance, distance)

	return min_distance

def distance_between_points(p1, p2):
	x1, y1, z1 = p1
	x2, y2, z2 = p2
	return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def get_nearest_loop_list(pdb_id, chain_id, pdb_chainwise_loops, loop1, spatial_proximity_data, distance_threshold):
	nearest_loop_list = []
	for loop2 in pdb_chainwise_loops[pdb_id][chain_id]:
		if loop1 == loop2:
			continue
		if spatial_proximity_data[pdb_id][chain_id][loop1][loop2] <= distance_threshold:
			nearest_loop_list.append(loop2)

	return nearest_loop_list

def get_nearest_loop_list_v2(pdb_id, chain_id, pdb_chainwise_loops, loop1, spatial_proximity_data, distance_threshold):
	nearest_loop_list = []
	residuewise_coords = pdb_chainwise_loops[pdb_id][chain_id]
	for loop2, coord in residuewise_coords:
		if loop1 == loop2:
			continue
		if spatial_proximity_data[pdb_id][chain_id][loop1][loop2] <= distance_threshold:
			nearest_loop_list.append(loop2)

	return nearest_loop_list

# def get_family_id(loop, families):
# 	for family_id in families:
# 		if loop in families[family_id]:
# 			return family_id
# 	return 'N/A'

def get_motif_module_frequencies(spatial_proximity_data, families):
	family_group_dict = {}
	for pdb_id in spatial_proximity_data:
		for chain_id in spatial_proximity_data[pdb_id]:
			for loop1 in spatial_proximity_data[pdb_id][chain_id]:
				if len(spatial_proximity_data[pdb_id][chain_id][loop1]) == 0:
					continue

				family_group = []
				loop_group = []
				
				family_group.append(get_family_id(loop1, families))
				loop_group.append(loop1)

				for loop2 in spatial_proximity_data[pdb_id][chain_id][loop1]:
					family_group.append(get_family_id(loop2, families))
					loop_group.append(loop2)

				sorted_family_group = sorted(family_group)
				sorted_loop_group = sorted(loop_group)

				sorted_family_group_tuple = tuple(sorted_family_group)

				if sorted_family_group_tuple not in family_group_dict:
					family_group_dict[sorted_family_group_tuple] = []
					family_group_dict[sorted_family_group_tuple].append(sorted_loop_group)
				else:
					if sorted_loop_group not in family_group_dict[sorted_family_group_tuple]:
						family_group_dict[sorted_family_group_tuple].append(sorted_loop_group)

	family_group_dict = dict(sorted(family_group_dict.items(), key=lambda item: len(item[1]), reverse=True))

	return family_group_dict

def update_dict_with_ids(family_group_dict):
	family_group_dict_with_ids = {}

	family_group_dict = dict(sorted(family_group_dict.items(), key=lambda item: len(item[1]), reverse=True))
	for module_id, family_group_tuple in enumerate(family_group_dict):
		family_group_dict_with_ids[(module_id+1, family_group_tuple)] = []
		for loop_group_id, loop_list in enumerate(family_group_dict[family_group_tuple]):
			family_group_dict_with_ids[(module_id+1, family_group_tuple)].append((loop_group_id+1, loop_list))

	return family_group_dict_with_ids

def generate_motif_module_info_from_pdb_chainwise_data(pdb_chainwise_loops, families, distance_threshold_to_be_nearest_residue):
	pdb_count = len(pdb_chainwise_loops)
	spatial_proximity_data = {}
	for i, pdb_id in enumerate(pdb_chainwise_loops):
		print('Processing ' + str(pdb_id) + ' (' + str(i+1) + '/' + str(pdb_count) + ')')
		chain_cnt = len(pdb_chainwise_loops[pdb_id])

		spatial_proximity_data[pdb_id] = {}
		for j, chain_id in enumerate(pdb_chainwise_loops[pdb_id]):
			print('processing chain ' + str(chain_id) + ' (' + str(j+1) + '/' + str(chain_cnt) + ')')

			spatial_proximity_data[pdb_id][chain_id] = {}
			loops_with_residues = pdb_chainwise_loops[pdb_id][chain_id]
			loop_info_dict = {}
			loops = []
			for loop, coords in loops_with_residues:
				loops.append(loop)
				loop_info_dict[loop] = coords

			for loop, coords in tqdm(loops_with_residues, ncols=75):
			# for loop, residues in loops_with_residues:
				spatial_proximity_data[pdb_id][chain_id][loop] = get_spatial_proximity_data_of_a_loop_from_coords(loop, loop_info_dict, loops)

	nearest_loop_data = {}
	for i, pdb_id in enumerate(pdb_chainwise_loops):
		nearest_loop_data[pdb_id] = {}
		for j, chain_id in enumerate(pdb_chainwise_loops[pdb_id]):
			nearest_loop_data[pdb_id][chain_id] = {}
			residuewise_coords = pdb_chainwise_loops[pdb_id][chain_id]
			for loop, coord in residuewise_coords:
				if loop in nearest_loop_data[pdb_id][chain_id]:
					nearest_loop_data[pdb_id][chain_id][loop] += get_nearest_loop_list_v2(pdb_id, chain_id, pdb_chainwise_loops, loop, spatial_proximity_data, distance_threshold_to_be_nearest_residue)
				else:
					nearest_loop_data[pdb_id][chain_id][loop] = get_nearest_loop_list_v2(pdb_id, chain_id, pdb_chainwise_loops, loop, spatial_proximity_data, distance_threshold_to_be_nearest_residue)

	family_group_dict = get_motif_module_frequencies(nearest_loop_data, families)
	# family_group_dict_with_ids = update_dict_with_ids(family_group_dict)

	# print(family_group_dict_with_ids.keys())
	# sys.exit()
	return family_group_dict