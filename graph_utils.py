import networkx as nx
import numpy

from functools import reduce
# from Bio import SeqIO, pairwise2
from Bio.PDB import *

import matplotlib.pyplot as plt

from utils import *

def load_pdb_res_map(directories, chain):
	"""load sequence index->pdb index"""
	"""{ref_index: (chain_id, pdb_index)}"""
	ret = {}
	# map_dir = '../nrPDBs_Old' # the directory for the mapping data
	fp = open(os.path.join(directories.pdb_fasta_mapping_dir, get_modified_chain_id_if_any_lowercase_letter(chain)+'.rmsx.nch'))
	for line in fp.readlines():
		decom = line.strip().split('\t')
		##################### for PDB #####################
		# if decom[0][0] == "'":
		#	 ret[int(decom[1])] = (decom[0][1], decom[0][3:].replace('.', ''))
		# else:
		#	 ret[int(decom[1])] = (decom[0][0], decom[0][1:].replace('.', ''))
		##################### for PDB #####################
		##################### for PDBx ####################
		if decom[0][0] == "'":
			chain_id = decom[0][1:].strip().split("'")[0]
			i = len(chain_id)+2
		else:
			chain_id = re.split('-?(\d+)',decom[0])[0]
			i = len(chain_id)

		if decom[0][-1].isalpha():
			icode = decom[0][-1]
			j = len(decom[0])-2
		else:
			icode = ''
			j = len(decom[0])

		seqnum = decom[0][i:j]
		ret[int(decom[1])] = (chain_id, seqnum, icode)
		##################### for PDBx ####################
	return ret

def pdb_pos_map(pdb_res_map, m):
	"""position in pdb index alignment"""
	ret = []
	for i in m:
		if i in pdb_res_map:
			ret.append(pdb_res_map[i])
		# if 
		else:
			# ret.append('na')
			ret.append(('', '', ''))
			logger.warning('!!!!!!!!!!!!!!!!!!!!!ALERT: APPENDING EMPTY TUPLE (NA) !!!!!!!!!!!!!!!!!!!!')

	return ret

def centroid(coord_list):
	if len(coord_list) > 0:
		return list(map(lambda z: 1.*z/len(coord_list), reduce(lambda x, y: (x[0]+y[0], x[1]+y[1], x[2]+y[2]), coord_list)))
	return None

def get_atom_coordinate(pdb_fn, residue_list):

	backbone_atoms, sugar_atoms = get_backbone_and_sugar_atoms()
	pdb_id = os.path.basename(pdb_fn)[:4]

	parser = FastMMCIFParser()
	structure = parser.get_structure('struct', pdb_fn)

	backbone = {}
	sugar = {}

	for chain_id, index, icd in residue_list:
		# if chain_id == 'n' and index == 'a':
		if chain_id == '':
			continue
		chain = structure[0][chain_id]
		residues = chain.get_residues()

		my_residues = {}

		for r in residues:
			hetflag, resseq, icode = r.get_id()
			my_residues[(resseq, icode)] = r

		i = int(index)
		icode = icd if len(icd) > 0 else ' '

		if (i, icode) not in my_residues:
			# ret.append(0)
			backbone[(pdb_id, chain_id, index, icd)] = [0., 0., 0.]
			sugar[(pdb_id, chain_id, index, icd)] = [0., 0., 0.]
		else:
			atom_coord = []
			for atom in backbone_atoms:
				if atom in my_residues[(i, icode)]:
					atom_coord.append(my_residues[(i, icode)][atom].get_vector())

			backbone[(pdb_id, chain_id, index, icd)] = centroid(atom_coord)

			atom_coord = []
			for atom in sugar_atoms:
				if atom in my_residues[(i, icode)]:
					atom_coord.append(my_residues[(i, icode)][atom].get_vector())

			sugar[(pdb_id, chain_id, index, icd)] = centroid(atom_coord)

	return backbone, sugar, structure

# def get_pdb_index_list(directories, lp):
# 	pdb_chain, regions = lp.split(':')
# 	# segments = regions.strip().split('_')
# 	# index_list = []
# 	r = list(map(lambda x: x.split('-'), regions.split('_')))
# 	index_list = reduce(lambda y, z: y+z, list(map(lambda x: list(range(int(x[0]), int(x[1])+1)), r)))
# 	pdb_res_map = load_pdb_res_map(directories, pdb_chain)
# 	pdb_index_list = pdb_pos_map(pdb_res_map, index_list)

# 	return pdb_index_list

# def get_residue_list_of_a_loop(loop, directories, chain_structure):
#     min_x = min_y = min_z = 99999
#     max_x = max_y = max_z = -99999

#     loop_all_res_list = []

#     # loop_cif_extension = 5
#     # partial_pdbx_dir = os.path.join(directories.data_dir, 'pdbx_extracted_ext' + str(loop_cif_extension))
#     # print(directories.partial_pdb_dir)
#     # directories.partial_pdb_dir = directories.partial_pdb_dir.replace('*', str(loop_cif_extension))
#     # print(directories.partial_pdb_dir)
#     # sys.exit()

#     pdb_chain = loop.strip().split(':')[0]
#     pdb_id, chain_id = pdb_chain.strip().split('_')

#     # pdb_fn = os.path.join(directories.partial_pdb_dir, loop.replace(':', '_') + '.cif')
#     index_list = get_pdb_index_list(directories, loop)
#     # print(index_list)
#     index_list = list(map(lambda x: (x[0], int(x[1]), x[2]), index_list))

#     # parser = FastMMCIFParser(QUIET=True)
#     # structure = parser.get_structure('struct', pdb_fn)

#     # chain = structure[0][chain_id]
#     residues = chain_structure.get_residues()

#     # print(index_list)
#     # residue_coord_dict = {}
#     for r in residues:
#         # print(r.get_id())
#         hetflag, resseq, icode = r.get_id()
#         icode = icode.strip()
#         ind = (chain_id, resseq, icode)
#         # print('Searching ')
#         # print(ind)
#         if ind in index_list:
#             # if ind not in residue_coord_dict:
#             #    residue_coord_dict[ind] = []
#             loop_all_res_list.append(r)

#             # for atom in r:
#             #   point = r[atom.name].get_vector()

#             #   min_x, min_y, min_z = min(point[0], min_x), min(point[1], min_y), min(point[2], min_z)
#             #   max_x, max_y, max_z = max(point[0], max_x), max(point[1], max_y), max(point[2], max_z)
#             #   # residue_coord_dict[ind].append(r[atom.name].get_vector())

#     return loop_all_res_list

def get_simplified_index(index_dict, index):
	if index not in index_dict:
		index_dict[index] = len(index_dict)

	return index_dict[index], index_dict

def distance_3d(x1, y1, z1, x2, y2, z2):
    d = math.sqrt(math.pow(x2 - x1, 2) + math.pow(y2 - y1, 2) + math.pow(z2 - z1, 2)* 1.0)
    # print("Distance is ")
    # print(d)
    return d

def distance_3d_from_2_points(point1 , point2):
	x1, y1, z1 = point1
	x2, y2, z2 = point2
	return distance_3d(x1, y1, z1, x2, y2, z2)

def create_motif_graph(loop, pdb_chain, loop_data, coord_backbone, coord_sugar, pdb_pm):
	# pass
	pdb_id, chain_id = pdb_chain.strip().split('_')
	joined_sequence, bps, bp_cnt, stks, stk_cnt, break_points = loop_data

	seq_ind_to_pdb_ind = dict(zip(list(range(len(loop_data[0]))), pdb_pm))
	coord_dict = {}
	for i in list(range(len(joined_sequence))):
		coords = []
		chain_id, seqnum, icode = seq_ind_to_pdb_ind[i]
		if (pdb_id, chain_id, seqnum, icode) in coord_backbone:
			coords.append(coord_backbone[(pdb_id, chain_id, seqnum, icode)])
		if (pdb_id, chain_id, seqnum, icode) in coord_sugar:
			coords.append(coord_sugar[(pdb_id, chain_id, seqnum, icode)])

		centroid = numpy.sum(coords, axis=0)/float(len(coords))
		# coord_list.append(centroid)
		coord_dict[i] = centroid

	G = nx.Graph()
	# node_attributes = {}
	# edge_attributes = {}
	for i, res in enumerate(joined_sequence):
		# G.add_node(i, residue_name=res, residue_id=i, xyz=coord_dict[i], is_terminal=(i==0 or i-1 in break_points))
		# node_attributes[i] = {}
		# node_attributes[i]['residue_name'] = res
		# node_attributes[i]['xyz'] = coord_dict[i]
		# node_attributes[i]['is_terminal'] = (i==0 or i-1 in break_points)
		# G.add_node(i, n_data=numpy.array([res, coord_dict[i], (i==0 or i-1 in break_points)]))
		# G.add_node(i, residue_name=res)
		G.add_node(i, residue_name=res, is_terminal=(i==0 or i-1 in break_points))

	k = 0
	for br_p in break_points:
		for i in list(range(k, br_p-1)):
			# G.add_edge(i, i+1, bond_type='backbone', distance=distance_3d_from_2_points(coord_dict[i], coord_dict[i+1]))
			# edge_attributes[(i, i+1)] = {}
			# edge_attributes[(i, i+1)]['bond_type'] = 'backbone'
			# edge_attributes[(i, i+1)]['hashed_intera'] = None
			# edge_attributes[(i, i+1)]['distance'] = distance_3d_from_2_points(coord_dict[i], coord_dict[i+1])
			# G.add_edge(i, i+1, e_data=numpy.array(['backbone', None, distance_3d_from_2_points(coord_dict[i], coord_dict[i+1])]))
			# G.add_edge(i, i+1, intera='backbone', distance=distance_3d_from_2_points(coord_dict[i], coord_dict[i+1]))
			G.add_edge(i, i+1, intera=0, distance=distance_3d_from_2_points(coord_dict[i], coord_dict[i+1]))
		k = br_p

	for s in bps:
		for e in bps[s]:
			for (ss, ee, ntd_pair, interaction) in bps[s][e]:
				# G.add_edge(ss, ee, paired_base1=ntd_pair[0], paired_base2=ntd_pair[1], interaction=interaction, hashed_intera=hash_basepair(interaction, ntd_pair), distance=distance_3d_from_2_points(coord_dict[ss], coord_dict[ee]))
				# edge_attributes[(ss, ee)] = {}
				# edge_attributes[(ss, ee)]['bond_type'] = 'base-pair'
				# edge_attributes[(ss, ee)]['hashed_intera'] = hash_basepair(interaction, ntd_pair)
				# edge_attributes[(ss, ee)]['distance'] = distance_3d_from_2_points(coord_dict[ss], coord_dict[ee])
				# G.add_edge(ss, ee, e_data=numpy.array(['basepair', hash_basepair(interaction, ntd_pair), distance_3d_from_2_points(coord_dict[ss], coord_dict[ee])]))
				G.add_edge(ss, ee, intera=hash_basepair(interaction, ntd_pair), distance=distance_3d_from_2_points(coord_dict[ss], coord_dict[ee]))

	# G.graph['sequence'] = joined_sequence

	# nx.set_node_attributes(G, node_attributes, 'n_attributes')
	# nx.set_edge_attributes(G, edge_attributes, 'e_attributes')

	nx.draw_circular(G, with_labels = True)
	plt.savefig(loop + ".png")
	plt.clf()
	# sys.exit()

	return G

def generate_graphs_from_loops(directories, loops):
	graphs = []
	# motifs = []
	not_found_list = []
	for loop in loops:
		# print(loop)
		pdb_chain, regions = loop.strip().split(':')
		file_name = os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf')
		if not os.path.exists(file_name):
			not_found_list.append(loop)
			continue
		loop_data = load_loop_data(os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf'), False)
		# print(loop_data)
		pdb_pm = get_pdb_index_list(directories, loop)
		# print(dict(zip(pdb_pm, list(range(len(loop_data[0]))))))
		coord_backbone, coord_sugar, pdb_structure = get_atom_coordinate(os.path.join(directories.partial_pdb_dir, loop.replace(':', '_') + '.cif'), pdb_pm)

		# index_dict = {}
		# for ind in sorted(coord_backbone):
		# 	simplified_index, index_dict = get_simplified_index(index_dict, ind)
		# print(coord_backbone)
		# print(index_dict)
		# sys.exit()

		g = create_motif_graph(loop, pdb_chain, loop_data, coord_backbone, coord_sugar, pdb_pm)
		graphs.append(g)

	return graphs