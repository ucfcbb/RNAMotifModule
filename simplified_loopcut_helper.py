import os
import sys
import time
import glob
import multiprocessing as mp
import numpy
from cif import *
# from Bio import SeqIO, pairwise2
from Bio import SeqIO, Align
from functools import reduce
import tempfile

# from ann_merge_helper import *
from ann_parser import *
# from classes import *
from utils import *

cWW_stack_threshold = 1
include_multi_chain_loops = False
wc_count_stat = {}
is_median = False
truncate_percentage = 0
extension = ".smf"

def create_directory(dir_to_create):
    if not os.path.exists(dir_to_create):
        os.makedirs(dir_to_create)

def write_loop_to_directory(input_dir, loops, output_dir, extension):
    for loop in loops:
        fr = open(os.path.join(input_dir, loop + extension))
        fw = open(os.path.join(output_dir, loop + extension), "w")
        for line in fr.readlines():
            fw.write(line)
        fr.close()
        fw.close()

def get_loop_length(segments):
    # will not work properly for pdb index
    loop_length = 0
    for segment in segments:
        pcs = segment.split("-")
        if pcs[0][-1].isalpha():
            s = int(pcs[0].strip().split(".")[0].strip())
        else:
            s = int(pcs[0].strip())
        if pcs[1][-1].isalpha():
            e = int(pcs[1].strip().split(".")[0].strip())
        else:
            e = int(pcs[1].strip())
        # s = int(pieces[0])
        # e = int(pieces[1])
        loop_length += (e-s+1)
    return loop_length

def get_modified_chain_id_if_any_lowercase_letter(chain_id):
    if chain_id != chain_id.upper():
        chain_id = chain_id + '_'
    return chain_id

def load_fasta_seq(pdb_id, chains, fasta_dir):
    fasta_seq_dict = {}
    fasta_fn = os.path.join(fasta_dir, pdb_id + '.fasta')
    for record in SeqIO.parse(fasta_fn, 'fasta'):
        # fasta_seq_dict[record.id.strip().split('|')[0].strip().split(':')[1]] = str(record.seq)
        # chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1].strip().split(',')
        # for chain_id in chain_ids:
        #     fasta_seq_dict[chain_id] = str(record.seq)
        chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1:]
        for chain_id in chain_ids:
            chain_id = chain_id.strip().strip(',')
            if '[' in chain_id:
                continue
                # chain_id = chain_id.split('[')[0].strip()
            elif ']' in chain_id:
                chain_id = chain_id.split(']')[0].strip()
            
            fasta_seq_dict[chain_id] = str(record.seq)

    return fasta_seq_dict

def base_abbreviation(fn):
    """get the abbreviation of the modified residues from the 3DNA baselist"""
    ret = {}
    # fn = os.path.join(os.path.dirname(os.path.abspath( __file__ )), fn)
    fp = open(fn)
    for line in fp.readlines():
        line = line.rstrip()
        if line.startswith("#") or len(line) == 0:
            continue
        else:
            three_letter = line[:3].strip()
            one_letter = line[8]
            if one_letter == "T" or one_letter == "t":
                one_letter = "U"
            ret[three_letter] = one_letter.upper()
    fp.close()
    return ret

def amino_acid_collection(fn):
    """get the list of the amino acids from the file"""
    ret = {}
    # fn = os.path.join(os.path.dirname(os.path.abspath( __file__ )), fn)
    fp = open(fn)
    for line in fp.readlines():
        line = line.rstrip()
        if line.startswith("#") or len(line) == 0:
            continue
        else:
            three_letter = line[:3].strip().upper()
            one_letter = line[4]
            ret[three_letter] = one_letter.upper()
    fp.close()
    return ret

def write_invalid_mapping_file(pdb_id, chain_id):
    chain_id = get_modified_chain_id_if_any_lowercase_letter(chain_id)
    mapping_fname = os.path.join(pdb_fasta_mapping_dir, pdb_id + '_' + chain_id + '_invalid.rmsx.nch')
    if not os.path.isfile(mapping_fname):
        fp = open(mapping_fname, 'w')
        fp.close()

def get_valid_chain_list(pdb_id, chains, residue_dict, modified_residue_dict, ref_seq_dict):
    valid_chains = []

    residue_abbreviation_dict = base_abbreviation('baselist.dat')     # baselist.bat downloaded from x3dna
    amino_acid_list = amino_acid_collection('aminoacidlist.dat')

    for chain_id in chains:
        if chain_id not in ref_seq_dict:
            print(pdb_id + '\tchain id: ' + chain_id + ' not found in ref_seq_dict. May be FASTA file is corrupted. Chain skipped.')
            write_invalid_mapping_file(pdb_id, chain_id)
            continue

        if chain_id not in residue_dict:
            print(pdb_id + '\tchain_id: ' + chain_id + ' residue_list None. Chain skipped.')
            write_invalid_mapping_file(pdb_id, chain_id)
            continue

        residue_dict[chain_id] = replace_modified_residue(pdb_id, chain_id, residue_dict[chain_id], modified_residue_dict, residue_abbreviation_dict, amino_acid_list)

        if len(residue_dict[chain_id]) == 0:
            print(pdb_id + '\tchain_id: ' + chain_id + ' residue_list Empty. Chain skipped.')
            write_invalid_mapping_file(pdb_id, chain_id)
            continue

        valid_chains.append(chain_id)       #chains that could be processed properly

    return valid_chains

def get_aln_mapping(aln_seq1, aln_seq2):
    """
    :return ret1: dict[i]->j  i in seq1; j in seq2
    :return ret2: dict[j]->i  j in seq2; i in seq1
    """
    if len(aln_seq1) != len(aln_seq2):
        return None

    i = j = 0

    ret1 = {}
    ret2 = {}
    for k in list(range(len(aln_seq1))):
        if aln_seq1[k] == "-" and aln_seq2[k] != "-":
            j += 1
        elif aln_seq2[k] == "-" and aln_seq1[k] != "-":
            i += 1
        elif aln_seq1[k] != "-" and aln_seq2[k] != "-":
            ret1[i] = j
            ret2[j] = i
            i += 1
            j += 1

    return ret1, ret2

def get_residue_reference_both_way_mapping_data_for_single_chain(residue_list, ref_seq):
    ind = len(residue_list) - 1
    while ind >= len(ref_seq) and residue_list[ind].symbol == 'X':
        ind -= 1

    residue_list = residue_list[:ind+1]

    ind = 0
    res_len = len(residue_list)
    while res_len - ind > len(ref_seq) and residue_list[ind].symbol == 'X':
        ind += 1

    residue_list = residue_list[ind:]
    
    residue_seq = ''.join(list(map(lambda x: x.symbol, residue_list)))    # get the residue sequence

    # aln = pairwise2.align.globalms(residue_seq, ref_seq, 5, -3, -10, -1)
    # aln = Align.PairwiseAligner.globalms(residue_seq, ref_seq, 5, -3, -10, -1)
    # (aln_residue, aln_ref, _, _, _) = aln[0]
    aln_residue, aln_ref = get_seq_alignment(residue_seq, ref_seq)

    ref_seq_replaced = replace_unknown_letter_in_ref(aln_ref, aln_residue)
    residue_to_ref_mapping, ref_to_residue_mapping = get_aln_mapping(aln_residue, ref_seq_replaced)

    return residue_list, ref_seq_replaced, residue_to_ref_mapping, ref_to_residue_mapping

def generate_pdbx_fasta_mapping_files_for_single_pdb(pdb_id, chains, pdbx_dir, fasta_dir, pdb_fasta_mapping_dir):
    pdb_fn = os.path.join(pdbx_dir, pdb_id + '.cif')
    fasta_fn = os.path.join(fasta_dir, pdb_id + '.fasta')

    # print(pdb_fn)
    # print(fasta_fn)

    if not os.path.isfile(pdb_fn) or not os.path.isfile(fasta_fn):
        print(pdb_id + '.cif' + ' or ' + pdb_id + '.fasta' + ' file is missing. Not generating the mapping file.')
        return

    residue_dict, missing_residue_dict, modified_residue_dict = load_pdb_data_from_file(pdb_fn)
    residue_dict = insert_missing_residue(residue_dict, missing_residue_dict)

    ref_seq_dict = load_fasta_seq(pdb_id, chains, fasta_dir)
        
    # for generating multi-chain sequence files
    multi_seq_dict = {}
    multi_chain_mapping = {}

    chains = get_valid_chain_list(pdb_id, chains, residue_dict, modified_residue_dict, ref_seq_dict)
    all_mapping_file_exists = True
    for chain_id in chains:
        chain_id = get_modified_chain_id_if_any_lowercase_letter(chain_id)
        mapping_fname = os.path.join(pdb_fasta_mapping_dir, pdb_id + '_' + chain_id + '.rmsx.nch')
        if not os.path.isfile(mapping_fname):
            all_mapping_file_exists = False
            break
    
    if all_mapping_file_exists == True:
        return

    if len(chains) > 1:
        multi_chain_mapping_fname = os.path.join(pdb_fasta_mapping_dir, pdb_id+'_'+'_'.join(chains)+'.rmsx.nch')
        fp_m = open(multi_chain_mapping_fname, 'w')

    for chain_id in chains:
        
        residue_dict[chain_id], ref_seq_dict[chain_id], residue_to_ref_mapping, ref_to_residue_mapping = get_residue_reference_both_way_mapping_data_for_single_chain(residue_dict[chain_id], ref_seq_dict[chain_id])
        
        multi_chain_mapping[chain_id] = {x.index: residue_to_ref_mapping[residue_dict[chain_id].index(x)] + sum(list(map(lambda x: len(multi_seq_dict[x]), multi_seq_dict))) for x in residue_dict[chain_id]}
        # multi_seq_dict[chain_id] = ref_seq_replaced
        multi_seq_dict[chain_id] = ref_seq_dict[chain_id]

        mapping_fname = os.path.join(pdb_fasta_mapping_dir, pdb_id + '_' + get_modified_chain_id_if_any_lowercase_letter(chain_id) + '.rmsx.nch')
        fp = open(mapping_fname, 'w')
        
        for m in residue_to_ref_mapping:
            fp.write(str(residue_dict[chain_id][m].index)+'\t'+str(residue_to_ref_mapping[m])+'\n')
            
            if len(chains) > 1:
                fp_m.write(str(residue_dict[chain_id][m].index)+'\t'+str(multi_chain_mapping[chain_id][residue_dict[chain_id][m].index])+'\n')

        fp.close()

    if len(chains) > 1:
        fp_m.close()

def _mapping_worker(p):
    generate_pdbx_fasta_mapping_files_for_single_pdb(*p)

def all_mapping_file_exists_for_single_pdb(pdb_id, chains, pdb_fasta_mapping_dir):
    valid_chains = []
    
    for chain_id in chains:
        chain_id = get_modified_chain_id_if_any_lowercase_letter(chain_id)
        mapping_fname = os.path.join(pdb_fasta_mapping_dir, pdb_id + '_' + chain_id + '_invalid.rmsx.nch')
        if not os.path.isfile(mapping_fname):
            valid_chains.append(chain_id)

    for chain_id in valid_chains:
        mapping_fname = os.path.join(pdb_fasta_mapping_dir, pdb_id + '_' + chain_id + '.rmsx.nch')
        if not os.path.isfile(mapping_fname):
            return False

    return True

def generate_pdbx_fasta_mapping_files(pdb_chains, pdb_dir, fasta_dir, pdb_fasta_mapping_dir):
    # if len(pdb_chains) == 0:
    #     print('Using existing PDBx-FASTA mapping files in ' + pdb_fasta_mapping_dir[base_path_len:] + '.')
    #     print('')
    #     return

    print('Generating PDBx-FASTA mapping files.')
    start_time = time.time()

    parameter_list = []
    for pdb_id in pdb_chains:
        # print(pdb_id, pdb_chains[pdb_id])
        if all_mapping_file_exists_for_single_pdb(pdb_id, pdb_chains[pdb_id], pdb_fasta_mapping_dir) == False:
            # generate_pdbx_fasta_mapping_files_for_single_pdb(pdb_id, pdb_chains[pdb_id], pdb_dir, fasta_dir, pdb_fasta_mapping_dir)
            parameter_list.append((pdb_id, pdb_chains[pdb_id], pdb_dir, fasta_dir, pdb_fasta_mapping_dir))

    pool = mp.Pool(8)
    pool.map(_mapping_worker, parameter_list)

    # wait_for_certain_time_according_to_wait_factor(len(parameter_list))
    print('Done')
    print('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

def get_pdbx_and_mapping_data(pdb_id, chains, pdbx_dir, fasta_dir):
    pdb_fn = os.path.join(pdbx_dir, pdb_id + '.cif')
    fasta_fn = os.path.join(fasta_dir, pdb_id + '.fasta')

    if not os.path.isfile(pdb_fn) or not os.path.isfile(fasta_fn):
        print(pdb_id + '.cif' + ' or ' + pdb_id + '.fasta' + ' file is missing.')
        return

    residue_dict, missing_residue_dict, modified_residue_dict = load_pdb_data_from_file(pdb_fn)
    residue_dict = insert_missing_residue(residue_dict, missing_residue_dict)

    ref_seq_dict = load_fasta_seq(pdb_id, chains, fasta_dir)

    chains = get_valid_chain_list(pdb_id, chains, residue_dict, modified_residue_dict, ref_seq_dict)

    res_to_ref = {}
    ref_to_res = {}
    for chain_id in chains:

        residue_dict[chain_id], ref_seq_dict[chain_id], residue_to_ref_mapping, ref_to_residue_mapping = get_residue_reference_both_way_mapping_data_for_single_chain(residue_dict[chain_id], ref_seq_dict[chain_id])

        res_to_ref[chain_id] = residue_to_ref_mapping
        ref_to_res[chain_id] = ref_to_residue_mapping

    return chains, residue_dict, ref_seq_dict, res_to_ref, ref_to_res, missing_residue_dict

def is_canonical_or_wobble(bp):
    if bp in ['AU', 'UA', 'GC', 'CG', 'GU', 'UG']:
        return True
    return False

def get_cWW_stack_length(wc_bp_dict, a, b):
    k = 1
    cWW_stack_len = 1
    while (a-k) in wc_bp_dict:
        found_stack = False
        for ind, _ in wc_bp_dict[a-k]:
            if ind - b == k:
                cWW_stack_len += 1
                found_stack = True
                break
        if found_stack == False:
            break
        k += 1

    k = 1
    while (a+k) in wc_bp_dict:
        found_stack = False
        for ind, _ in wc_bp_dict[a+k]:
            if ind - b == -1*k:
                cWW_stack_len += 1
                found_stack = True
                break
        if found_stack == False:
            break
        k += 1

    return cWW_stack_len

def remove_redundant_ww_bp(pdb_id, chain_id, wc_bp_dict, log_file_name):
    # stat generating code
    if "pool" not in log_file_name:
        for a in sorted(wc_bp_dict):
            if len(wc_bp_dict[a]) not in wc_count_stat:
                wc_count_stat[len(wc_bp_dict[a])] = 0
            wc_count_stat[len(wc_bp_dict[a])] += 1

    # cWW interaction conflict resolve
    
    for a in sorted(wc_bp_dict):
        filtered_bp_list = list(filter(lambda x: is_canonical_or_wobble(x[1]), wc_bp_dict[a]))
        if len(filtered_bp_list) > 0:
            for b, bp in wc_bp_dict[a]:
                if (b, bp) not in filtered_bp_list:
                    wc_bp_dict[b].remove((a, bp[::-1]))
                    # print "removing bp " + str(a) + ", " + str(b)
            wc_bp_dict[a] = filtered_bp_list

    new_wc_bp_dict = {}
    for a in sorted(wc_bp_dict):
        if len(wc_bp_dict[a]) == 1:
            new_wc_bp_dict[a] = wc_bp_dict[a][0][0]
            # if wc_bp_dict[a][0][0] == (1421, 'CG'):
            #     print('here1')
            #     sys.exit()
        elif len(wc_bp_dict[a]) > 1:
            selected_bp_ind = -1
            selected_bp = ""
            max_cWW_stack_len = 0
            for b, bp in sorted(wc_bp_dict[a]):
                cWW_stack_len = get_cWW_stack_length(wc_bp_dict, a, b)
                if cWW_stack_len > 1:
                    if cWW_stack_len > max_cWW_stack_len or (cWW_stack_len == max_cWW_stack_len and len(wc_bp_dict[b]) < len(wc_bp_dict[selected_bp_ind])):
                        if max_cWW_stack_len > 0:   # entering for the second time
                            fp = open(log_file_name, "a")
                            print(pdb_id + "\tchain id: " + chain_id + ", stack conflict (" + str(cWW_stack_len) + ", " + str(max_cWW_stack_len) + ") resolving for " + str(a) + ", " + str(b))
                            fp.write(pdb_id + "\tchain id: " + chain_id + ", stack conflict (" + str(cWW_stack_len) + ", " + str(max_cWW_stack_len) + ") resolving for " + str(a) + ", " + str(b) + "\n")
                            fp.close()
                        max_cWW_stack_len = cWW_stack_len
                        selected_bp_ind = b
                        selected_bp = bp

            if selected_bp_ind == -1:
                for b, bp in sorted(wc_bp_dict[a]):
                    if len(wc_bp_dict[b]) == 1:
                        selected_bp_ind = b

            if selected_bp_ind == -1:
                selected_bp_ind = sorted(wc_bp_dict[a])[0][0]

            for b, bp in wc_bp_dict[a]:
                if b != selected_bp_ind:
                    if (a, bp[::-1]) in wc_bp_dict[b]:
                    	wc_bp_dict[b].remove((a, bp[::-1]))
                    if b < a:
                        new_wc_bp_dict.pop(b)

                    # print "removing bp " + str(a) + ", " + str(b)

            new_wc_bp_dict[a] = selected_bp_ind
            # if selected_bp_ind == (1421, 'CG'):
            #     print('here2')
            #     sys.exit()
            wc_bp_dict[a] = [(selected_bp_ind, selected_bp)]


    return new_wc_bp_dict

def get_knot_free_struct(pdb_id, chain_id, ref_seq, residue_list, residue_to_ref_mapping, intera_list, log_file_name, is_detailed_ann):
    """
    get the pseudoknot-free secondary structure for the reference sequence
    ref_seq: the reference sequence
    residue_to_ref_mapping: the mapping from <residue sequence index> to <reference sequence index>
    residue_list: the list of residue objects
    intera_list: the base pair list (results of RNAview and MCAnnotate)
    :return: the wc base pairs in the knot free secondary structure
    """
    bp_item_len = 4
    if is_detailed_ann:
        bp_item_len = 5

    # get the watson-crick base pairs
    wc_bp_dict = {}
    residue_index_list = list(map(lambda x: x.index, residue_list))

    for intera in intera_list:
        # print(chain_id)
        # print(intera[0].chain_id)
        # print(intera[1].chain_id)
        # print(len(intera))
        # print(intera[2])
        # print(intera[3])
        # sys.exit()
        if (intera[0].chain_id == chain_id and intera[1].chain_id == chain_id) and (intera[0] in residue_index_list and intera[1] in residue_index_list) and len(intera) == bp_item_len and intera[2] == 'W/W' and intera[3] == 'cis':
            # print('here')
            # sys.exit()
            # index in residue sequence
            i = residue_index_list.index(intera[0])
            j = residue_index_list.index(intera[1])

            # index at the reference sequence
            mi = residue_to_ref_mapping[i]
            mj = residue_to_ref_mapping[j]

            # get the letters in reference sequence
            s1 = ref_seq[mi]
            s2 = ref_seq[mj]

            # find the watson-crick base pairs
            #if s1+s2 in ['AU', 'UA', 'GC', 'CG', 'GU', 'UG']:
            if mi not in wc_bp_dict:
                wc_bp_dict[mi] = []
            if mj not in wc_bp_dict:
                wc_bp_dict[mj] = []

            wc_bp_dict[mi].append((mj, s1+s2))
            wc_bp_dict[mj].append((mi, s2+s1))

            # if len(intera) == bp_item_len and intera[2] == 'W/W' and intera[3] == 'cis':
            #     wc_bp_dict[mi] = mj
            #     wc_bp_dict[mj] = mi

    # remove the inconsistent index from wc base pairs (i--->j; j-X->i)
    # map(lambda y: wc_bp_dict.pop(y), filter(lambda x: x != wc_bp_dict[wc_bp_dict[x]], wc_bp_dict.keys()))
    # print('before')
    # print(wc_bp_dict)
    wc_bp_dict = remove_redundant_ww_bp(pdb_id, chain_id, wc_bp_dict, log_file_name)
    # print('after')
    # print(wc_bp_dict)

    # if wc base pairs
    if len(wc_bp_dict) != 0:
        input_fp, input_path = tempfile.mkstemp()
        output_fp, output_path = tempfile.mkstemp()

        os.write(input_fp, str.encode(str(len(ref_seq))+"\n"))
        for i in list(range(len(ref_seq))):
            if i in wc_bp_dict:
                # print(wc_bp_dict[i])
                # print(str(i+1)+'\t'+ref_seq[i]+'\t'+str(i)+'\t'+str(i+2)+'\t'+str(wc_bp_dict[i]+1) + '\t'+str(i+1)+'\n')
                os.write(input_fp, str.encode(str(i+1)+'\t'+ref_seq[i]+'\t'+str(i)+'\t'+str(i+2)+'\t'+str(wc_bp_dict[i]+1) + '\t'+str(i+1)+'\n'))

            else:
                os.write(input_fp, str.encode(str(i+1)+'\t'+ref_seq[i]+'\t'+str(i)+'\t'+str(i+2)+'\t'+'0'+'\t'+str(i+1)+'\n'))
        os.close(input_fp)
        os.close(output_fp)
        
        # remove the pseudoknot by using RNAstructure
        # the path should be changed if deployed
        # os.environ["DATAPATH"] = "/home/pge/Workspace/RMDB/tools/data_tables"
        # os.system("/home/pge/Workspace/RMDB/tools/RemovePseudoknots -m %s %s 1>/dev/null" % (input_path, output_path))

        # if flag_remove_pseudoknots is True:
        #     os.system("../lib/k2n_standalone/knotted2nested.py %s > %s" % (input_path, output_path))
        os.system("k2n_standalone/knotted2nested.py %s > %s" % (input_path, output_path))

        # ret = []
        # lineno = 0
        # if flag_remove_pseudoknots is True:
        #     fp = open(output_path)
        # else:
        #     fp = open(input_path)

        # fp = open(output_path)
        # for line in fp.readlines():
        #     if lineno == 0:
        #         lineno += 1
        #         continue
        #     line = line.strip().split()
        #     if line[4] != '0' and int(line[0]) < int(line[4]):
        #         ret.append((int(line[0])-1, int(line[4])-1))
        # fp.close()

        # if flag_remove_pseudoknots is False:
        ret, pseudoknots = extract_pseudoknot(input_path, output_path)

        # print pseudoknots
        # sys.exit()
        # print input_path
        # print output_path
        os.remove(input_path)
        os.remove(output_path)

        return ret, pseudoknots
    else:
        return [], []

def extract_pseudoknot(input_path, output_path):

    with_knot = []
    lineno = 0
    fp = open(input_path)
    for line in fp.readlines():
        if lineno == 0:
            lineno += 1
            continue
        line = line.strip().split()
        if line[4] != '0' and int(line[0]) < int(line[4]):
            with_knot.append((int(line[0])-1, int(line[4])-1))
    fp.close()

    without_knot = []
    lineno = 0
    fp = open(output_path)
    for line in fp.readlines():
        if lineno == 0:
            lineno += 1
            continue
        line = line.strip().split()
        if line[4] != '0' and int(line[0]) < int(line[4]):
            without_knot.append((int(line[0])-1, int(line[4])-1))
    fp.close()

    knots = []
    for item in with_knot:
        if item not in without_knot:
            knots.append(item)

    return without_knot, knots

def get_loop_in_ss(ss_bp, seq_len, stack_size_cutoff):
    lifo_stack = []
    root = StackTNode(-1, seq_len, 0, [])
    lifo_stack.append(root)

    for i, j, bp in ss_bp:
        is_cano = False
        if bp in ["AU", "UA", "GC", "CG"]:
            is_cano = True
        # continous base pairs
        if lifo_stack[-1].i+lifo_stack[-1].size == i and lifo_stack[-1].j-lifo_stack[-1].size == j:
            if is_cano == True:
                lifo_stack[-1].set_cano_status()
            lifo_stack[-1].size += 1
        # juxtaposing or enclosing
        else:
            # juxtaposing
            if lifo_stack[-1].j < i:
                while True:
                    lifo_stack.pop()
                    if lifo_stack[-1].j > i:
                        break

            node = StackTNode(i, j, 1, [], is_cano)
            lifo_stack[-1].children.append(node)
            lifo_stack.append(node)

    root.filter_cWW_stack(stack_size_cutoff)
    # root.display()
    ret = []
    root.loop(ret)
    return ret

def get_loop_end(loop):
    ret = []
    for i, r in enumerate(loop):
        ret.append((r[0], loop[i-1][1]))
    return ret

def rotate(l, x):
    return l[-x:] + l[:-x]

def get_loop_intera(loop, loop_end, residue_list, ref_to_residue_mapping, residue_to_ref_mapping, intera_list, is_detailed_ann):
    intera_in_loop = []
    residue_index_in_loop = []
    residue_index_to_loop_index ={}
    loop_end_residue = []
    loop_i = 0
    bp_item_len = 4
    if is_detailed_ann:
        bp_item_len = 5

    for i, j in loop_end:
        loop_end_residue.append((residue_list[ref_to_residue_mapping[i]].index, residue_list[ref_to_residue_mapping[j]].index))

    for r in loop:
        for v in list(range(r[0], r[1]+1)):
            # get the mapped reference index
            if v in ref_to_residue_mapping:
                residue_index_in_loop.append(residue_list[ref_to_residue_mapping[v]].index)
                residue_index_to_loop_index[residue_list[ref_to_residue_mapping[v]].index] = loop_i
            loop_i += 1

    for intera in intera_list:
        if ((intera[0], intera[1]) in loop_end_residue or (intera[1], intera[0]) in loop_end_residue) and \
                (len(intera) == bp_item_len and intera[2] == 'W/W' and intera[3] == 'cis'):
            continue

        if intera[0] in residue_index_in_loop and intera[1] in residue_index_in_loop:
            i1 = residue_index_to_loop_index[intera[0]]     # index in loop
            i2 = residue_index_to_loop_index[intera[1]]     # index in loop
            intera_in_loop.append((i1, i2, intera))

    return intera_in_loop

def generate_missing_residue_statistics(loop_i_index, missing_residue_list, residue_list, ref_to_residue_mapping):
    # loop_i_residue = []
    segment_lengths = []
    count = 0
    flag = False
    for item in loop_i_index:
        if item in ref_to_residue_mapping:# and ref_to_residue_mapping[item] in residue_list:
            residue = residue_list[ref_to_residue_mapping[item]]
            # loop_i_residue.append(residue)
            if residue in missing_residue_list:
                count += 1
                flag = True
            else:
                if flag == True:
                    segment_lengths.append(count)
                flag = False
                count = 0
    if count > 0:
        segment_lengths.append(count)

    # if len(segment_lengths) > 0:
    #     print segment_lengths
    #     sys.exit()
    return segment_lengths

def remove_unknown_bp(lines):
    ret = []
    dot_count = 0
    flag_read = False
    for line in lines:
        if line.startswith("#info=basepair"):
            flag_read = True
        elif line.startswith("#info=stacking"):
            flag_read = False
        elif flag_read is True:
            edges = line.strip().split(",")[1]
            parts = edges.strip().split("/")
            if parts[0] == "." or parts[1] == ".":
                dot_count += 1
                continue
        ret.append(line)
    # print str(dot_count) + " base pairs removed."
    return ret

def get_mean(a_list):
    if is_median:
        return round(numpy.median(a_list), 1)
    return round(numpy.mean(a_list), 1)

def filter_outlier_loops(junction_wise_loops, outlier_zscore_threshold):
    junction_wise_filtered_loops = {}
    junction_wise_outlier_loops = {}

    for junc_cnt in junction_wise_loops:
        junction_wise_loops[junc_cnt] = sorted(junction_wise_loops[junc_cnt])

    for junc_cnt in junction_wise_loops:
        total_item = len(junction_wise_loops[junc_cnt])
        truncate_count = int(round(total_item * truncate_percentage / 100))
        start = truncate_count
        end = total_item - truncate_count

        mean = get_mean(list(map(lambda x: x[0], junction_wise_loops[junc_cnt][start:end])))
        std = numpy.std(list(map(lambda x: x[0], junction_wise_loops[junc_cnt][start:end])))

        length_threshold = int(round(mean + std * outlier_zscore_threshold))
        # length_threshold_n = int(round(mean - std * outlier_zscore_threshold))
        # print junc_cnt, length_threshold, mean, std

        for loop_length, loop in junction_wise_loops[junc_cnt]:
            # z_score = zscore(float(loop_length), float(mean), float(std))
            # if z_score > outlier_zscore_threshold:
            if loop_length > length_threshold:
                if junc_cnt not in junction_wise_outlier_loops:
                    junction_wise_outlier_loops[junc_cnt] = []
                junction_wise_outlier_loops[junc_cnt].append(loop)

            # elif loop_length < length_threshold_n:
            #     if junc_cnt not in junction_wise_outlier_loops:
            #         junction_wise_outlier_loops[junc_cnt] = []
            #     junction_wise_outlier_loops[junc_cnt].append(loop)

            else:
                if junc_cnt not in junction_wise_filtered_loops:
                    junction_wise_filtered_loops[junc_cnt] = []
                junction_wise_filtered_loops[junc_cnt].append(loop)

    return junction_wise_filtered_loops, junction_wise_outlier_loops

def filter_empty_loops(loop_dir, loop_subdirs, flag_remove_empty_loops):
    filtered_loops = {}
    empty_loops = {}
    for subdir in loop_subdirs:
        # filtered_loops[subdir] = []
        # empty_loops[subdir] = []

        for fn in glob.glob(os.path.join(loop_dir + "/" + subdir, "*.smf")):
            
            loop = os.path.basename(fn)[:-4]
            loop_segments = loop.strip().split(":")[1].strip().split("_")
            junction_count = len(loop_segments)
            loop_length = get_loop_length(loop_segments)

            fp = open(fn)
            lines = fp.readlines()
            lines = remove_unknown_bp(lines)
            fp.close()

            if flag_remove_empty_loops:
                if lines[3].startswith("#info=stacking"):
                    if junction_count not in empty_loops:
                        empty_loops[junction_count] = []
                    empty_loops[junction_count].append(loop)
                    continue

            #line_count = sum(1 for line in open(fn)) # get the count of lines in the file
            # if line_count == 4:     #3 for pmf, 4 for rmf and smf
            #     continue

            if junction_count not in filtered_loops:
                filtered_loops[junction_count] = []
            filtered_loops[junction_count].append((loop_length, loop))

    return filtered_loops, empty_loops

def loopcut(pdbx_dir, fasta_dir, pdb_fasta_mapping_dir, annotation_source, annotation_dir, loop_dir, pdb_id, chains, log_file_name):
    is_detailed_ann = True
    bp_item_len = 4
    stk_item_len = 3
    if is_detailed_ann:
        bp_item_len = 5
        stk_item_len = 4
    # print "processing %s with chains %s ..." % (pdb_id + ".cif", ",".join(chains))
    print('Processing ' + pdb_id + ' with chains ' + ','.join(chains))

    annotation_list = parseMergedAnnotation(os.path.join(annotation_dir, annotation_source.lower(), pdb_id + '.' + annotation_source.lower()), is_detailed_ann)

    if len(annotation_list) == 0:
        print(pdb_id + "\tNo annotation found!")
        # return
    
    if os.path.isfile(os.path.join(annotation_dir, "knot_free/" + pdb_id + ".kf")):
        annotation_list = parseMergedAnnotation(os.path.join(annotation_dir, "knot_free/" + pdb_id + ".kf"), is_detailed_ann)
    else:
        print("WARNING! Knot-free annotation for " + pdb_id + " NOT FOUND.")

    # pdb_fn = pdbx_dir + pdb_id + ".cif"
    pdb_fn = os.path.join(pdbx_dir, pdb_id + ".cif")

    # residue_list_dict = {}    #it was in previous code but I (MMR) found it is not used further.

    # for generating multi-chain sequence files
    multi_seq_dict = {}
    multi_chain_mapping = {}

    if len(chains) > 1:
        fp_m = open(os.path.join(pdb_fasta_mapping_dir, pdb_id+"_"+"_".join(chains)+".rmsx.nch"), "a")
    # generate loop files for multiple chain

    valid_chains = []
    loop_with_missing_res_dict = {}
    knot_free_cWW_annotations = {}
    # print(pdb_id, chains, pdbx_dir)
    chains, residue_dict, ref_seq_dict, res_to_ref, ref_to_res, missing_residue_dict = get_pdbx_and_mapping_data(pdb_id, chains, pdbx_dir, fasta_dir)
    residue_with_missing_base_dict = missing_residue_dict

    for chain_id in chains:

        residue_list = residue_dict[chain_id]
        ref_to_residue_mapping = ref_to_res[chain_id]
        residue_to_ref_mapping = res_to_ref[chain_id]
        ref_seq_replaced = ref_seq_dict[chain_id]
        # print(residue_list, ref_to_residue_mapping, residue_to_ref_mapping, ref_seq_replaced)
        # print(annotation_list)
        # sys.exit()
        ss_bp, knots = get_knot_free_struct(pdb_id, chain_id, ref_seq_replaced, residue_list, residue_to_ref_mapping, annotation_list, log_file_name, True)
        # print(ss_bp, knots)
        knot_free_cWW_annotations[chain_id] = []
        ss_bp_with_symbol = []
        for s, e in ss_bp:
            res1 = residue_list[ref_to_residue_mapping[s]]
            res2 = residue_list[ref_to_residue_mapping[e]]
            ss_bp_with_symbol.append((s, e, res1.symbol+res2.symbol))
            knot_free_cWW_annotations[chain_id].append((res1.index, res2.index))

        # get all the loops in the pseudoknot free structure
        loops = get_loop_in_ss(ss_bp_with_symbol, len(ref_seq_replaced), cWW_stack_threshold)
        # print(loops)
        # sys.exit()
        # loop_with_missing_res_dict[chain_id] = 0
        for loop, type in loops:
            loop_end = get_loop_end(loop)
            for i in list(range(len(loop))):
                loop_i = rotate(loop, i)
                loop_i_index = reduce(lambda y, z: y+z, list(map(lambda x: list(range(x[0], x[1]+1)), loop_i)))
                
                # loop_intera = util.get_loop_intera(loop_i, residue_list, ref_to_residue_mapping, residue_to_ref_mapping, annotation_list)
                loop_intera = get_loop_intera(loop_i, loop_end, residue_list, ref_to_residue_mapping, residue_to_ref_mapping, annotation_list, is_detailed_ann)
                loop_bp_intera = list(filter(lambda x: len(x[2]) == bp_item_len, loop_intera))
                loop_stack_intera = list(filter(lambda x: len(x[2]) == stk_item_len, loop_intera))   #activated to get stackings
                loop_bp_intera = sorted(loop_bp_intera, key=lambda x: x[0])
                loop_stack_intera = sorted(loop_stack_intera, key=lambda x: x[0])   #activated to get stackings

                # filter the loop without base pairing iteractions
                # if len(filter(lambda x: len(x[2]) == 4, loop_intera)) == 0:
                #    continue
                create_directory(os.path.join(loop_dir, type))
                loop_fn = os.path.join(loop_dir, type, pdb_id+"_%s:%s.pmf" % (chain_id, "_".join(list(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i)))))  #for scan
                loop_fn_rmf = os.path.join(loop_dir, type, pdb_id+"_%s:%s.rmf" % (chain_id, "_".join(list(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i)))))  #for scanx without stacking
                loop_fn_smf = os.path.join(loop_dir, type, pdb_id+"_%s:%s.smf" % (chain_id, "_".join(list(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i)))))  #for scanx with stacking

                fms = open(log_file_name + "_missing_residue_stat.txt", "a")
                if chain_id in missing_residue_dict or chain_id in residue_with_missing_base_dict:
                    missing_res_segment_lengths = []
                    if chain_id in residue_with_missing_base_dict:
                        missing_res_segment_lengths += generate_missing_residue_statistics(loop_i_index, residue_with_missing_base_dict[chain_id], residue_list, ref_to_residue_mapping)
                        # if len(missing_res_segment_lengths) > 0:
                        #     print missing_res_segment_lengths
                    if chain_id in missing_residue_dict:
                        missing_res_segment_lengths += generate_missing_residue_statistics(loop_i_index, missing_residue_dict[chain_id], residue_list, ref_to_residue_mapping)
                    if len(missing_res_segment_lengths) > 0:
                        loop_with_missing_res_dict[pdb_id, chain_id, type, str(loop_i)] = (len(loop_i_index), missing_res_segment_lengths)
                        fms.write(pdb_id + "\t" + chain_id + "\t" + type + "\t" + str(loop) + "\t" + str(len(loop_i_index)))
                        # print pdb_id + "\t" + chain_id + "\t" + type + "\t" + str(loop) + "\t" + str(len(loop_i_index))
                        missing_res_length = 0
                        for item in missing_res_segment_lengths:
                            # print "\t" + item
                            fms.write("\t" + str(item))
                            missing_res_length += item
                        fms.write("\n")

                        if missing_res_length > 5 or 2 * missing_res_length > len(loop_i_index):
                            create_directory(os.path.join(loop_dir, 'missing_res', type))
                            loop_fn = os.path.join(loop_dir, "missing_res/" + type, pdb_id+"_%s:%s.pmf" % (chain_id, "_".join(list(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i)))))  #for scan
                            loop_fn_rmf = os.path.join(loop_dir, "missing_res/" + type, pdb_id+"_%s:%s.rmf" % (chain_id, "_".join(list(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i)))))  #for scanx without stacking
                            loop_fn_smf = os.path.join(loop_dir, "missing_res/" + type, pdb_id+"_%s:%s.smf" % (chain_id, "_".join(list(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i)))))  #for scanx with stacking
                fms.close()

                # fp = open(loop_fn, 'w')
                # fr = open(loop_fn_rmf, 'w')
                fs = open(loop_fn_smf, 'w')

                # fp.write(">%s_%s:%s\n" % (pdb_id, chain_id, "_".join(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i))))
                # fp.write("...".join(map(lambda x: ref_seq_replaced[x[0]:x[1]+1], loop_i))+"\n")

                # fr.write(">%s_%s:%s\n" % (pdb_id, chain_id, "_".join(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i))))
                # fr.write("...".join(map(lambda x: ref_seq_replaced[x[0]:x[1]+1], loop_i))+"\n")

                fs.write(">%s_%s:%s\n" % (pdb_id, chain_id, "_".join(list(map(lambda x: str(x[0])+"-"+str(x[1]), loop_i)))))
                fs.write("...".join(list(map(lambda x: ref_seq_replaced[x[0]:x[1]+1], loop_i)))+"\n")

                # fp.write("#info=basepair\n")
                # fr.write("#info=basepair\n")
                fs.write("#info=basepair\n")

                # fp.write("//\n")
                # fr.write("//\n")
                # fs.write("//\n")

                for i, j, intera in loop_bp_intera:
                    if intera[3] == "hbond":
                        continue
                    if i < j:
                        # fp.write("#\t%d-%d : %s-%s,%s %s\n" % (i, j, intera[2], intera[3], intera[0], intera[1]))
                        # fp.write("#\t%d-%d : %s-%s %s %s\n" % (i, j, ref_seq_replaced[loop_i_index[i]], ref_seq_replaced[loop_i_index[j]], intera[2], intera[3]))

                        # fr.write("%d-%d,%s,%s\n" % (i, j, intera[2], intera[3]))

                        fs.write("%d-%d,%s,%s\n" % (i, j, intera[2], intera[3]))
                    else:
                        # fp.write("#\t%d-%d : %s-%s %s %s\n" % (j, i, intera[2][::-1], intera[3], intera[1], intera[0]))
                        # fp.write("#\t%d-%d : %s-%s %s %s\n" % (j, i, ref_seq_replaced[loop_i_index[j]], ref_seq_replaced[loop_i_index[i]], intera[2][::-1], intera[3]))

                        # fr.write("%d-%d,%s,%s\n" % (j, i, intera[2][::-1], intera[3]))

                        fs.write("%d-%d,%s,%s\n" % (j, i, intera[2][::-1], intera[3]))

                fs.write("#info=stacking\n")
                for i, j, intera in loop_stack_intera:
                    if i < j:
                        fs.write("%d-%d,%s\n" % (i, j, intera[2]))
                    else:
                        if intera[2] == "upward":
                            fs.write("%d-%d,%s\n" % (j, i, "downward"))
                        elif intera[2] == "downward":
                            fs.write("%d-%d,%s\n" % (j, i, "upward"))
                        else:
                            fs.write("%d-%d,%s\n" % (j, i, intera[2]))

                # fp.close()
                # fr.close()
                fs.close()

        valid_chains.append(chain_id)       #chains that could be processed properly

    # write_knot_free_annotations(pdb_id, annotation_list, knot_free_cWW_annotations, annotation_dir, is_detailed_ann)

    # print chains
    chains = valid_chains
    if len(chains) == 0:    #chains may become empty after removing a chain that could not be processed
        return

    if len(chains) > 1:
        fp_m.close()

    if include_multi_chain_loops == True:
        cross_bp_intera = list(filter(lambda x: (x[0].chain_id != x[1].chain_id) and (x[0].chain_id in chains)
                                            and (x[1].chain_id in chains) and (len(x) == bp_item_len), annotation_list))
        # cross_stack_intera = filter(lambda x: (x[0].chain_id != x[1].chain_id) and (x[0].chain_id in chains)
        #                                    and (x[1].chain_id in chains) and (len(x) == stk_item_len), annotation_list)

        # multi_chain_seq_ref = "".join(map(lambda x: multi_seq_dict[x], chains))
        multi_chain_fn = os.path.join(loop_dir, "CC", pdb_id+"_%s.in" % "_".join(chains))
        # print "multi: " + str(multi_seq_dict)
        # sys.exit()
        fp = open(multi_chain_fn, "w")
        multi_chain_seq_sep = "...".join(list(map(lambda x: multi_seq_dict[x], chains)))

        fp.write(">%s_%s\n" % (pdb_id, "_".join(chains)))
        fp.write(multi_chain_seq_sep+'\n')
        fp.write("#info=basepair\n")
        for intera in cross_bp_intera:
            if intera[3] == "hbond":
                continue
            if (intera[0] in multi_chain_mapping[intera[0].chain_id]) and (intera[1] in multi_chain_mapping[intera[1].chain_id]):
                i = multi_chain_mapping[intera[0].chain_id][intera[0]]
                j = multi_chain_mapping[intera[1].chain_id][intera[1]]
                if i < j:
                    fp.write("%d-%d,%s,%s,%s-%s\n" % (i, j, intera[2], intera[3], intera[0], intera[1]))
                else:
                    fp.write("%d-%d,%s,%s,%s-%s\n" % (j, i, intera[2][::-1], intera[3], intera[1], intera[0]))
        fp.write("#info=stacking\n")
        fp.close()

        # make sure the path of RNAMotifScanX is current
        #os.system("/home/pge/Workspace/MotifCluster/tools/RNAMotifScanX-release/bin/cut %s --out_dir %s" % (multi_chain_fn, os.path.join(loop_dir, "CC")))

        # print multi_chain_fn
        # sys.exit()
        # os.system("../lib/RNAMotifScanX-release/bin/cut %s --out_dir %s" % (multi_chain_fn, os.path.join(loop_dir, "CC")))

        for fn in glob.glob(os.path.join(loop_dir, "CC", "%s*.rmf" % pdb_id)):
            fn_base = os.path.basename(fn)[:-4]
            chain_id, regions = fn_base.split(':')
            output_fn = os.path.join(loop_dir, "CC", chain_id+":"+regions+".pmf")
            fp = open(output_fn, 'w')
            fp.write(">%s:%s\n" % (chain_id, regions))

            multi_chain_loop_sep = []
            for r in regions.split('_'):
                i, j = r.split('-')
                multi_chain_loop_sep.append(multi_chain_fragment(multi_chain_seq_sep, (int(i), int(j))))
            fp.write("...".join(multi_chain_loop_sep)+'\n')

            fp.write("//\n")

            line_no = 0
            fp_i = open(fn)
            for line in fp_i.readlines():
                if line.startswith("#info=stacking"):
                    break
                elif line_no > 2:
                    decom = line.strip().split(',')
                    i, j = decom[0].split('-')
                    c1 = multi_chain_loop_ref[int(i)]
                    c2 = multi_chain_loop_ref[int(j)]
                    fp.write("#\t%s-%s : %s-%s %s %s\n" % (i, j, c1, c2, decom[1], decom[2]))
                elif line_no == 1:
                    line = line.strip()
                    multi_chain_loop_ref = line.replace(".", "")
                line_no += 1
            fp_i.close()
            fp.close()

def _loopcut_worker(p):
    loopcut(*p)

def cut_loops(pdb_chain_dict, directories, ann_source, remove_empty_loops, outlier_zscore_threshold):
    if ann_source == None:
        ann_source = 'merged'

    pdb_dir = directories.pdbx_dir
    fasta_dir = directories.fasta_dir
    pdb_fasta_mapping_dir = directories.pdb_fasta_mapping_dir
    ann_dir = directories.annotation_dir
    # loop_dir = directories.loop_dir

    # create_directory(pdb_fasta_mapping_dir)
    # generate_pdbx_fasta_mapping_files(pdb_chain_dict, pdb_dir, fasta_dir, pdb_fasta_mapping_dir)
    # create_directory(loop_dir)
    temp_loop_dir = os.path.join(directories.base_dir, 'temp_loops')
    if os.path.exists(temp_loop_dir):
        remove_all_from_dir(temp_loop_dir)
        delete_directory(temp_loop_dir)
    create_directory(temp_loop_dir)

    log_file_name = "main_log.log"
    fpt = open(log_file_name, "w")
    fpt.write("\n")
    fpt.close()

    parameter_list = []
    for pdb_id in pdb_chain_dict:
        chains = pdb_chain_dict[pdb_id]
        loopcut(pdb_dir, fasta_dir, pdb_fasta_mapping_dir, ann_source, ann_dir, temp_loop_dir, pdb_id, chains, log_file_name)
        # parameter_list.append((pdb_dir, fasta_dir, pdb_fasta_mapping_dir, ann_source, ann_dir, temp_loop_dir, pdb_id, chains, log_file_name))
    
    # print parameter_list
    pool = mp.Pool(8)
    pool.map(_loopcut_worker, parameter_list)

    print("loopcut done.")

    # print(loop_dir)
    # sys.exit()
    loop_subdirs = ["HL", "IL", "ML"]#, "CC"]
    # loop_subdirs = ["HL", "IL", "ML", "CC", "HL_nc-helix", "IL_nc-helix", "ML_nc-helix"]
    # loop_subdirs = ["HL", "IL", "ML"]
    for subdir in loop_subdirs:
        if not os.path.exists(os.path.join(temp_loop_dir, subdir)):
            os.makedirs(os.path.join(temp_loop_dir, subdir))
        if not os.path.exists(os.path.join(temp_loop_dir, "missing_res/" + subdir)):
            os.makedirs(os.path.join(temp_loop_dir, "missing_res/" + subdir))
        # if not os.path.exists(os.path.join(loop_dir, "outlier_loops/" + subdir)):
        #     os.makedirs(os.path.join(loop_dir, "outlier_loops/" + subdir))
        # if not os.path.exists(os.path.join(loop_dir, "empty_loops/" + subdir)):
        #     os.makedirs(os.path.join(loop_dir, "empty_loops/" + subdir))

    filtered_loops, empty_loops = filter_empty_loops(temp_loop_dir, loop_subdirs, remove_empty_loops)
    junction_wise_filtered_loops, junction_wise_outlier_loops = filter_outlier_loops(filtered_loops, outlier_zscore_threshold)

    remove_all_from_dir(temp_loop_dir)
    delete_directory(temp_loop_dir)

    return junction_wise_filtered_loops

    # for junc_cnt in junction_wise_filtered_loops:

    #     loop_type_index = junc_cnt-1
    #     if loop_type_index > 2:
    #         loop_type_index = 2

    #     input_dir = os.path.join(temp_loop_dir, loop_subdirs[loop_type_index])
    #     output_dir = os.path.join(temp_loop_dir, loop_subdirs[loop_type_index] + "D")
    #     outlier_dir = os.path.join(temp_loop_dir, "outlier_loops/" + loop_subdirs[loop_type_index])
    #     empty_loop_dir = os.path.join(temp_loop_dir, "empty_loops/" + loop_subdirs[loop_type_index])

    #     if not os.path.exists(output_dir):
    #         os.makedirs(output_dir)
    #     if not os.path.exists(outlier_dir):
    #         os.makedirs(outlier_dir)
    #     if not os.path.exists(empty_loop_dir):
    #         os.makedirs(empty_loop_dir)

    #     if junc_cnt > 2:
    #         output_dir2 = os.path.join(output_dir, loop_subdirs[loop_type_index] + "D_" + str(junc_cnt) + "way")
    #         if not os.path.exists(output_dir2):
    #             os.makedirs(output_dir2)
    #         write_loop_to_directory(input_dir, junction_wise_filtered_loops[junc_cnt], output_dir2, extension)

    #     if not os.path.exists(output_dir):
    #         os.makedirs(output_dir)

    #     write_loop_to_directory(input_dir, junction_wise_filtered_loops[junc_cnt], output_dir, extension)

    #     if junc_cnt in junction_wise_outlier_loops:
    #         write_loop_to_directory(input_dir, junction_wise_outlier_loops[junc_cnt], outlier_dir, extension)

    #     if junc_cnt in empty_loops:
    #         write_loop_to_directory(input_dir, empty_loops[junc_cnt], empty_loop_dir, extension)

    # # filter_loops(loop_dir, outlier_dir, loop_subdirs, "D", flag_remove_empty_loops, junction_wise_outlier_loop_threshold)
    # print("loop filtering done.")