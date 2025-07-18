import sys
import os
import logging
import time
# from Bio import SeqIO, pairwise2
from Bio import SeqIO, Align
import multiprocessing as mp

# python 3 compatibility
from functools import reduce
from past.builtins import map
from future.standard_library import install_aliases
install_aliases()
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request, urlretrieve
from urllib.error import HTTPError

# sys.path.append('../../')
# from config import *
# sys.path.append(scripts_dir)
from logger import *
from config import *
from utils import *
from cif import *
from ann_generator import *
from ann_parser import *
# from simplified_loopcut_helper import *

def load_fasta_seq(pdb_id, chains, directories):
    fasta_seq_dict = {}
    fasta_fn = os.path.join(directories.fasta_dir, pdb_id + '.fasta')
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

def get_pdbx_and_fasta_files(pdb_list, pdbx_url, fasta_url, directories, content_download_attempts):
    success = True

    print('')
    if len(pdb_list) == 0:
        base_path_len = len(os.path.dirname(directories.base_dir)) + 1
        logger.info('Using existing PDBx and FASTA files in ' + directories.pdbx_dir[base_path_len:] + ' and ' + directories.fasta_dir[base_path_len:] + ' respectively.')
        print('')
        return

    logger.info('Getting PDBx and FASTA files.')
    start_time = time.time()
    # parameter_list = []
    for pdb in pdb_list:
        pdb_fname = pdb + '.cif'
        fasta_fname = pdb + '.fasta'
        if not os.path.isfile(os.path.join(directories.pdbx_dir, pdb_fname)):
            # parameter_list.append((pdb_fname, 'cif'))
            success &= download_single_pdbx_or_fasta_file(pdb_fname, 'cif', pdbx_url, directories, content_download_attempts)
        if not os.path.isfile(os.path.join(directories.fasta_dir, fasta_fname)):
            # parameter_list.append((fasta_fname, 'fasta'))
            success &= download_single_pdbx_or_fasta_file(fasta_fname, 'fasta', fasta_url, directories, content_download_attempts)

    # pool = mp.Pool(number_of_multiprocess)
    # pool.map(_pdb_and_fasta_download_worker, parameter_list)

    if success == True:
        logger.info('Done')
        logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')
    else:
        logger.error('Error downloading PDBx or FASTA files, process terminating.')
        sys.exit()

# def _pdb_and_fasta_download_worker(p):
    # download_single_pdbx_or_fasta_file(*p)

def download_single_pdbx_or_fasta_file(fname, file_ext, url, directories, download_attempts):
    status = False
    number_of_tries = 1
    wait_time = 1
    print(fname)
    while status == False and number_of_tries <= download_attempts:
        if number_of_tries > 1:
            logger.info('Repeating download: ' + fname + ' (Attempt ' + str(number_of_tries) + ').')

        if file_ext == 'cif':
            urlretrieve(url + fname, os.path.join(directories.pdbx_dir, fname))
            wait_for_certain_time_according_to_wait_factor(wait_time)
            fp = open(os.path.join(directories.pdbx_dir, fname))
            lines = fp.readlines()
            fp.close()
            if lines[0].startswith('data_' + fname.strip().split('.')[0]) and has_complete_ATOM_portion(lines) and lines[-1].strip() == '#':
                status = True

        elif file_ext == 'fasta':
            urlretrieve(url + fname.strip().split('.')[0], os.path.join(directories.fasta_dir, fname))
            wait_for_certain_time_according_to_wait_factor(wait_time)
            fp = open(os.path.join(directories.fasta_dir, fname))
            lines = fp.readlines()
            fp.close()
            if lines[0].startswith('>' + fname.strip().split('.')[0]):
                status = True
        else:
            logger.error('Invalid file_type provided to download.')
            break

        wait_time += number_of_tries
        number_of_tries += 1

    if status == False:
        logger.error(fname + ' download unsuccessful.')
        if file_ext == 'cif':
            if os.path.isfile(os.path.join(directories.pdbx_dir, fname)):
                os.remove(os.path.join(directories.pdbx_dir, fname))
        elif file_ext == 'fasta':
            if os.path.isfile(os.path.join(directories.fasta_dir, fname)):
                os.remove(os.path.join(directories.fasta_dir, fname))

        return False
    return True

def has_complete_ATOM_portion(lines):
    is_complete = False
    in_section = False
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            in_section = True
        elif line.startswith('#') and in_section == True:
            is_complete = True
            break

    return is_complete

def all_mapping_file_exists_for_single_pdb(pdb_id, chains, directories):
    valid_chains = []
    
    for chain_id in chains:
        chain_id = get_modified_chain_id_if_any_lowercase_letter(chain_id)
        mapping_fname = os.path.join(directories.pdb_fasta_mapping_dir, pdb_id + '_' + chain_id + '_invalid.rmsx.nch')
        if not os.path.isfile(mapping_fname):
            valid_chains.append(chain_id)

    for chain_id in valid_chains:
        mapping_fname = os.path.join(directories.pdb_fasta_mapping_dir, pdb_id + '_' + chain_id + '.rmsx.nch')
        if not os.path.isfile(mapping_fname):
            return False

    return True

def write_invalid_mapping_file(pdb_id, chain_id, directories):
    chain_id = get_modified_chain_id_if_any_lowercase_letter(chain_id)
    mapping_fname = os.path.join(directories.pdb_fasta_mapping_dir, pdb_id + '_' + chain_id + '_invalid.rmsx.nch')
    if not os.path.isfile(mapping_fname):
        fp = open(mapping_fname, 'w')
        fp.close()

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
    for k in range(len(aln_seq1)):
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

def get_valid_chain_list(pdb_id, chains, residue_dict, modified_residue_dict, ref_seq_dict, directories):
    valid_chains = []

    residue_abbreviation_dict = base_abbreviation(os.path.join(directories.lib_dir, 'baselist.dat'))     # baselist.bat downloaded from x3dna
    amino_acid_list = amino_acid_collection(os.path.join(directories.lib_dir, 'aminoacidlist.dat'))

    for chain_id in chains:
        if chain_id not in ref_seq_dict:
            logger.warning(pdb_id + '\tchain id: ' + chain_id + ' not found in ref_seq_dict. May be FASTA file is corrupted. Chain skipped.')
            write_invalid_mapping_file(pdb_id, chain_id, directories)
            continue

        if chain_id not in residue_dict:
            logger.warning(pdb_id + '\tchain_id: ' + chain_id + ' residue_list None. Chain skipped.')
            write_invalid_mapping_file(pdb_id, chain_id, directories)
            continue

        residue_dict[chain_id] = replace_modified_residue(pdb_id, chain_id, residue_dict[chain_id], modified_residue_dict, residue_abbreviation_dict, amino_acid_list)

        if len(residue_dict[chain_id]) == 0:
            logger.warning(pdb_id + '\tchain_id: ' + chain_id + ' residue_list Empty. Chain skipped.')
            write_invalid_mapping_file(pdb_id, chain_id, directories)
            continue

        valid_chains.append(chain_id)       #chains that could be processed properly

    return valid_chains

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
    
    residue_seq = ''.join(map(lambda x: x.symbol, residue_list))    # get the residue sequence

    # print(residue_seq)
    # print(ref_seq)
    # print('before aligning')

    # deprecated
    # aln = pairwise2.align.globalms(residue_seq, ref_seq, 5, -3, -10, -1)
    # (aln_residue, aln_ref, _, _, _) = aln[0]

    # new code
    aln_residue, aln_ref = get_seq_alignment(residue_seq, ref_seq)

    # print(aln_ref)
    # print(aln_residue)
    ref_seq_replaced = replace_unknown_letter_in_ref(aln_ref, aln_residue)
    # print(aln_residue)
    # print(ref_seq_replaced)
    residue_to_ref_mapping, ref_to_residue_mapping = get_aln_mapping(aln_residue, ref_seq_replaced)
    # print(residue_to_ref_mapping)
    # print(ref_to_residue_mapping)

    return residue_list, ref_seq_replaced, residue_to_ref_mapping, ref_to_residue_mapping

def get_pdbx_and_fasta_mapping_data(pdb_id, chains, directories):
    pdb_fn = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
    fasta_fn = os.path.join(directories.fasta_dir, pdb_id + '.fasta')

    if not os.path.isfile(pdb_fn) or not os.path.isfile(fasta_fn):
        logger.error(pdb_id + '.cif' + ' or ' + pdb_id + '.fasta' + ' file is missing.')
        return

    residue_dict, missing_residue_dict, modified_residue_dict = load_pdb_data_from_file(pdb_fn)
    residue_dict = insert_missing_residue(residue_dict, missing_residue_dict)

    # ref_seq_dict = {}
    # for record in SeqIO.parse(fasta_fn, 'fasta'):
    #     # ref_seq_dict[record.id.strip().split('|')[0].strip().split(':')[1]] = str(record.seq)
    #     # chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1].strip().split(',')
    #     # for chain_id in chain_ids:
    #     #     ref_seq_dict[chain_id] = str(record.seq)
    #     chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1:]
    #     for chain_id in chain_ids:
    #         chain_id = chain_id.strip().strip(',')
    #         if '[' in chain_id:
    #             continue
    #             # chain_id = chain_id.split('[')[0].strip()
    #         elif ']' in chain_id:
    #             chain_id = chain_id.split(']')[0].strip()
            
    #         ref_seq_dict[chain_id] = str(record.seq)
    ref_seq_dict = load_fasta_seq(pdb_id, chains, directories)

    chains = get_valid_chain_list(pdb_id, chains, residue_dict, modified_residue_dict, ref_seq_dict, directories)

    res_to_ref = {}
    ref_to_res = {}
    for chain_id in chains:

        residue_dict[chain_id], ref_seq_dict[chain_id], residue_to_ref_mapping, ref_to_residue_mapping = get_residue_reference_both_way_mapping_data_for_single_chain(residue_dict[chain_id], ref_seq_dict[chain_id])

        res_to_ref[chain_id] = residue_to_ref_mapping
        ref_to_res[chain_id] = ref_to_residue_mapping

    return chains, residue_dict, ref_seq_dict, res_to_ref, ref_to_res, missing_residue_dict

def generate_pdbx_fasta_mapping_files_for_single_pdb(pdb_id, chains, directories):
    pdb_fn = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
    fasta_fn = os.path.join(directories.fasta_dir, pdb_id + '.fasta')

    if not os.path.isfile(pdb_fn) or not os.path.isfile(fasta_fn):
        logger.error(pdb_id + '.cif' + ' or ' + pdb_id + '.fasta' + ' file is missing. Not generating the mapping file.')
        return

    residue_dict, missing_residue_dict, modified_residue_dict = load_pdb_data_from_file(pdb_fn)
    residue_dict = insert_missing_residue(residue_dict, missing_residue_dict)
    # residue_with_missing_base_dict = get_missing_residue(pdbx_dir + 'residue_with_missing_base_all/' + pdb_id + '.cif')

    # ref_seq_dict = {}
    # for record in SeqIO.parse(fasta_fn, 'fasta'):
    #     # ref_seq_dict[record.id.strip().split('|')[0].strip().split(':')[1]] = str(record.seq)
    #     # chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1].strip().split(',')
    #     chain_ids = record.description.strip().split('|')[1].strip().split(' ')[1:]
    #     for chain_id in chain_ids:
    #         chain_id = chain_id.strip().strip(',')
    #         if '[' in chain_id:
    #             continue
    #             # chain_id = chain_id.split('[')[0].strip()
    #         elif ']' in chain_id:
    #             chain_id = chain_id.split(']')[0].strip()
            
    #         ref_seq_dict[chain_id] = str(record.seq)
    ref_seq_dict = load_fasta_seq(pdb_id, chains, directories)
        
    # for generating multi-chain sequence files
    multi_seq_dict = {}
    multi_chain_mapping = {}

    chains = get_valid_chain_list(pdb_id, chains, residue_dict, modified_residue_dict, ref_seq_dict, directories)
    all_mapping_file_exists = True
    for chain_id in chains:
        chain_id = get_modified_chain_id_if_any_lowercase_letter(chain_id)
        mapping_fname = os.path.join(directories.pdb_fasta_mapping_dir, pdb_id + '_' + chain_id + '.rmsx.nch')
        if not os.path.isfile(mapping_fname):
            all_mapping_file_exists = False
            break
    
    if all_mapping_file_exists == True:
        return

    if len(chains) > 1:
        multi_chain_mapping_fname = os.path.join(directories.pdb_fasta_mapping_dir, pdb_id+'_'+'_'.join(chains)+'.rmsx.nch')
        fp_m = open(multi_chain_mapping_fname, 'w')

    for chain_id in chains:
        
        residue_dict[chain_id], ref_seq_dict[chain_id], residue_to_ref_mapping, ref_to_residue_mapping = get_residue_reference_both_way_mapping_data_for_single_chain(residue_dict[chain_id], ref_seq_dict[chain_id])
        
        multi_chain_mapping[chain_id] = {x.index: residue_to_ref_mapping[residue_dict[chain_id].index(x)] + sum(list(map(lambda x: len(multi_seq_dict[x]), multi_seq_dict))) for x in residue_dict[chain_id]}
        # multi_seq_dict[chain_id] = ref_seq_replaced
        multi_seq_dict[chain_id] = ref_seq_dict[chain_id]

        mapping_fname = os.path.join(directories.pdb_fasta_mapping_dir, pdb_id + '_' + get_modified_chain_id_if_any_lowercase_letter(chain_id) + '.rmsx.nch')
        fp = open(mapping_fname, 'w')
        
        for m in residue_to_ref_mapping:
            fp.write(str(residue_dict[chain_id][m].index)+'\t'+str(residue_to_ref_mapping[m])+'\n')
            
            if len(chains) > 1:
                fp_m.write(str(residue_dict[chain_id][m].index)+'\t'+str(multi_chain_mapping[chain_id][residue_dict[chain_id][m].index])+'\n')

        fp.close()

    if len(chains) > 1:
        fp_m.close()

def generate_pdbx_fasta_mapping_files(pdb_chains, mp_number_of_process, directories):
    if len(pdb_chains) == 0:
        base_path_len = len(os.path.dirname(directories.base_dir)) + 1
        logger.info('Using existing PDBx-FASTA mapping files in ' + directories.pdb_fasta_mapping_dir[base_path_len:] + '.')
        print('')
        return

    logger.info('Generating PDBx-FASTA mapping files.')
    start_time = time.time()

    parameter_list = []
    for pdb_id in pdb_chains:
        # print(pdb_id, pdb_chains[pdb_id])
        if all_mapping_file_exists_for_single_pdb(pdb_id, pdb_chains[pdb_id], directories) == False:
            # generate_pdbx_fasta_mapping_files_for_single_pdb(pdb_id, pdb_chains[pdb_id])
            parameter_list.append((pdb_id, pdb_chains[pdb_id], directories))

    pool = mp.Pool(mp_number_of_process)
    pool.map(_mapping_worker, parameter_list)

    wait_for_certain_time_according_to_wait_factor(len(parameter_list))
    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

def _mapping_worker(p):
    generate_pdbx_fasta_mapping_files_for_single_pdb(*p)

def _dssr_annotation_worker(p):
    generate_dssr_annotation_for_single_pdb(*p)

def _fr3d_annotation_worker(p):
    download_fr3d_annotation_for_single_pdb(*p)

def _dssr_and_fr3d_merged_annotation_worker(p):
    get_dssr_and_fr3d_merged_annotation_for_single_pdb(*p)

def get_annotation_files(pdb_list, directories, annotation_source, fr3d_url, dssr_url, mp_number_of_process, env='global'):
    selected_annotation_dir = os.path.join(directories.annotation_dir, annotation_source)
    create_directory(selected_annotation_dir)

    if len(pdb_list) == 0:
        base_path_len = len(os.path.dirname(directories.base_dir)) + 1
        logger.info('Using existing annotation files in ' + selected_annotation_dir[base_path_len:] + '.')
        print('')
        return

    annotation_file_ext = ''
    if annotation_source not in ['merged', 'dssr', 'fr3d']:
        logger.error('Invalid annotation_source provided in config file. Given "' + annotation_source + '", expected "merged/dssr/fr3d".')
        sys.exit()
        return
    
    parameter_list = []
    start_time = time.time()

    if annotation_source.lower() == 'dssr':
        logger.info('Getting DSSR annotation files.')
        
        for pdb_id in pdb_list:
            # logger.info('Getting DSSR annotation file for ' + pdb_id)
            annotation_fname = os.path.join(selected_annotation_dir, pdb_id + '.' + annotation_source)
            if os.path.isfile(annotation_fname):
                # logger.info('File exists. Skipping to next item.')
                continue
            else:
                parameter_list.append((selected_annotation_dir, pdb_id, dssr_url, directories, annotation_source, env))

        pool = mp.Pool(mp_number_of_process)
        pool.map(_dssr_annotation_worker, parameter_list)

    elif annotation_source.lower() == 'fr3d':
        logger.info('Getting FR3D annotation files.')

        for pdb_id in pdb_list:
            # logger.info('Getting FR3D annotation file for ' + pdb_id)
            annotation_fname = os.path.join(selected_annotation_dir, pdb_id + '.' + annotation_source)
            if os.path.isfile(annotation_fname):
                # logger.info('File exists. Skipping to next item.')
                continue
            else:
                parameter_list.append((selected_annotation_dir, pdb_id, fr3d_url, annotation_source))

        pool = mp.Pool(mp_number_of_process)
        pool.map(_fr3d_annotation_worker, parameter_list)

    else:
        create_directory(os.path.join(directories.annotation_dir, 'dssr'))
        create_directory(os.path.join(directories.annotation_dir, 'fr3d'))

        # Should we re-generate the stat? If so, proceed to this function. A flag can be used here.
        regenerate_ann_stats = False
        if regenerate_ann_stats:
            generate_ann_stats_before_merging(directories)

        logger.info('Getting DSSR and FR3D merged annotation files.')

        for pdb_id in pdb_list:
            # logger.info('Getting FR3D annotation file for ' + pdb_id)
            annotation_fname = os.path.join(selected_annotation_dir, pdb_id + '.' + annotation_source)
            if os.path.isfile(annotation_fname):
                # logger.info('File exists. Skipping to next item.')
                continue
            else:
                parameter_list.append((selected_annotation_dir, pdb_id, directories, annotation_source, fr3d_url, dssr_url, env))
                # print(pdb_id)
                # get_dssr_and_fr3d_merged_annotation_for_single_pdb(selected_annotation_dir, pdb_id, directories, annotation_source, fr3d_url, dssr_url, env)

        pool = mp.Pool(mp_number_of_process)
        pool.map(_dssr_and_fr3d_merged_annotation_worker, parameter_list)

    wait_for_certain_time_according_to_wait_factor(len(parameter_list))
    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

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

    # print(ref_to_residue_mapping)
    # print(residue_to_ref_mapping)
    for i, j in loop_end:
        loop_end_residue.append((residue_list[ref_to_residue_mapping[i]].index, residue_list[ref_to_residue_mapping[j]].index))

    for r in loop:
        for v in range(r[0], r[1]+1):
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

def get_loop_end(loop):
    ret = []
    for i, r in enumerate(loop):
        ret.append((r[0], loop[i-1][1]))
    return ret

def get_loop_type(loop):
    _, regions = loop.strip().split(':')
    regions = regions.strip().split('_')
    region_cnt = len(regions)
    if region_cnt == 1:
        return 'HL'
    elif region_cnt == 2:
        return 'IL'
    elif region_cnt > 2:
        return 'ML'
    logger.error('Invalid loop')
    return ''

def generate_loop_files_for_single_pdb(pdb_id, chains, loops, directories, annotation_source, env):
    chains, residue_dict, ref_seq_dict, res_to_ref, ref_to_res, missing_residue_dict = get_pdbx_and_fasta_mapping_data(pdb_id, chains, directories)

    is_detailed_ann = True
    bp_item_len = 5
    stk_item_len = 4
    if annotation_source.lower() == 'merged':
        annotation_list = parseMergedAnnotation(os.path.join(directories.annotation_dir, annotation_source.lower(), pdb_id + '.' + annotation_source.lower()), is_detailed_ann)
    elif annotation_source.lower() == 'dssr':
        annotation_list = parseDSSR(os.path.join(directories.annotation_dir, annotation_source.lower(), pdb_id + '.' + annotation_source.lower()), is_detailed_ann)
    elif annotation_source.lower() == 'fr3d':
        annotation_list = parseFR3D(os.path.join(directories.annotation_dir, annotation_source.lower(), pdb_id + '.' + annotation_source.lower()), is_detailed_ann)
    else:
        logger.error('Wrong annotation source provided. Please choose from Merged/DSSR/FR3D.')
        sys.exit()

    for item in loops:
        pdb_chain, regions = item.strip().split(':')
        _, chain_id = pdb_chain.strip().split('_')
        regions = regions.strip().split('_')

        # if chain_id not in chains:
        #     logger.error('Chain id ' + str(chain_id) + ' for pdb id ' + pdb_id + ' is somehow invalid. Skipping ' + item)

        loop = []
        for region in regions:
            s, e = region.strip().split('-')
            loop.append((int(s), int(e)))

        type = get_loop_type(item)
        residue_list = residue_dict[chain_id]
        ref_to_residue_mapping = ref_to_res[chain_id]
        residue_to_ref_mapping = res_to_ref[chain_id]

        loop_end = get_loop_end(loop)
        for i in range(len(loop)):
            loop_i = rotate(loop, i)
            loop_i_index = reduce(lambda y, z: y+z, map(lambda x: list(range(x[0], x[1]+1)), loop_i))
            
            # loop_intera = util.get_loop_intera(loop_i, residue_list, ref_to_residue_mapping, residue_to_ref_mapping, annotation_list)
            loop_intera = get_loop_intera(loop_i, loop_end, residue_list, ref_to_residue_mapping, residue_to_ref_mapping, annotation_list, is_detailed_ann)
            loop_bp_intera = filter(lambda x: len(x[2]) == bp_item_len, loop_intera)
            loop_stack_intera = filter(lambda x: len(x[2]) == stk_item_len, loop_intera)
            loop_bp_intera = sorted(loop_bp_intera, key=lambda x: x[0])
            loop_stack_intera = sorted(loop_stack_intera, key=lambda x: x[0])

            # loop_fn_smf = os.path.join(directories.loop_dir, pdb_id+'_%s:%s.smf' % (chain_id, '_'.join(map(lambda x: str(x[0])+'-'+str(x[1]), loop_i))))  #for scanx with stacking
            loop_fn_smf = os.path.join(directories.loop_dir, pdb_id+'_%s_%s.smf' % (chain_id, '_'.join(map(lambda x: str(x[0])+'-'+str(x[1]), loop_i))))  #for scanx with stacking

            if env == 'local':
                fms = open('missing_residue_stat.txt', 'a')

            if chain_id in missing_residue_dict:# or chain_id in residue_with_missing_base_dict:
                missing_res_segment_lengths = []
                # if chain_id in residue_with_missing_base_dict:
                #     missing_res_segment_lengths += generate_missing_residue_statistics(loop_i_index, residue_with_missing_base_dict[chain_id], residue_list, ref_to_residue_mapping)
                    # if len(missing_res_segment_lengths) > 0:
                    #     print missing_res_segment_lengths
                if chain_id in missing_residue_dict:
                    missing_res_segment_lengths += generate_missing_residue_statistics(loop_i_index, missing_residue_dict[chain_id], residue_list, ref_to_residue_mapping)
                if len(missing_res_segment_lengths) > 0:
                    # loop_with_missing_res_dict[pdb_id, chain_id, type, str(loop_i)] = (len(loop_i_index), missing_res_segment_lengths)
                    if env == 'local':
                        fms.write(pdb_id + '\t' + chain_id + '\t' + type + '\t' + str(loop) + '\t' + str(len(loop_i_index)))
                    missing_res_length = 0
                    for item in missing_res_segment_lengths:
                        if env == 'local':
                            fms.write('\t' + str(item))
                        missing_res_length += item
                    if env == 'local':
                        fms.write('\n')

                    if missing_res_length > 5 or 2 * missing_res_length > len(loop_i_index):
                        create_directory(os.path.join(directories.loop_dir, 'missing_res/' + type))
                        # loop_fn_smf = os.path.join(directories.loop_dir, 'missing_res/' + type, pdb_id+'_%s:%s.smf' % (chain_id, '_'.join(map(lambda x: str(x[0])+'-'+str(x[1]), loop_i))))  #for scanx with stacking
                        loop_fn_smf = os.path.join(directories.loop_dir, 'missing_res/' + type, pdb_id+'_%s_%s.smf' % (chain_id, '_'.join(map(lambda x: str(x[0])+'-'+str(x[1]), loop_i))))  #for scanx with stacking
            if env == 'local':
                fms.close()

            fs = open(loop_fn_smf, 'w')

            fs.write('>%s_%s:%s\n' % (pdb_id, chain_id, '_'.join(map(lambda x: str(x[0])+'-'+str(x[1]), loop_i))))
            fs.write('...'.join(map(lambda x: ref_seq_dict[chain_id][x[0]:x[1]+1], loop_i))+'\n')
            fs.write('#info=basepair\n')

            for i, j, intera in loop_bp_intera:
                if intera[3] == 'hbond':
                    continue
                if i < j:
                    fs.write('%d-%d,%s,%s\n' % (i, j, intera[2], intera[3]))
                else:
                    fs.write('%d-%d,%s,%s\n' % (j, i, intera[2][::-1], intera[3]))
            
            fs.write('#info=stacking\n')
            for i, j, intera in loop_stack_intera:
                if i < j:
                    fs.write('%d-%d,%s\n' % (i, j, intera[2]))
                else:
                    if intera[2] == 'upward':
                        fs.write('%d-%d,%s\n' % (j, i, 'downward'))
                    elif intera[2] == 'downward':
                        fs.write('%d-%d,%s\n' % (j, i, 'upward'))
                    else:
                        fs.write('%d-%d,%s\n' % (j, i, intera[2]))

            fs.close()

def _generate_loop_files_worker(p):
    generate_loop_files_for_single_pdb(*p)

def generate_loop_files(pdb_loops, pdb_chains, directories, annotation_source, mp_number_of_process, env):
    if len(pdb_chains) == 0:
        base_path_len = len(os.path.dirname(directories.base_dir)) + 1
        logger.info('Using existing loop files (loop.smf) in ' + directories.loop_dir[base_path_len:] + '.')
        print('')
        return

    logger.info('Generating loop files (loop.smf) with base-pairs and stackings according to annotation.')

    start_time = time.time()
    if env == 'local':
        fms = open('missing_residue_stat.txt', 'w')
        fms.close()

    parameter_list = []
    for pdb_id in pdb_loops:
        chains = pdb_chains[pdb_id]
        loops = pdb_loops[pdb_id]
        # print(pdb_id, chains, len(loops))
        # generate_loop_files_for_single_pdb(pdb_id, chains, loops, directories, annotation_source, env)
        parameter_list.append((pdb_id, chains, loops, directories, annotation_source, env))
        
    pool = mp.Pool(mp_number_of_process)
    pool.map(_generate_loop_files_worker, parameter_list)

    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

def is_all_data_available(pdb, chain, directories, annotation_source):

    if not os.path.isfile(os.path.join(directories.pdbx_dir, pdb + '.cif')):
        return False
    if not os.path.isfile(os.path.join(directories.fasta_dir, pdb + '.fasta')):
        return False
    chain = get_modified_chain_id_if_any_lowercase_letter(chain)
    if not os.path.isfile(os.path.join(directories.pdb_fasta_mapping_dir, pdb + '_' + chain + '.rmsx.nch')):
        if not os.path.isfile(os.path.join(directories.pdb_fasta_mapping_dir, pdb + '_' + chain + '_invalid.rmsx.nch')):
            return False
    if not os.path.isfile(os.path.join(directories.annotation_dir, annotation_source, pdb + '.' + annotation_source)):
        return False

    return True

def get_all_loop_combination(loop):
    loop_combinations = []
    pdb_chain, regions = loop.strip().split(':')
    regions = regions.strip().split('_')
    loop = []
    for region in regions:
        s, e = region.strip().split('-')
        loop.append((s, e))
    
    for i in range(len(loop)):
        loop_i = rotate(loop, i)
        loop_combinations.append(pdb_chain + ':' + '_'.join(list(map(lambda x: '-'.join(x), loop_i))))

    # print(loop_combinations)
    return loop_combinations

def prepare_loop_files(loop_node_list_str, directories, annotation_source, mp_number_of_process, env):

    loop_list_to_generate = []
    for node in loop_node_list_str:
        loop_rotations = get_all_loop_combination(str(node))
        has_all_rotations = True
        for loop in loop_rotations:
            loop_fn_smf = os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf')
            if not os.path.isfile(loop_fn_smf):
                has_all_rotations = False
                break
        if has_all_rotations == False:
            loop_list_to_generate.append(str(node))

    pdb_chains = {}
    pdb_loops = {}
    
    for loop in loop_list_to_generate:
        pdb_chain, loop_region = loop.strip().split(':')
        pdb, chain = pdb_chain.strip().split('_')
        if pdb not in pdb_chains:
            pdb_chains[pdb] = []
        if pdb not in pdb_loops:
            pdb_loops[pdb] = []
        if chain not in pdb_chains[pdb]:
            pdb_chains[pdb].append(chain)
        if loop not in pdb_loops[pdb]:
            pdb_loops[pdb].append(loop)

    pdb_count = len(pdb_chains)

    generate_loop_files(pdb_loops, pdb_chains, directories, annotation_source, mp_number_of_process, env)  # generate mapping file first and the proceed to generating loop files
    wait_for_certain_time_according_to_wait_factor(pdb_count)

def prepare_data(families, directories, annotation_source, content_download_attempts, mp_number_of_process):
    pdb_chains = {}
    # pdb_loops = {}
    for family in families:
        loops = families[family]
        for loop in loops:
            pdb_chain, loop_region = loop.strip().split(':')
            pdb, chain = pdb_chain.strip().split('_')
            if is_all_data_available(pdb, chain, directories, annotation_source):
                continue
            if pdb not in pdb_chains:
                pdb_chains[pdb] = []
            if chain not in pdb_chains[pdb]:
                pdb_chains[pdb].append(chain)

    pdb_list = pdb_chains.keys()
    pdb_count = len(pdb_list)

    pdbx_url, fasta_url, fr3d_url, dssr_url = get_urls()

    get_pdbx_and_fasta_files(pdb_list, pdbx_url, fasta_url, directories, content_download_attempts)
    get_annotation_files(pdb_list, directories, annotation_source, fr3d_url, dssr_url, mp_number_of_process, get_env())
    generate_pdbx_fasta_mapping_files(pdb_chains, mp_number_of_process, directories)
    # print('Done till here')
    # sys.exit()

def prepare_data_for_a_chain(pdb_chain_to_search, directories, annotation_source, content_download_attempts, mp_number_of_process):
    pdb_id, chain_id = pdb_chain_to_search.strip().split('_')
    # chain_id = get_modified_chain_id_if_any_lowercase_letter(chain_id)
    pdb_chains = {}
    pdb_chains[pdb_id] = []
    pdb_chains[pdb_id].append(chain_id)
    if is_all_data_available(pdb_id, chain_id, directories, annotation_source) == False:
        pdb_list = []
        pdb_list.append(pdb_id)
        pdbx_url, fasta_url, fr3d_url, dssr_url = get_urls()
        get_pdbx_and_fasta_files(pdb_list, pdbx_url, fasta_url, directories, content_download_attempts)
        get_annotation_files(pdb_list, directories, annotation_source, fr3d_url, dssr_url, mp_number_of_process, get_env())
        generate_pdbx_fasta_mapping_files(pdb_chains, mp_number_of_process, directories)

def get_loops_from_a_chain(pdb_chain_to_search, directories, ann_source, junction_count):
    pdb_id, chain_id = pdb_chain_to_search.strip().split('_')
    outlier_zscore_threshold = 2.5
    remove_empty_loops = True

    pdb_chain_dict = {}
    pdb_chain_dict[pdb_id] = []
    pdb_chain_dict[pdb_id].append(chain_id)

    junction_wise_filtered_loops = cut_loops(pdb_chain_dict, directories, ann_source, remove_empty_loops, outlier_zscore_threshold)

    return junction_wise_filtered_loops[junction_count]
