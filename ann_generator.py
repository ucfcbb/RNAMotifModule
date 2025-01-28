import sys
import os
import logging
import requests

from future.standard_library import install_aliases
install_aliases()
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request, urlretrieve
from urllib.error import HTTPError

# sys.path.append('../../')
# from config import *
# sys.path.append(scripts_dir)
from ann_parser import *

def generate_dssr_annotation_for_single_pdb(selected_annotation_dir, pdb_id, dssr_url, directories, annotation_source, env):
    output_env = env
    if output_env == 'global':
        download_dssr_annotation_for_single_pdb(selected_annotation_dir, pdb_id, dssr_url)
        return
        # if annotation_source == 'dssr':
        #     logger.error('DSSR annotation file not exists for ' + pdb_id + '.')
        #     sys.exit()
        # else:
        #     # logger.warning('DSSR annotation file not exists for ' + pdb_id + '.')
        #     return

    pdb_fname = os.path.join(directories.pdbx_dir, pdb_id + '.cif')
    t_dir = os.path.join(directories.temp_dir, pdb_id)
    create_directory(t_dir)
    os.chdir(t_dir)
    # os.system("./x3dna-dssr -i=%s -o=%s --non-pair" % (pdb_fname, os.path.join(selected_annotation_dir, pdb_id + '.dssr')))
    os.system("%s -i=%s -o=%s --non-pair" % (os.path.join(directories.dssr_dir, "x3dna-dssr"), pdb_fname, os.path.join(selected_annotation_dir, pdb_id + '.dssr')))

def download_dssr_annotation_for_single_pdb(selected_annotation_dir, pdb_id, dssr_url):
    try:
        response = requests.get(dssr_url.replace('XXXX', pdb_id.lower()), allow_redirects=False)

        if response.status_code == 200:
            fp = open(os.path.join(selected_annotation_dir, pdb_id+'.dssr'), 'wb')
            fp.write(response.content)
            fp.close()
        else:
            logger.error('Error downloading dssr annotation file for ' + pdb_id + '. ' + e.Read())
            logger.info('Try generating DSSR annotation for ' + pdb_id + ' using DSSR tool or from DSSR website \'http://skmatic.x3dna.org/\'.')
        # logger.info('File downloaded successfully: ' + os.path.join(selected_annotation_dir, pdb_id+'.fr3d'))
    except HTTPError as e:
        if annotation_source == 'dssr':
            logger.error('Error downloading dssr annotation file for ' + pdb_id + '. ' + e.Read())
            logger.info('Try generating DSSR annotation for PDB: ' + pdb_id + ' using DSSR tool or from DSSR website \'http://skmatic.x3dna.org/\'.')
            sys.exit()
        else:
            # logger.warning('Error downloading fr3d annotation file for ' + pdb_id + '. ' + e.Read())
            return

def download_fr3d_annotation_for_single_pdb(selected_annotation_dir, pdb_id, fr3d_url, annotation_source):
    try:
        csvfile = urlopen(fr3d_url.replace('XXXX', pdb_id))
        fp = open(os.path.join(selected_annotation_dir, pdb_id+'.fr3d'), 'wb')
        fp.write(csvfile.read())
        fp.close()
        # logger.info('File downloaded successfully: ' + os.path.join(selected_annotation_dir, pdb_id+'.fr3d'))
    except HTTPError as e:
        if annotation_source == 'fr3d':
            logger.error('Error downloading fr3d annotation file for ' + pdb_id + '. ' + e.Read())
            sys.exit()
        else:
            # logger.warning('Error downloading fr3d annotation file for ' + pdb_id + '. ' + e.Read())
            return

def get_highest_freq_interactions(interactions):
    interact_dict = {}
    for interact in interactions:
        if interact not in interact_dict:
            interact_dict[interact] = 0
        interact_dict[interact] += 1

    max_freq = 0
    for interact in interact_dict:
        max_freq = max(interact_dict[interact], max_freq)

    max_freq_interactions = []
    for interact in interact_dict:
        if interact_dict[interact] == max_freq:
            max_freq_interactions.append(interact)

    return max_freq_interactions

def get_resolved_dict(pdb_id, pairwise_dict, interaction_category_rank):
    resolved_dict = {}
    for ind_pair in pairwise_dict:
        interactions = get_highest_freq_interactions(pairwise_dict[ind_pair])
        if len(interactions) > 1:
            best_rank = 1000
            best_rank_interaction = ("", "")
            for bp, interact in interactions:
                if bp not in interaction_category_rank:
                    logger.error('ERROR: Problem in ' + pdb_id + '. Check bp ' + bp + '. Skipping.')
                    continue
                    sys.exit()
                if interact not in interaction_category_rank[bp]:
                    rank = len(interaction_category_rank[bp])
                else:
                    rank = interaction_category_rank[bp].index(interact)

                if rank < best_rank:
                    best_rank = rank
                    best_rank_interaction = (bp, interact)
            interactions = []
            interactions.append((best_rank_interaction[0], best_rank_interaction[1]))

        resolved_dict[ind_pair] = interactions

    return resolved_dict

def write_merged_annotation_to_file(pdb_id, ann_list, merged_annotation_dir, detailed):
    if detailed:
        bp_item_len = 5
        stk_item_len = 4
    else:
        logger.error('Detailed BP info not available.')
        sys.exit()

    # for pdb_id in ann_dict:
    fp = open(os.path.join(merged_annotation_dir, pdb_id + ".merged"), "w")
    for annotation in ann_list:
        if len(annotation) == bp_item_len:
            index1, index2, edges, orientation, bp = annotation
            fp.write(str(index1) + '\t' + str(index2) + '\t' + bp + '\t' + edges + '\t' + orientation + '\n')

        elif len(annotation) == stk_item_len:
            index1, index2, direction, bp = annotation
            fp.write(str(index1) + '\t' + str(index2) + '\t' + bp + '\t' + direction + '\n')

        else:
            logger.error('ERROR: invalid interact_info length in ' + pdb_id)
            sys.exit()
    fp.close()

def resolve_pairwise_annotation_conflict_helper(pdb_id, ann_list, bp_interaction_category_rank, stk_interaction_category_rank, merged_annotation_dir, detailed):
    if detailed:
        bp_item_len = 5
        stk_item_len = 4
    else:
        logger.error('Detailed BP info not available.')
        sys.exit()

    pairwise_bp_dict = {}
    pairwise_stk_dict = {}

    logger.info('Resolving ' + pdb_id)
    for interact_info in ann_list:
        if len(interact_info) == bp_item_len:
            index1, index2, edges, orientation, bp = interact_info
            interaction = orientation[0] + edges[0] + edges[2]
            if (index1, index2) not in pairwise_bp_dict:
                pairwise_bp_dict[(index1, index2)] = []
            pairwise_bp_dict[(index1, index2)].append((bp, interaction))

        elif len(interact_info) == stk_item_len:
            index1, index2, direction, bp = interact_info
            if (index1, index2) not in pairwise_stk_dict:
                pairwise_stk_dict[(index1, index2)] = []
            pairwise_stk_dict[(index1, index2)].append((bp, direction))

        else:
            logger.error('ERROR: invalid interact_info length in ' + pdb_id)
            sys.exit()

    pairwise_bp_dict = get_resolved_dict(pdb_id, pairwise_bp_dict, bp_interaction_category_rank)
    pairwise_stk_dict = get_resolved_dict(pdb_id, pairwise_stk_dict, stk_interaction_category_rank)

    # if pdb_id == '4V5G':
    # print('Creating resolved_annotation_list')

    resolved_annotation_list = []
    for (index1, index2) in pairwise_bp_dict:
        (bp, interaction) = pairwise_bp_dict[(index1, index2)][0]
        orientation = 'cis' if interaction[0] == 'c' else 'trans'
        edges = interaction[1] + '/' + interaction[2]
        resolved_annotation_list.append((index1, index2, edges, orientation, bp))

    for index1, index2 in pairwise_stk_dict:
        bp, direction = pairwise_stk_dict[(index1, index2)][0]
        resolved_annotation_list.append((index1, index2, direction, bp))

    resolved_annotation_list = sorted(resolved_annotation_list)
    # if pdb_id == '4V5G':
    # print('writing to file')
    write_merged_annotation_to_file(pdb_id, resolved_annotation_list, merged_annotation_dir, detailed)

def load_category_rank_from_file(fname):
    category_rank = {}

    fp = open(fname)
    rank_list = csv_to_list(fp.readlines())
    fp.close()

    for item in rank_list:
        category_rank[item[0]] = item[1:]

    return category_rank

def get_dssr_and_fr3d_merged_annotation_for_single_pdb(selected_annotation_dir, pdb_id, directories, annotation_source, fr3d_url, dssr_url, env):
    dssr_ann_dir = os.path.join(directories.annotation_dir, 'dssr')
    fr3d_ann_dir = os.path.join(directories.annotation_dir, 'fr3d')
    dssr_fname = os.path.join(dssr_ann_dir, pdb_id + '.dssr')
    fr3d_fname = os.path.join(fr3d_ann_dir, pdb_id + '.fr3d')
    if not os.path.isfile(dssr_fname):
        generate_dssr_annotation_for_single_pdb(dssr_ann_dir, pdb_id, dssr_url, directories, annotation_source, env)
    if not os.path.isfile(fr3d_fname):
        download_fr3d_annotation_for_single_pdb(fr3d_ann_dir, pdb_id, fr3d_url, annotation_source)

    detailed_info = True
    ann_list = parseDSSR(dssr_fname, detailed_info)
    
    if os.path.isfile(fr3d_fname):
        ann_list.extend(parseFR3D(fr3d_fname, detailed_info))

    if len(ann_list) == 0:
        logger.error('None of the DSSR and FR3D annotation files were found. Cannot generate the merged annotation file for ' + pdb_id + '.')
        sys.exit()

    bp_interaction_category_rank = load_category_rank_from_file(os.path.join(directories.lib_dir, 'bp_category_rank.csv'))
    stk_interaction_category_rank = load_category_rank_from_file(os.path.join(directories.lib_dir, 'stk_category_rank.csv'))

    resolve_pairwise_annotation_conflict_helper(pdb_id, ann_list, bp_interaction_category_rank, stk_interaction_category_rank, selected_annotation_dir, detailed_info)
    # if pdb_id == "4V5G":
    # print('resolved')

def initialize_interaction_dict(bp_category_dict, stk_category_dict):
    bp_list = []
    interaction_list = []

    res_list = ['A', 'C', 'G', 'U']
    edge_list = ['W', 'H', 'S']
    orientation_list = ['c', 't']
    stack_direction_list = ["upward", "downward", "inward", "outward"]

    for res1 in res_list:
        for res2 in res_list:
            bp_list.append(res1 + "-" + res2)

    for orientation in orientation_list:
        for edge1 in edge_list:
            for edge2 in edge_list:
                interaction_list.append(orientation + edge1 + edge2)

    for bp in bp_list:
        bp_category_dict[bp] = {}
        stk_category_dict[bp] = {}
        for interaction in interaction_list:
            bp_category_dict[bp][interaction] = 0
        for direction in stack_direction_list:
            stk_category_dict[bp][direction] = 0

    # print_a_dict_sorted(bp_category_dict)
    # print_a_dict_sorted(stk_category_dict)
    # sys.exit()
    return bp_list, interaction_list

def add_to_category_dict(category_dict, bp, interaction):
    if bp not in category_dict:
        category_dict[bp] = {}
    if interaction not in category_dict[bp]:
        category_dict[bp][interaction] = 0

    category_dict[bp][interaction] += 1

def generate_category_stat(pdb_id, annotation_list, bp_category_dict, stk_category_dict, detailed):
    # global bp_category_dict
    # global stk_category_dict
    bp_cnt = 0
    stk_cnt = 0

    annotation_list = list(set(annotation_list))
    if detailed:
        bp_item_len = 5
        stk_item_len = 4
    else:
        print("Detailed BP info not available.")
        sys.exit()

    for interact_info in annotation_list:
        if len(interact_info) == bp_item_len:
            (index1, index2, edges, orientation, bp) = interact_info
            rev_edges, rev_bp = get_inverse_bp_info(edges, bp)
            edges = orientation[0] + edges[0] + edges[2]
            rev_edges = orientation[0] + rev_edges[0] + rev_edges[2]
            add_to_category_dict(bp_category_dict, bp, edges)
            # if not (bp == rev_bp and edges == rev_edges):
            add_to_category_dict(bp_category_dict, rev_bp, rev_edges)
                
            bp_cnt += 1

        elif len(interact_info) == stk_item_len:
            (index1, index2, direction, bp) = interact_info
            rev_direction, rev_bp = get_inverse_stk_info(direction, bp)
            add_to_category_dict(stk_category_dict, bp, direction)
            # if not (bp == rev_bp and direction == rev_direction):
            add_to_category_dict(stk_category_dict, rev_bp, rev_direction)

            stk_cnt += 1

        else:
            print("ERROR: invalid item length in " + pdb_id)
            sys.exit()

    return bp_cnt, stk_cnt

def write_interaction_category_rank_file(category_dict, fn):
    # print_a_dict_sorted(category_dict)
    fp = open(fn, "w")
    for bp in sorted(category_dict):
        # print bp + ": ",
        fp.write(bp)
        for interaction in sorted(category_dict[bp], key=lambda x: category_dict[bp][x], reverse=True):
            if category_dict[bp][interaction] > 0:
                # print interaction + "-" + str(category_dict[bp][interaction]),
                fp.write("," + interaction)
        # print ""
        fp.write("\n")
    fp.close()

def write_interaction_dict_in_file(category_dict, fn):
    # print_a_dict_sorted(category_dict)
    fp = open(fn, "w")
    for bp in sorted(category_dict):
        print(bp + ": ", end=",")
        fp.write(bp + ": ")
        for interaction in sorted(category_dict[bp], key=lambda x: category_dict[bp][x], reverse=True):
            if category_dict[bp][interaction] > 0:
                print(interaction + "-" + str(category_dict[bp][interaction]), end=",")
                fp.write(interaction + "-" + str(category_dict[bp][interaction]) + " ")
        print("")
        fp.write("\n")
    fp.close()

def write_bp_matrix(bp_dict, fn, bp_list, interaction_list):
    # matrix_fn = fn[:-4] + "_matrix.txt"
    fx = open(fn, "w")
    for interaction in interaction_list:
        fx.write("\t" + interaction)
    fx.write("\n")
    for bp in bp_list:
        fx.write(bp)
        for interaction in interaction_list:
            fx.write("\t" + str(bp_dict[bp][interaction]))
        fx.write("\n")
    fx.close()

def generate_ranking_from_interaction_dict(category_dict):
    interaction_category_rank = {}
    for bp in category_dict:
        interaction_category_rank[bp] = sorted(category_dict[bp], key=lambda x: category_dict[bp][x], reverse=True)

    return interaction_category_rank

def generate_ann_stats_before_merging(directories):

    bp_rank_fn = os.path.join(directories.lib_dir, "bp_category_rank.csv")
    stk_rank_fn = os.path.join(directories.lib_dir, "stk_category_rank.csv")

    rename_filename_if_exists(bp_rank_fn)
    rename_filename_if_exists(stk_rank_fn)

    dssr_dir = os.path.join(directories.annotation_dir, 'dssr')
    fr3d_dir = os.path.join(directories.annotation_dir, 'fr3d')
    # merged_annotation_dir = "../../Annotation/dssr_fr3d_merged/"

    ann_dict = {}
    bp_category_dict = {}
    stk_category_dict = {}
    total_bp = 0
    total_stk = 0
    # bp_list = []
    # interaction_list = []

    total_merged_interact = 0
    total_dssr_interact = 0
    total_fr3d_interact = 0

    bp_list, interaction_list = initialize_interaction_dict(bp_category_dict, stk_category_dict)
    detailed_info = True
    # test_list = ["4V9F.cif.out", "4WF9.cif.out"]
    # for fn in test_list:
    for fn in glob.glob(os.path.join(dssr_dir, "*.dssr")):
        # fn = os.path.join(dssr_dir, "4V9F.cif.out")
        # fn = os.path.join(dssr_dir, "5K8H.cif.out")
        # fn = os.path.join(dssr_dir, "5AFI.cif.out")
        # fn = os.path.join(dssr_dir, "4R4V.cif.out")
        # fn = os.path.join(dssr_dir, "4WF9.cif.out")
        # fn = os.path.join(dssr_dir, fn)
        
        print("parsing %s" % fn)
        pdb_id = os.path.basename(fn)[:-5]
        ann_dict[pdb_id] = parseDSSR(fn, detailed_info)
        total_dssr_interact += len(ann_dict[pdb_id])
        # bp_ann_list, stk_ann_list = generate_annotation_list(ann_dict[pdb_id], detailed_info)

        # if prefix != "dssr":
        fn = os.path.join(fr3d_dir, pdb_id + ".fr3d")
        if os.path.isfile(fn):

            fr3d_ann_list = parseFR3D(fn, detailed_info)
            total_fr3d_interact += len(fr3d_ann_list)
            # bp_fr3d_ann_list, stk_fr3d_ann_list = generate_annotation_list(fr3d_ann_dict, detailed_info)
            ann_dict[pdb_id].extend(fr3d_ann_list)

            # if prefix == "fr3d":
            #     # bp_ann_list = bp_fr3d_ann_list
            #     # stk_ann_list = stk_fr3d_ann_list
            #     ann_dict[pdb_id] = fr3d_ann_list
            # else:
            #     # bp_ann_list.extend(bp_fr3d_ann_list)
            #     # stk_ann_list.extend(stk_fr3d_ann_list)
            #     ann_dict[pdb_id].extend(fr3d_ann_list)

        else:
            print("FR3D annotation not found for " + pdb_id)

        
        # ann_dict[pdb_id] = list(set(ann_dict[pdb_id]))
        # ann_dict[pdb_id] = sorted(ann_dict[pdb_id])

        bp_cnt, stk_cnt = generate_category_stat(pdb_id, ann_dict[pdb_id], bp_category_dict, stk_category_dict, detailed_info)
        total_bp += bp_cnt
        total_stk += stk_cnt
        total_merged_interact += len(ann_dict[pdb_id])
        
        # break


    # if prefix == "fr3d":
    #     # print "Total PDB\t" + str(len(glob.glob(os.path.join(fr3d_dir, "*.fr3d"))))
    #     total_bp /= 2
    #     total_stk /= 2
    # else:
    print("Total PDB\t" + str(len(glob.glob(os.path.join(dssr_dir, "*.dssr")))))
    # print "Total bp\t" + str(total_bp)
    # print "Total stack\t" + str(total_stk)
    print("Total DSSR interaction\t" + str(total_dssr_interact))
    print("Total FR3D interaction\t" + str(total_fr3d_interact))
    print("Total merged interaction\t" + str(total_merged_interact))
    print("Total merged (bp, stack)\t" + str(total_bp) + "\t" + str(total_stk))
    print("\n")
    # for e in edge_dict:
    #     print e
    #     print_a_dict(edge_dict[e])
    #     print ""
    write_interaction_dict_in_file(bp_category_dict, os.path.join(directories.annotation_dir, "bp_interaction_stat.txt"))
    print("")
    write_interaction_dict_in_file(stk_category_dict, os.path.join(directories.annotation_dir, "stk_interaction_stat.txt"))
    write_bp_matrix(bp_category_dict, os.path.join(directories.annotation_dir, "bp_interaction_stat_matrix.txt"), bp_list, interaction_list)
    bp_interaction_category_rank = generate_ranking_from_interaction_dict(bp_category_dict)
    stk_interaction_category_rank = generate_ranking_from_interaction_dict(stk_category_dict)

    write_interaction_category_rank_file(bp_category_dict, bp_rank_fn)
    write_interaction_category_rank_file(stk_category_dict, stk_rank_fn)