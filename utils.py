import sys
import os
import math
import copy
import logging
import time
import datetime
import shutil

from logger import *
from classes import *
from interaction_utils import *
# from scoring_utils import *

Infinity = float('inf')
max_score = -Infinity
maximal_clique = []
lower_bound = -Infinity

def get_backbone_and_sugar_atoms():
    
    backbone_atoms = ["C3'", "C4'", "C5'", "O3'", "O5'", "P"]
    sugar_atoms = ["C1'", "C2'", "C3'", "C4'", "O4'"]
    # base_atoms = ["C1'", "C2", "C4", "C5", "C6", "C8", "N1", "N3", "N7", "N9"]

    return backbone_atoms, sugar_atoms
    
def is_indices_across_break_points(break_points, ind1, ind2):
    for i in range(len(break_points)):
        if (ind1 < break_points[i] and ind2 >= break_points[i]) or (ind2 < break_points[i] and ind1 >= break_points[i]):
            return True
    return False

def get_break_point_ind(break_points, ind):
    break_points = sorted(break_points)
    for i in range(len(break_points)):
        if ind < break_points[i]:
            return i
    return -1

def get_three_letter_string_for_bp(interaction, orientation):
    return orientation[0] + ''.join(interaction.split('/'))

def get_reversed_interaction(interaction, interaction_type):
    if interaction_type == 'bp':
        return interaction[0] + interaction[1:][::-1]
    elif interaction_type == 'stack':
        if interaction == 'upward':
            return 'downward'
        elif interaction == 'downward':
            return 'upward'
        elif interaction == 'inward':
            # return 'outward'
            return interaction
        elif interaction == 'outward':
            # return 'inward'
            return interaction
        else:
            logger.error('Invalid stack interaction')
            sys.exit()
    else:
        # logger.error('Invalid interaction type')
        sys.exit()

def load_loop_data(loop_fn, include_reversed=True):
    
    fp = open(loop_fn)
    lines = fp.readlines()
    fp.close()

    

    sequence = lines[1].strip()
    sequence_pieces = sequence.strip().split('...')
    joined_sequence = ''.join(sequence_pieces)

    break_points = []
    prev = 0
    for i in range(len(sequence_pieces)):
        break_points.append(len(sequence_pieces[i]) + prev)
        prev += len(sequence_pieces[i])

    # bps = []
    # stks = []
    bps = {}
    bp_cnt = 0
    stks = {}
    stk_cnt = 0

    reading_bp = True
    for line in lines[3:]:
        line = line.strip()
        pieces = line.split(',')

        if line == '#info=stacking':
            reading_bp = False
            continue

        s, e = pieces[0].strip().split('-')
        s = int(s)
        e = int(e)
        ntd_pair = joined_sequence[s] + joined_sequence[e]
        rev_ntd_pair = joined_sequence[e] + joined_sequence[s]

        if reading_bp == True:
            if s not in bps:
                bps[s] = {}
            if e not in bps[s]:
                bps[s][e] = []
            interaction = get_three_letter_string_for_bp(pieces[1].strip(), pieces[2].strip())
            bps[s][e].append((s, e, ntd_pair, interaction))
            # bps[s][e].append((s, e, hash_basepair(interaction, ntd_pair)))
            if include_reversed == True:
                if e not in bps:
                    bps[e] = {}
                if s not in bps[e]:
                    bps[e][s] = []
                # rev_ntd_pair = ntd_pair[1] + ntd_pair[0]
                rev_interaction = get_reversed_interaction(interaction, 'bp')
                bps[e][s].append((e, s, rev_ntd_pair, rev_interaction))
                # bps[e][s].append((e, s, hash_basepair(rev_interaction, rev_ntd_pair)))
            bp_cnt += 1
        else:
            if s not in stks:
                stks[s] = {}
            if e not in stks[s]:
                stks[s][e] = []
            interaction = pieces[1].strip()
            stks[s][e].append((s, e, ntd_pair, interaction))
            # stks[s][e].append((s, e, hash_stacking(interaction, ntd_pair)))
            if include_reversed == True:
                if e not in stks:
                    stks[e] = {}
                if s not in stks[e]:
                    stks[e][s] = []
                rev_interaction = get_reversed_interaction(interaction, 'stack')
                stks[e][s].append((e, s, rev_ntd_pair, rev_interaction))
                # stks[e][s].append((e, s, hash_stacking(rev_interaction, rev_ntd_pair)))
            stk_cnt += 1

    return joined_sequence, bps, bp_cnt, stks, stk_cnt, break_points

def load_profile_data(profile_fn, profile_name):
    
    sequenceP = []
    bpsP = {}
    stksP = {}
    break_points = []
    loopCount = 0

    fp = open(profile_fn)
    lines = fp.readlines()
    fp.close()

    reading_label = None
    for line in lines:
        if line.startswith('#'):
            reading_label = line.strip().split('=')[1]
            continue

        if reading_label == 'sequence':
            nucl_dict_in_cur_location = {}
            
            nucl_freq_pairs = line.strip().split(', ')
            for nucl_freq_pair in nucl_freq_pairs:
                nucl, freq = nucl_freq_pair.strip().split(' ')

                # nucl, freq = line.strip().split(' ')
                nucl = nucl.strip(':')
                freq = int(freq.strip())

                hashed_nucl = hash_nucleotide(nucl)
            
                if hashed_nucl not in nucl_dict_in_cur_location:
                    nucl_dict_in_cur_location[hashed_nucl] = 0
                nucl_dict_in_cur_location[hashed_nucl] += freq

            sequenceP.append(nucl_dict_in_cur_location)

        elif reading_label == 'basepair':
            location, interaction_with_freq_list = line.strip().split('\t')
            s, e = location.strip().split('-')
            s = int(s)
            e = int(e)
            if (s, e) not in bpsP:
                bpsP[(s, e)] = {}

            interaction_with_freq_list = interaction_with_freq_list.strip().split(', ')
            for interaction_with_freq in interaction_with_freq_list:
                intera_with_ntd_pair, freq = interaction_with_freq.strip().split(':')
                interaction, ntd_pair = intera_with_ntd_pair.strip().split(' - ')
                interaction = interaction.strip()
                ntd_pair = ntd_pair.strip()
                freq = int(freq.strip())
                hashed_interaction = hash_basepair(interaction, ntd_pair)
            
                if hashed_interaction not in bpsP[(s, e)]:
                    bpsP[(s, e)][hashed_interaction] = 0
                bpsP[(s, e)][hashed_interaction] += freq

        elif reading_label == 'stacking':
            location, interaction_with_freq_list = line.strip().split('\t')
            s, e = location.strip().split('-')
            s = int(s)
            e = int(e)
            if (s, e) not in stksP:
                stksP[(s, e)] = {}

            interaction_with_freq_list = interaction_with_freq_list.strip().split(', ')
            for interaction_with_freq in interaction_with_freq_list:    
                interaction, freq = interaction_with_freq.strip().split(':')
                interaction = interaction.strip()
                freq = int(freq.strip())
                hashed_interaction = hash_stacking(interaction, '')
                
                if hashed_interaction not in stksP[(s, e)]:
                    stksP[(s, e)][hashed_interaction] = 0
                stksP[(s, e)][hashed_interaction] += freq

        elif reading_label == 'break_points':
            break_points = list(map(int, line.strip().split(',')))

        elif reading_label == 'total_loops':
            loopCount = int(line.strip())

        else:
            print('Invalid reading label.')
            sys.exit()

    return Profile(sequenceP, bpsP, stksP, break_points, loopCount, profile_name)

def add_interactions_from_src_ind_to_interaction_dict(source, interaction_dict, interaction_profile_dict):
    if source in interaction_dict:
        for target in interaction_dict[source]:
            for item in interaction_dict[source][target]:
                # s, e, ntd_pair, interaction = item
                s, e, hashed_interaction = item
                if (s, e) not in interaction_profile_dict:
                    interaction_profile_dict[(s, e)] = {}

                # if interaction not in interaction_profile_dict[(s, e)]:
                #     interaction_profile_dict[(s, e)][interaction] = []

                if hashed_interaction not in interaction_profile_dict[(s, e)]:
                    interaction_profile_dict[(s, e)][hashed_interaction] = 0

                # interaction_profile_dict[(s, e)][interaction].append(ntd_pair)
                interaction_profile_dict[(s, e)][hashed_interaction] += 1

def MotifToProfile(motif):
    sequenceP = []
    bpsP = {}
    stksP = {}
    for i in range(len(motif.sequence)):
        # nucl_list_in_cur_location = []
        # nucl_list_in_cur_location.append(motif.sequence[i])
        # sequenceP.append(nucl_list_in_cur_location)

        hashed_nucl = hash_nucleotide(motif.sequence[i])
        nucl_dict_in_cur_location = {}
        if hashed_nucl not in nucl_dict_in_cur_location:
            nucl_dict_in_cur_location[hashed_nucl] = 0
        nucl_dict_in_cur_location[hashed_nucl] += 1
        sequenceP.append(nucl_dict_in_cur_location)

        add_interactions_from_src_ind_to_interaction_dict(i, motif.bps, bpsP)
        add_interactions_from_src_ind_to_interaction_dict(i, motif.stks, stksP)

    # print(sequenceP)
    # print(bpsP)
    # print(stksP)
    return Profile(sequenceP, bpsP, stksP, motif.break_points, 1, motif.loopName)

def get_interaction_ind_pair_list(profile):
    interaction_ind_pair_list = []
    for (x, y) in profile.bpsP:
        # hashed_interaction_list = []
        hashed_interaction_list = list(profile.bpsP[(x, y)].keys())
        # for interaction in profile.bpsP[(x, y)]:
        #     for ntd_pair in profile.bpsP[(x, y)][interaction]:
        #         hashed_interaction_list.append(hash_basepair(interaction, ntd_pair))
        interaction_ind_pair_list.append((x, y, 'b', hashed_interaction_list))

    for (x, y) in profile.stksP:
        # hashed_interaction_list = []
        hashed_interaction_list = list(profile.stksP[(x, y)].keys())
        # for interaction in profile.stksP[(x, y)]:
        #     for ntd_pair in profile.stksP[(x, y)][interaction]:
        #         hashed_interaction_list.append(hash_stacking(interaction, ntd_pair))
        interaction_ind_pair_list.append((x, y, 's', hashed_interaction_list))

    # return list(profile.bpsP.keys())
     # + list(profile.stksP.keys())

    return interaction_ind_pair_list

def sort_interaction_ind_pairs(interaction_ind_pair_list):
    n = len(interaction_ind_pair_list)
    for i in range(n - 1):
        for j in range(i+1, n):
            if (interaction_ind_pair_list[i][1] > interaction_ind_pair_list[j][1]) \
            or \
            (interaction_ind_pair_list[i][1] == interaction_ind_pair_list[j][1] \
            and interaction_ind_pair_list[i][0] < interaction_ind_pair_list[j][0]):
                temp = interaction_ind_pair_list[i]
                interaction_ind_pair_list[i] = interaction_ind_pair_list[j]
                interaction_ind_pair_list[j] = temp
    return interaction_ind_pair_list

def compute_interaction_compatibility(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, include_stackings=True):
    l1 = len(interaction_ind_pair_list1)
    l2 = len(interaction_ind_pair_list2)
    rows, cols = (l1 * l2, l1 * l2)
    interaction_compatibility_matrix = [[0 for i in range(cols)] for j in range(rows)]
    # rank_dict = {}
    # for i in range(l1 * l2):
    #     rank_dict[i] = 0
    
    for i in range(l1-1):
        for j in range(l2-1):
            ind1 = i * l2 + j

            if interaction_ind_pair_list1[i][2] != interaction_ind_pair_list2[j][2]:
                #to ensure that a basepair cannot match with a stacking
                continue

            ind_across_break_points1 = is_indices_across_break_points(profile1.break_points, interaction_ind_pair_list1[i][0], interaction_ind_pair_list1[i][1])
            ind_across_break_points2 = is_indices_across_break_points(profile2.break_points, interaction_ind_pair_list2[j][0], interaction_ind_pair_list2[j][1])
            if ind_across_break_points1 != ind_across_break_points2:
                continue

            if ind_across_break_points1 == False:
                b1 = get_break_point_ind(profile1.break_points, interaction_ind_pair_list1[i][0])
                b2 = get_break_point_ind(profile2.break_points, interaction_ind_pair_list2[j][0])
                if b1 != b2:
                    continue

            for k in range(i+1, l1):
                for l in range(j+1, l2):
                    ind2 = k * l2 + l
                    if interaction_ind_pair_list1[k][2] != interaction_ind_pair_list2[l][2]:
                        #to ensure that a basepair cannot match with a stacking
                        continue

                    ind_across_break_points1 = is_indices_across_break_points(profile1.break_points, interaction_ind_pair_list1[k][0], interaction_ind_pair_list1[k][1])
                    ind_across_break_points2 = is_indices_across_break_points(profile2.break_points, interaction_ind_pair_list2[l][0], interaction_ind_pair_list2[l][1])
                    if ind_across_break_points1 != ind_across_break_points2:
                        continue

                    if ind_across_break_points1 == False:
                        b1 = get_break_point_ind(profile1.break_points, interaction_ind_pair_list1[k][0])
                        b2 = get_break_point_ind(profile2.break_points, interaction_ind_pair_list2[l][0])
                        if b1 != b2:
                            continue

                    # if interaction_ind_pair_list1[i][2] != interaction_ind_pair_list1[k][2] or \
                    # interaction_ind_pair_list2[j][2] != interaction_ind_pair_list2[l][2]:
                    #     continue

                    if include_stackings == False:
                        if interaction_ind_pair_list1[i][2] == 's' or interaction_ind_pair_list1[k][2] == 's' or \
                        interaction_ind_pair_list2[j][2] == 's' or interaction_ind_pair_list2[l][2] == 's':
                            continue

                    if (interaction_ind_pair_list1[k][1] < interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list2[l][1] < interaction_ind_pair_list2[j][0]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][0] > interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list2[l][0] > interaction_ind_pair_list2[j][1]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][1] == interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list2[l][1] == interaction_ind_pair_list2[j][0]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][0] == interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list2[l][0] == interaction_ind_pair_list2[j][1]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][0] == interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list1[k][1] < interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list2[l][0] == interaction_ind_pair_list2[j][0] \
                    and interaction_ind_pair_list2[l][1] < interaction_ind_pair_list2[j][1]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][0] == interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list1[k][1] > interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list2[l][0] == interaction_ind_pair_list2[j][0] \
                    and interaction_ind_pair_list2[l][1] > interaction_ind_pair_list2[j][1]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][1] == interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list1[k][0] < interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list2[l][1] == interaction_ind_pair_list2[j][1] \
                    and interaction_ind_pair_list2[l][0] < interaction_ind_pair_list2[j][0]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][1] == interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list1[k][0] > interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list2[l][1] == interaction_ind_pair_list2[j][1] \
                    and interaction_ind_pair_list2[l][0] > interaction_ind_pair_list2[j][0]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][0] > interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list1[k][1] < interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list2[l][0] > interaction_ind_pair_list2[j][0] \
                    and interaction_ind_pair_list2[l][1] < interaction_ind_pair_list2[j][1]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][0] < interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list1[k][1] > interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list2[l][0] < interaction_ind_pair_list2[j][0] \
                    and interaction_ind_pair_list2[l][1] > interaction_ind_pair_list2[j][1]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][0] < interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list1[k][1] > interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list1[k][1] < interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list2[l][0] < interaction_ind_pair_list2[j][0] \
                    and interaction_ind_pair_list2[l][1] > interaction_ind_pair_list2[j][0] \
                    and interaction_ind_pair_list2[l][1] < interaction_ind_pair_list2[j][1]) \
                    \
                    or \
                    \
                    (interaction_ind_pair_list1[k][0] > interaction_ind_pair_list1[i][0] \
                    and interaction_ind_pair_list1[k][0] < interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list1[k][1] > interaction_ind_pair_list1[i][1] \
                    and interaction_ind_pair_list2[l][0] > interaction_ind_pair_list2[j][0] \
                    and interaction_ind_pair_list2[l][0] < interaction_ind_pair_list2[j][1] \
                    and interaction_ind_pair_list2[l][1] > interaction_ind_pair_list2[j][1]):

                        interaction_compatibility_matrix[ind1][ind2] = 1
                        interaction_compatibility_matrix[ind2][ind1] = 1
                        # print('#######################################yes##################################')
                        # sys.exit()
                        # rank_dict[ind1] += 1
                        # rank_dict[ind2] += 1
                
    # return interaction_compatibility_matrix, rank_dict
    return interaction_compatibility_matrix

def is_acceptable_interaction_for_triple_bonus(interaction_ind_pair):
    # if interaction_ind_pair[2] == 'b' and interpret_hash_value(interaction_ind_pair[3].keys()[0]) != 'cWW':
    if interaction_ind_pair[2] == 'b' and interpret_hash_value(max(set(interaction_ind_pair[3]), key = interaction_ind_pair[3].count))[0] != 'cWW':
        return True
    return False

def significant_match(interaction_ind_pair1, interaction_ind_pair2):
    # if len(interaction_ind_pair1[3]) == 1 and len(interaction_ind_pair2[3]) == 1 and \
    if is_acceptable_interaction_for_triple_bonus(interaction_ind_pair1) and \
    is_acceptable_interaction_for_triple_bonus(interaction_ind_pair2) and \
    max(set(interaction_ind_pair1[3]), key = interaction_ind_pair1[3].count) == max(set(interaction_ind_pair2[3]), key = interaction_ind_pair2[3].count):
        # interaction_ind_pair1[3].keys()[0] == interaction_ind_pair2[3].keys()[0]:
        #  similar type of non-canonical base pairs, same edges orientation
        return True
    else:  
        return False

# def add_a_region(aligned_regions, start, end):
#     aligned_regions.append((start, end))

# def extend_a_regions(aligned_regions, start, end):
#     last_ind = aligned_regions[-1][1]
#     if last_ind + 1 != start:
#         print("Error")
#     aligned_regions[-1][1] = end

# def get_aligned_regions(aligned_nucl_indices1, aligned_regions1, aligned_nucl_indices2, aligned_regions2):
def get_aligned_regions(aligned_nucl_indices1, aligned_nucl_indices2):
    aligned_regions1 = []
    aligned_regions2 = []

    if len(aligned_nucl_indices1) < 1 or len(aligned_nucl_indices2) < 1:
        return

    aligned_nucl_indices1.sort()
    aligned_nucl_indices2.sort()

    #add_a_region()
    add_a_region(aligned_regions1, aligned_nucl_indices1[0], aligned_nucl_indices1[0])
    add_a_region(aligned_regions2, aligned_nucl_indices2[0], aligned_nucl_indices2[0])
    # aligned_regions1.append((aligned_nucs1[0], aligned_nucs1[0]))
    # aligned_regions2.append((aligned_nucs2[0], aligned_nucs2[0]))

    for i in range(1, len(aligned_nucl_indices1)):
        if aligned_nucl_indices1[i] > aligned_nucl_indices1[i-1] or aligned_nucl_indices2[i] > aligned_nucl_indices2[i-1]:
            if aligned_nucl_indices1[i] > aligned_nucl_indices1[i - 1]:
                # print('from utils')
                extend_a_region(aligned_regions1, aligned_nucl_indices1[i - 1] + 1, aligned_nucl_indices1[i])
                # aligned_regions1[-1][1] = aligned_nucl_indices1[i]
            if aligned_nucl_indices2[i] > aligned_nucl_indices2[i - 1]:
                # print('from utils')
                extend_a_region(aligned_regions2, aligned_nucl_indices2[i - 1] + 1, aligned_nucl_indices2[i])
                # aligned_regions2[-1][1] = aligned_nucl_indices2[i]

    return aligned_regions1, aligned_regions2

def get_aligned_nucs(clique, interaction_ind_pair_list1, interaction_ind_pair_list2):
    # const std::vector<int>& clique,
    # std::vector<int>& nucs_first,
    # std::vector<int>& nucs_second)

    interaction_count1 = len(interaction_ind_pair_list1)
    interaction_count2 = len(interaction_ind_pair_list2)
    # print('interaction_count')
    # print(interaction_count1, interaction_count2)
    clique_size = len(clique)

    #cout << "get_aligned_nucs:  clique size:  " << clique.size() << endl;
    #  construct the nucleotides that particitate in selected interactions
    nucleotide_indices1 = []
    nucleotide_indices2 = []
    for i in range(clique_size):
        #  the two base interactions idx_i and idx_j match
        ind_i = clique[i] // interaction_count2
        ind_j = clique[i] % interaction_count2
        #  copy the matched indices

        # print(clique[i])
        # print(ind_i, ind_j)

        nucleotide_indices1.append(interaction_ind_pair_list1[ind_i][0])
        nucleotide_indices1.append(interaction_ind_pair_list1[ind_i][1])
        nucleotide_indices2.append(interaction_ind_pair_list2[ind_j][0])
        nucleotide_indices2.append(interaction_ind_pair_list2[ind_j][1])

        # print(nucleotide_indices1)
        # print(nucleotide_indices2)

    #  assert that the sizes of the matched nucleotides are the same
    if len(nucleotide_indices1) != len(nucleotide_indices2):
        print('The sizes of matched nucleotides do not match: MotifGraphMatching::compute_alignment_score')
        return [], []

    if len(nucleotide_indices1) == 0:
        return [], []

    #  sort the matched nucleotides
    nucleotide_indices1.sort()
    nucleotide_indices2.sort()

    return nucleotide_indices1, nucleotide_indices2

def is_better_alignment(current_score, max_score):
    return current_score >= max_score

def clique_based_extended_sequence_alignment(estimated_clique, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, aligned_regions1, aligned_regions2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2):

    interaction_count1 = len(interaction_ind_pair_list1)
    interaction_count2 = len(interaction_ind_pair_list2)
    estimated_clique_size = len(estimated_clique)

    pairs_record1 = [0 for i in range(interaction_count1)]
    pairs_record2 = [0 for i in range(interaction_count2)]

    #  flag the base pairs that can be used to identify the missing BP for adding penalty
    for i in range(len(estimated_clique)):
        ind_i = estimated_clique[i] // interaction_count2
        ind_j = estimated_clique[i] % interaction_count2

        pairs_record1[ind_i] = 1
        pairs_record2[ind_j] = 1

    
    aligned_nucleotide_indices1, aligned_nucleotide_indices2 = get_aligned_nucs(estimated_clique, interaction_ind_pair_list1, interaction_ind_pair_list2)

    evaluate_score = 0
    # aligned_seq1 = "";
    # aligned_seq2 = "";
    if len(aligned_nucleotide_indices1) == 0:
        # Best sequence alignment, considering adjacent stack match and missing BP/Stack penalty
        # evaluate_score += local_loop_align(aligned_seq1, aligned_seq2, aligned_nucleotide_indices1, aligned_nucleotide_indices2);
        evaluate_score += local_loop_align(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2)
        # evaluate_score += loop_align_sequence_only(0, 0, hold_alignment_1, hold_alignment_2, M1_seq[seq_index].length() - 1, M2_seq[seq_index].length() - 1);
    print(evaluate_score)
    sys.exit()
    # else:
    
    #     # sort the matched nucleotides
    #     # sort(M1_nuc.begin(), M1_nuc.end());
    #     # sort(M2_nuc.begin(), M2_nuc.end());
    #     # vector<int> M1_nuc_extended(M1_nuc);
    #     # vector<int> M2_nuc_extended(M2_nuc);

    #     aligned_nucleotide_indices1_extended = copy.deepcopy(aligned_nucleotide_indices1)
    #     aligned_nucleotide_indices2_extended = copy.deepcopy(aligned_nucleotide_indices2)

    #     # get the region defined by the clique
    #     get_aligned_regions(aligned_nucleotide_indices1, aligned_regions1, aligned_nucleotide_indices2, aligned_regions2)
    #     #cout << "After calling get_aligned_regions" << endl;
    #     # If global alignment, add start and end of the sequence
    #     # if(!is_local_align)
    #     # {
    #     string M1_hold_alignment;
    #     string M2_hold_alignment;
    #     evaluate_score += loop_align_beginning(aligned_nucleotide_indices1, aligned_nucleotide_indices2, M1_hold_alignment, M2_hold_alignment);
    #     aligned_seq1 += M1_hold_alignment;
    #     aligned_seq2 += M2_hold_alignment;
    #     # }
    #     # else
    #     # {
    #     #   if(EXTEND_LOOP)
    #     #   {
    #     #     int prefix_score = 0;
    #     #     string M1_hold_alignment_prefix;
    #     #     string M2_hold_alignment_prefix;
    #     #     vector<int> M1_nuc_prefix_match, M2_nuc_prefix_match;
    #     #     # Extend to the left of the segment
    #     #     prefix_score = end_fixed_beginning_open_loop_align(0, 0, M1_hold_alignment_prefix, M2_hold_alignment_prefix, M1_nuc[0], M2_nuc[0], M1_regions, M2_regions, M1_pairs_record, M2_pairs_record, M1_nuc_prefix_match, M2_nuc_prefix_match);
    #     #     extend_segments(M1_nuc_extended, M1_nuc_prefix_match, M1_regions, M2_nuc_extended, M2_nuc_prefix_match, M2_regions);
    #     #     evaluate_score += prefix_score;
    #     #     aligned_seq1 += M1_hold_alignment_prefix;
    #     #     aligned_seq2 += M2_hold_alignment_prefix;

    #     #     # cout << "check1: aligning " << M1_hold_alignment_prefix << " and " << M2_hold_alignment_prefix << endl;
    #     #   }
    #     # }
    #     # int i;
    #     # utility_functions::print_formatted_int_vector(M1_nuc);
    #     # utility_functions::print_formatted_int_vector(M2_nuc);
    #     # cout << endl;

    #     # for(i = 0; i < (int) M1_nuc.size() - 1; ++ i)  {
    #     for i in range(len(aligned_nucleotide_indices1)-1):
    #         if M1_nuc[i] != M1_nuc[i + 1] || M2_nuc[i] != M2_nuc[i + 1]:
    #             if(M1_nuc[i] != M1_nuc[i + 1] && M2_nuc[i] != M2_nuc[i + 1])
    #             {
    #               aligned_seq1.append(1, M1_seq[seq_index][M1_nuc[i]]);
    #               aligned_seq2.append(1, M2_seq[seq_index][M2_nuc[i]]);
    #             }

    #             bool is_segment1_across_breakpoint = is_indices_across_break_points(M1_break_points[0], M1_nuc[i], M1_nuc[i + 1]);
    #             bool is_segment2_across_breakpoint = is_indices_across_break_points(M2_break_points[0], M2_nuc[i], M2_nuc[i + 1]);
    #             # Segments not across breakpoint
    #             if(!is_segment1_across_breakpoint && !is_segment2_across_breakpoint)
    #             {
    #               string M1_hold_alignment;
    #               string M2_hold_alignment;
    #               # evaluate the global sequence alignment score of the segment
    #               evaluate_score += loop_align(M1_nuc[i], M2_nuc[i], M1_hold_alignment, M2_hold_alignment, M1_nuc[i + 1], M2_nuc[i + 1]);
    #               aligned_seq1 += M1_hold_alignment;
    #               aligned_seq2 += M2_hold_alignment;
    #             }
    #             else 
    #             {
    #               if(EXTEND_LOOP)
    #               {
    #               string M1_hold_alignment_inside1, M1_hold_alignment_inside2;
    #               string M2_hold_alignment_inside1, M2_hold_alignment_inside2;
    #               vector<int> M1_nuc_inside_match1, M2_nuc_inside_match1, M1_nuc_inside_match2, M2_nuc_inside_match2;

    #               # Extend to the right of the segment
    #               int score1 = beginning_fixed_end_open_loop_align(M1_nuc[i], M2_nuc[i], M1_hold_alignment_inside1, M2_hold_alignment_inside1, M1_nuc[i + 1] - 1, M2_nuc[i + 1] - 1, M1_regions, M2_regions, M1_pairs_record, M2_pairs_record, M1_nuc_inside_match1, M2_nuc_inside_match1);

    #               extend_segments(M1_nuc_extended, M1_nuc_inside_match1, M1_regions, M2_nuc_extended, M2_nuc_inside_match1, M2_regions);
    #               aligned_seq1 += M1_hold_alignment_inside1;
    #               aligned_seq2 += M2_hold_alignment_inside1;
    #               evaluate_score += score1;

    #               # cout << "check2: aligning " << M1_hold_alignment_inside1 << " and " << M2_hold_alignment_inside1 << endl;

    #               int M1_inside2_start = M1_nuc[i] + 1;
    #               int M2_inside2_start = M2_nuc[i] + 1;

    #               #If inside1 is extended, ignore that for the next alignment
    #               if(M1_nuc_inside_match1.size() > 0)
    #               {
    #                 sort(M1_nuc_inside_match1.begin(), M1_nuc_inside_match1.end());
    #                 sort(M2_nuc_inside_match1.begin(), M2_nuc_inside_match1.end());
    #                 M1_inside2_start = M1_nuc_inside_match1.back() + 1;
    #                 M2_inside2_start = M2_nuc_inside_match1.back() + 1;
    #               }

    #               # Extend to the left of the segment
    #               int score2 = end_fixed_beginning_open_loop_align(M1_inside2_start, M2_inside2_start, M1_hold_alignment_inside2, M2_hold_alignment_inside2, M1_nuc[i + 1], M2_nuc[i + 1], M1_regions, M2_regions, M1_pairs_record, M2_pairs_record, M1_nuc_inside_match2, M2_nuc_inside_match2);
    #               extend_segments(M1_nuc_extended, M1_nuc_inside_match2, M1_regions, M2_nuc_extended, M2_nuc_inside_match2, M2_regions);
    #               aligned_seq1 += M1_hold_alignment_inside2;
    #               aligned_seq2 += M2_hold_alignment_inside2;
    #               evaluate_score += score2;

    #               # cout << "check3: aligning " << M1_hold_alignment_inside2 << " and " << M2_hold_alignment_inside2 << endl;
    #             }
    #             if(!is_local_align)
    #             {
    #               string M1_hold_alignment;
    #               string M2_hold_alignment;
    #               # evaluate the global sequence alignment score of the segment
    #               # evaluate_score += loop_align(M1_nuc[i], M2_nuc[i], M1_hold_alignment, M2_hold_alignment, M1_nuc[i + 1], M2_nuc[i + 1]);
    #               evaluate_score += loop_align_sequence_only(M1_nuc[i]+1, M2_nuc[i]+1, M1_hold_alignment, M2_hold_alignment, M1_nuc[i + 1]-1, M2_nuc[i + 1]-1);
    #               aligned_seq1 += M1_hold_alignment;
    #               aligned_seq2 += M2_hold_alignment;

    #               # string M1_hold_alignment;
    #               # string M2_hold_alignment;
    #               # loop_align_sequence_only(
    #               #   M1_nuc[i] + 1, M2_nuc[i] + 1,
    #               #   M1_hold_alignment, M2_hold_alignment,
    #               #   M1_nuc[i + 1] - 1, M2_nuc[i + 1] - 1
    #               # );
    #               # cout << "Before adding hold alignment" << endl;
    #               # cout << aligned_seq1 << endl;
    #               # cout << aligned_seq2 << endl;
    #               # aligned_seq1 += M1_hold_alignment;
    #               # aligned_seq2 += M2_hold_alignment;
    #               # cout << "After adding hold alignment" << endl;
    #               # cout << aligned_seq1 << endl;
    #               # cout << aligned_seq2 << endl;
    #             }
    #             }
        

    #     aligned_seq1.append(1, M1_seq[seq_index][M1_nuc[i]]);
    #     aligned_seq2.append(1, M2_seq[seq_index][M2_nuc[i]]);


    #     # if(EXTEND_LOOP)
    #     # {
    #     #   int suffix_score = 0;
    #     #   string M1_hold_alignment_suffix;
    #     #   string M2_hold_alignment_suffix;
    #     #   vector<int> M1_nuc_suffix_match, M2_nuc_suffix_match;

    #     #   # Extend to the right of the segment
    #     #   suffix_score = beginning_fixed_end_open_loop_align(M1_nuc[M1_nuc.size() - 1], M2_nuc[M2_nuc.size() - 1], M1_hold_alignment_suffix, M2_hold_alignment_suffix, M1_seq[seq_index].length() - 1, M2_seq[seq_index].length() - 1, M1_regions, M2_regions, M1_pairs_record, M2_pairs_record, M1_nuc_suffix_match, M2_nuc_suffix_match);
    #     #   extend_segments(M1_nuc_extended, M1_nuc_suffix_match, M1_regions, M2_nuc_extended, M2_nuc_suffix_match, M2_regions);
    #     #   aligned_seq1 += M1_hold_alignment_suffix;
    #     #   aligned_seq2 += M2_hold_alignment_suffix;
    #     #   evaluate_score += suffix_score;

    #     #   # cout << "check4: aligning " << M1_hold_alignment_suffix << " and " << M2_hold_alignment_suffix << endl;
    #     # }
    #     # if(!is_local_align)
    #     # {
    #     string M1_hold_alignment;
    #     string M2_hold_alignment;
    #     evaluate_score += loop_align_end(M1_nuc, M2_nuc, M1_hold_alignment, M2_hold_alignment);
    #     aligned_seq1 += M1_hold_alignment;
    #     aligned_seq2 += M2_hold_alignment;
    #     # }

    #     # place_break_point_indicators(aligned_seq1, aligned_seq2, M1_break_points[0], M2_break_points[0]);
    
    return evaluate_score


def compute_alignment_score(estimated_clique, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2):

    interaction_count1 = len(interaction_ind_pair_list1)
    interaction_count2 = len(interaction_ind_pair_list2)
    estimated_clique_size = len(estimated_clique)

    pairs_record1 = [0 for i in range(interaction_count1)]
    pairs_record2 = [0 for i in range(interaction_count2)]

    matched_pairs = []
    nucleotide_indices1 = []
    nucleotide_indices2 = []
    matched_positions1 = [0 for i in range(len(profile1.sequenceP))]
    matched_positions2 = [0 for i in range(len(profile2.sequenceP))]

    evaluate_score = 0.0

    for i in range(estimated_clique_size):
        ind_i = estimated_clique[i] // interaction_count2
        ind_j = estimated_clique[i] % interaction_count2
        pairs_record1[ind_i] = 1
        pairs_record2[ind_j] = 1

        # if interaction_ind_pair_list1[ind_i][2] == 'b':
        matched_pairs.append((ind_i, ind_j))

        if (interaction_ind_pair_list1[ind_i][2] == 'b' and interaction_ind_pair_list2[ind_j][2] == 'b'):
            evaluate_score += scores.weight_isosteric * match_isosteric_basepair_profile(interaction_ind_pair_list1[ind_i][3], interaction_ind_pair_list2[ind_j][3])
        #  match of two non-adjacent stackings. hbond may present
        elif (interaction_ind_pair_list1[ind_i][2] == 's' and interaction_ind_pair_list2[ind_j][2] == 's'):
            evaluate_score += scores.weight_nonadjacent_stacking * match_stacking_profile(interaction_ind_pair_list1[ind_i][3], interaction_ind_pair_list2[ind_j][3])
        # match of two hydrogen bonds
        # elif (M1_interaction_union[idx_i][2] == 0 and M2_interaction_union[idx_j][2] == 0):
        #     evaluate_score += scores.weight_isosteric * (scores.hbond_match_base + scores.bonuses.hbond_match_bonus)

        # Score for the matching of nucleotide in the basepairs
        evaluate_score += match_nucleotide_profile(profile1.sequenceP[interaction_ind_pair_list1[ind_i][0]], profile2.sequenceP[interaction_ind_pair_list2[ind_j][0]])
        evaluate_score += match_nucleotide_profile(profile1.sequenceP[interaction_ind_pair_list1[ind_i][1]], profile2.sequenceP[interaction_ind_pair_list2[ind_j][1]])
        
        nucleotide_indices1.append(interaction_ind_pair_list1[ind_i][0])
        nucleotide_indices1.append(interaction_ind_pair_list1[ind_i][1])
        nucleotide_indices2.append(interaction_ind_pair_list2[ind_j][0])
        nucleotide_indices2.append(interaction_ind_pair_list2[ind_j][1])


        if significant_match(interaction_ind_pair_list1[ind_i], interaction_ind_pair_list2[ind_j]):
            #  check for bonus score of triple interaction
            if matched_positions1[interaction_ind_pair_list1[ind_i][0]] > 0 and matched_positions2[interaction_ind_pair_list2[ind_j][0]] > 0:
                evaluate_score += scores.weight_isosteric * scores.bonuses.triple_interaction_bonus

            if matched_positions1[interaction_ind_pair_list1[ind_i][1]] > 0 and matched_positions2[interaction_ind_pair_list2[ind_j][1]] > 0:
                evaluate_score += scores.weight_isosteric * scores.bonuses.triple_interaction_bonus

            matched_positions1[interaction_ind_pair_list1[ind_i][0]] = 1
            matched_positions1[interaction_ind_pair_list1[ind_i][1]] = 1

            matched_positions2[interaction_ind_pair_list2[ind_j][0]] = 1
            matched_positions2[interaction_ind_pair_list2[ind_j][1]] = 1

    # Asymmetric Penalty (for segment size difference between bases of pairs)
    evaluate_score += asymmetric_penalty(matched_pairs, interaction_ind_pair_list1, interaction_ind_pair_list2, scores)

    # print(evaluate_score)
    # sys.exit()

    if len(nucleotide_indices1) != len(nucleotide_indices2):
        return 0
    else:
        # If empty clique
        if len(nucleotide_indices1) != 0:
            # Apply penalty for unmatched base pairs in the region defined by the clique
            # Note: further missing BP penalty will be given below for the extended regions of loops
            aligned_regions1 = []
            aligned_regions2 = []
            # list<AlignedRegionType> M1_regions, M2_regions;

            get_aligned_regions(nucleotide_indices1, aligned_regions1, nucleotide_indices2, aligned_regions2)
            evaluate_score += get_missing_bp_and_stack_penalty(interaction_ind_pair_list1, interaction_ind_pair_list2, aligned_regions1, aligned_regions2, pairs_record1, pairs_record2, scores)

        # print('here')
        # print(evaluate_score)
        # sys.exit()
        aligned_regions1 = []
        aligned_regions2 = []
        # list<AlignedRegionType> M1_regions, M2_regions;
        # # list<AlignedRegionType> M1_align_nucs, M2_align_nucs;
        # string aligned_seq1, aligned_seq2;
        evaluate_score += clique_based_extended_sequence_alignment(estimated_clique, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, aligned_regions1, aligned_regions2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2)

    return evaluate_score

#Return 1 if successful, 0 otherwise
def extend_clique_recursive(current_clique, candidates, iteration, depth, branch, begin_time, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2):
  # (*iteration)++;
    alignment_score = compute_alignment_score(current_clique, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2)
    if alignment_score > max_score:
        max_score = alignment_score
        maximal_clique = current_clique
    
#   if max_score > lower_bound:
#     lower_bound = max_score
  
#   if candidates.empty():
#     return 1;

#   if(use_branch_and_bound)
#   {
#     int alignment_score_upper_bound = compute_alignment_score_upper_bound(current_clique, candidates);

#     if(alignment_score_upper_bound <= lower_bound)  {
#       #cout << "return condition 2" << endl;
#       return 1;
#     }
#   }

#   int i, j;
#   # cout << "Here 4" << endl;
#   for(i = 0; i < (int) candidates.size(); ++ i)  {
#     vector<int> next_clique = current_clique;
#     next_clique.push_back(candidates[i]);
#     #  pruning candidates
#     vector<int> next_candidates;
#     for(j = 0; j < (int) candidates.size(); ++ j)  {
#       if(j != i && interaction_compatibility[candidates[i]][candidates[j]] == 1)  {
#         next_candidates.push_back(candidates[j]);
#       }
#     }

#     # cout << "Here 5" << endl;

#     if(depth == 0)
#         branch = i;

#     if(branchDepth[branch] < depth)
#         branchDepth[branch] = depth;

#     #3 Cases
#     #Don't stop if stop_time_consuming_scanx == 0. Otherwise
#     #Stop if *iteration >= ITERATION_LIMIT and stop_time_consuming_scanx_early == 1
#     #Stop if *iteration >= ITERATION_LIMIT_MAX
#     if(stop_time_consuming_scanx == 1 && ((stop_time_consuming_scanx_early == 1 && *iteration >= ITERATION_LIMIT) || (*iteration >= ITERATION_LIMIT_MAX)))
#     {
#         if(*iteration >= ITERATION_LIMIT_MAX)
#             cout << endl << "!!! MAX ITERATION LIMIT CROSSED !!!" << endl << endl;
#         return 0;
#     }

#     if(depth == 0)
#     {
#         if(test_code && *iteration - prevIteration > 10)
#         {
#             double current_time = float( clock () - begin_time ) / CLOCKS_PER_SEC;
#             cout << "Time " << current_time << " sec" << "\t";
#             cout << "Step " << i << "\tIteration " << *iteration << "\tdiff " << *iteration - prevIteration << endl;
#         }
#         prevIteration = *iteration;
#     }
#     if(extend_clique_recursive(next_clique, next_candidates, iteration, depth + 1, branch, begin_time) == 0)
#         return 0;
#   }
#   return 1;
# }

def find_maximal_clique_optimal(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2):
    l1 = len(interaction_ind_pair_list1)
    l2 = len(interaction_ind_pair_list2)

    initial_candidates = []
    for i in range(l1):
        for j in range(l2):
            if interaction_ind_pair_list1[i][2] == interaction_ind_pair_list2[j][2]:
                initial_candidates.append(i * l2 + j)

    print(initial_candidates)
    # print(scores.weight_sequence)

    # aGA = GA()

    # print(aGA)
    # degrees = aGA.get_degrees(initial_candidates, interaction_compatibility_matrix);
    # print(degrees)

    # n = len(initial_candidates)
    # branchDepth = [0 for i in range(n)]

    # totalDegree = 0
    # log_product = 0.0
    # for i in range(n):
    #     totalDegree += degree[i]
    #     if degree[i] > 0:
    #         log_product += log2(degree[i])
    
    # if has_more_edges_than_acceptable(totalDegree):
    #     return 0

    # const clock_t begin_time = std::clock();
    begin_time = ''
    iteration = 0;
    initial_clique = []
    success = extend_clique_recursive(initial_clique, initial_candidates, iteration, 0, 0, begin_time, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2);


    return success;




#############################################################################
# Old code end                                                              #
# New code start                                                            #
#############################################################################


def interactions_are_compatible_based_on_breakpoints(interaction_ind_pair1, interaction_ind_pair2, profile1, profile2):
    break_point_ind_s1 = get_break_point_ind(profile1.break_points, interaction_ind_pair1[0])
    break_point_ind_e1 = get_break_point_ind(profile1.break_points, interaction_ind_pair1[1])

    break_point_ind_s2 = get_break_point_ind(profile2.break_points, interaction_ind_pair2[0])
    break_point_ind_e2 = get_break_point_ind(profile2.break_points, interaction_ind_pair2[1])

    if break_point_ind_s1 == break_point_ind_s2 and break_point_ind_e1 == break_point_ind_e2:
        return True
    return False

clique_score_container = {}

# def find_maximal_optimal_clique(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, rank_dict, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG):
def find_maximal_optimal_clique(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG):
    l1 = len(interaction_ind_pair_list1)
    l2 = len(interaction_ind_pair_list2)

    # print(l1, l2)
    # sys.exit()
    initial_candidates = []
    for i in range(l1):
        for j in range(l2):
            # len1 = len(interaction_ind_pair_list1[i][3])
            # len2 = len(interaction_ind_pair_list2[j][3])
            # mn = min(len1, len2)
            # mx = max(len1, len2)
            # if mn != mx:
            #     pc = mn * 100.0 / mx
            #     if pc <= 35.0:
            #         continue
            # if interaction_compatibility_matrix[candidates[i]][candidates[j]] == 1:
            # if interaction_ind_pair_list1[i][2] == 'b' and interaction_ind_pair_list1[i][2] == interaction_ind_pair_list2[j][2]:
            # print(interaction_ind_pair_list1[i][3], interaction_ind_pair_list2[j][3])
            # if interaction_ind_pair_list1[i][2] == interaction_ind_pair_list2[j][2]:
                # if interaction_ind_pair_list1[i][2] == 's' and len1 != len2:
                #     continue
                # print(interaction_ind_pair_list1[i][3], interaction_ind_pair_list2[j][3])
                # if get_break_point_ind(profile1.break_points, i) != get_break_point_ind(profile2.break_points, j):
                    # print('putting indices from different break points together as candidates, hault!')
                    # sys.exit()
                # if get_break_point_ind(profile1.break_points, i) == get_break_point_ind(profile2.break_points, j):
                # initial_candidates.append(i * l2 + j)

            if not interactions_are_compatible_based_on_breakpoints(interaction_ind_pair_list1[i], interaction_ind_pair_list2[j], profile1, profile2):
                continue
            if include_stackings:
                if interaction_ind_pair_list1[i][2] == interaction_ind_pair_list2[j][2]:
                    initial_candidates.append(i * l2 + j)
            else:
                if interaction_ind_pair_list1[i][2] == 'b' and interaction_ind_pair_list1[i][2] == interaction_ind_pair_list2[j][2]:
                    initial_candidates.append(i * l2 + j)

    # print('initial_candidates')
    # print(initial_candidates)
    global clique_score_container
    clique_score_container = {}

    global max_score
    max_score = -Infinity

    global maximal_clique
    maximal_clique = []

    global lower_bound
    lower_bound = -Infinity

    initial_clique = []
    # success = recursive_extend_clique(initial_clique, initial_candidates, 0, 0, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, interaction_compatibility_matrix);
    # recursive_extend_clique(initial_clique, initial_candidates, 0, 0, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, interaction_compatibility_matrix, rank_dict, initialized_basepair_index, use_branch_and_bound, dBG);
    recursive_extend_clique(initial_clique, initial_candidates, 0, 0, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, interaction_compatibility_matrix, initialized_basepair_index, use_branch_and_bound, dBG);
    if dBG:
        logger.debug('maximal_clique: ' + str(maximal_clique))
        # print('maximal_clique')
        # print(maximal_clique)
        logger.debug('max_score: ' + str(max_score))
        # print('max_score')
        # print(max_score)

    return max_score, maximal_clique


# def recursive_extend_clique(current_clique, candidates, depth, branch, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, interaction_compatibility_matrix, rank_dict, initialized_basepair_index, branch_and_bound, dBG):
def recursive_extend_clique(current_clique, candidates, depth, branch, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, interaction_compatibility_matrix, initialized_basepair_index, branch_and_bound, dBG):
    global max_score
    global maximal_clique
    global lower_bound
    global clique_score_container

    # print('current_clique')
    # print(current_clique)
    # print('candidates')
    # print(candidates)
    # (*iteration)++;

    # if len(candidates) > 120:
    #     print(candidates)
    #     sys.exit()
    sorted_current_clique = copy.deepcopy(current_clique)
    sorted_current_clique = tuple(sorted(sorted_current_clique))

    if sorted_current_clique in clique_score_container:
        alignment_score = clique_score_container[sorted_current_clique]
        # print('Found in dict, saving time.')
    else:
        alignment_score = calculate_alignment_score(current_clique, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, initialized_basepair_index, dBG)
        clique_score_container[sorted_current_clique] = alignment_score

    # print('current alignment_score')
    # print(alignment_score)
    if alignment_score > max_score:
        max_score = alignment_score
        maximal_clique = copy.deepcopy(current_clique)


    ########################################################
    # For Branch-and-bound
    ########################################################
    # branch_and_bound = True

    if branch_and_bound == True:
        if max_score > lower_bound:
            lower_bound = max_score

        if len(candidates) == 0:
            return

        # alignment_score_upper_bound = calculate_alignment_score_upper_bound(current_clique, candidates, profile1, profile2, interaction_compatibility_matrix, rank_dict, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, initialized_basepair_index, dBG)
        alignment_score_upper_bound = calculate_alignment_score_upper_bound(current_clique, candidates, profile1, profile2, interaction_compatibility_matrix, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, initialized_basepair_index, dBG)
        if alignment_score_upper_bound <= lower_bound:
            return

    ########################################################
    # For Branch-and-bound end
    ########################################################

    # if len(candidates) == 0:
    #     # return 1
    #     return

    # avg_rank = sum(rank_dict.values()) / len(rank_dict)
    for i in range(len(candidates)):
        next_clique = copy.deepcopy(current_clique)
        next_clique.append(candidates[i])
        next_candidates = []
        for j in range(len(candidates)):
            # if j != i and interaction_compatibility_matrix[candidates[i]][candidates[j]] == 1 and (rank_dict[candidates[i]] >= avg_rank and rank_dict[candidates[j]] >= avg_rank):
            if j != i and interaction_compatibility_matrix[candidates[i]][candidates[j]] == 1:
                next_candidates.append(candidates[j])
            # else:
            #     print('i=j or not compatible')
            #     print(i, j)
            #     print(candidates[i], candidates[j])

        # recursive_extend_clique(next_clique, next_candidates, depth, branch, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, interaction_compatibility_matrix, rank_dict, initialized_basepair_index, branch_and_bound, dBG)
        recursive_extend_clique(next_clique, next_candidates, depth, branch, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, interaction_compatibility_matrix, initialized_basepair_index, branch_and_bound, dBG)
        # if recursive_extend_clique(next_clique, next_candidates, depth, branch, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, interaction_compatibility_matrix, initialized_basepair_index) == 0:
        #     return 0

    # return 1

def calculate_alignment_score(estimated_clique, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, basepair_index, dBG):
    evaluate_score = 0

    interaction_count1 = len(interaction_ind_pair_list1)
    interaction_count2 = len(interaction_ind_pair_list2)

    matched_pairs = []
    nucleotide_indices1 = []
    nucleotide_indices2 = []
    matched_positions1 = [0 for i in range(len(profile1.sequenceP))]
    matched_positions2 = [0 for i in range(len(profile2.sequenceP))]

    for i in range(len(estimated_clique)):
        ind_i = estimated_clique[i] // interaction_count2
        ind_j = estimated_clique[i] % interaction_count2

        # print(ind_i, ind_j)
        matched_pairs.append((ind_i, ind_j))

        if (interaction_ind_pair_list1[ind_i][2] == 'b' and interaction_ind_pair_list2[ind_j][2] == 'b'):
            evaluate_score += scores.weight_isosteric * match_isosteric_basepair_profile(interaction_ind_pair_list1[ind_i][3], interaction_ind_pair_list2[ind_j][3], basepair_index, scores, bp_scoring_data)
            # print('\tmatching bp')
            # print(evaluate_score)
        #  match of two non-adjacent stackings. hbond may present
        elif (interaction_ind_pair_list1[ind_i][2] == 's' and interaction_ind_pair_list2[ind_j][2] == 's'):
            evaluate_score += scores.weight_nonadjacent_stacking * match_stacking_profile(interaction_ind_pair_list1[ind_i][3], interaction_ind_pair_list2[ind_j][3], scores, stk_scoring_data)
            # print('\tmatching stk')
            # print(evaluate_score)

        evaluate_score += match_nucleotide_profile(profile1.sequenceP[interaction_ind_pair_list1[ind_i][0]], profile2.sequenceP[interaction_ind_pair_list2[ind_j][0]], nucl_scoring_data)
        evaluate_score += match_nucleotide_profile(profile1.sequenceP[interaction_ind_pair_list1[ind_i][1]], profile2.sequenceP[interaction_ind_pair_list2[ind_j][1]], nucl_scoring_data)

        # print('\tmatching nucleotides')
        # print(evaluate_score)

        nucleotide_indices1.append(interaction_ind_pair_list1[ind_i][0])
        nucleotide_indices1.append(interaction_ind_pair_list1[ind_i][1])
        nucleotide_indices2.append(interaction_ind_pair_list2[ind_j][0])
        nucleotide_indices2.append(interaction_ind_pair_list2[ind_j][1])

        if significant_match(interaction_ind_pair_list1[ind_i], interaction_ind_pair_list2[ind_j]):
            #  check for bonus score of triple interaction
            if matched_positions1[interaction_ind_pair_list1[ind_i][0]] > 0 and matched_positions2[interaction_ind_pair_list2[ind_j][0]] > 0:
                evaluate_score += scores.weight_isosteric * scores.bonuses.triple_interaction_bonus

            if matched_positions1[interaction_ind_pair_list1[ind_i][1]] > 0 and matched_positions2[interaction_ind_pair_list2[ind_j][1]] > 0:
                evaluate_score += scores.weight_isosteric * scores.bonuses.triple_interaction_bonus

            # print('\tsignificant_match')
            # print(evaluate_score)

            matched_positions1[interaction_ind_pair_list1[ind_i][0]] = 1
            matched_positions1[interaction_ind_pair_list1[ind_i][1]] = 1

            matched_positions2[interaction_ind_pair_list2[ind_j][0]] = 1
            matched_positions2[interaction_ind_pair_list2[ind_j][1]] = 1

    # Asymmetric Penalty (for segment size difference between bases of pairs)
    evaluate_score += asymmetric_penalty(matched_pairs, interaction_ind_pair_list1, interaction_ind_pair_list2, scores)

    # print('\tasymmetric_penalty')
    # print(evaluate_score)

    # print(evaluate_score)
    # sys.exit()
    # # evaluate_score += local_loop_align(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2)
    # evaluate_score += global_align_segments(profile1, profile2, scores, nucl_scoring_data)

    if len(nucleotide_indices1) != len(nucleotide_indices2):
        return 0
    # elif len(nucleotide_indices1) == 0:
    # did not consider the missing bp and stack penalty at this point in case of len(nucleotide_indices1) != 0
    else:
        # if len(nucleotide_indices1) != 0:
        #     aligned_regions1, aligned_regions2 = get_aligned_regions(nucleotide_indices1, nucleotide_indices2)
        #     evaluate_score += get_missing_bp_and_stack_penalty(interaction_ind_pair_list1, interaction_ind_pair_list2, aligned_regions1, aligned_regions2, pairs_record1, pairs_record2, scores)

        # evaluate_score += clique_based_extended_sequence_alignment()






        # If empty clique
        evaluate_score += global_align_segments(profile1, profile2, scores, nucl_scoring_data)
        # print('\tbasic alignment')
        # print(evaluate_score)


    # print(evaluate_score)
    return evaluate_score


#############################################################
# Code to improve efficiency start                          #
#############################################################

def maximal_seq_score(start, end, profile, scores, nucleotide_substitution_matrix, stacking_substitution_matrix):
    max_score = 0
    for i in range(start, end):
        max_score += get_match_nucleotide_score(i, i, profile.sequenceP, profile.sequenceP, scores, nucleotide_substitution_matrix)


    for (b, e) in profile.stksP:
        if (start <= b and b <= end) and (start <= e and e <= end):
            max_score += get_adjacent_stacking_score(b, e, b, e, profile.sequenceP, profile.sequenceP, scores, stacking_substitution_matrix)

    return max_score

def calculate_alignment_score_max(estimated_clique, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, basepair_index, dBG):
    evaluate_score = 0

    interaction_count1 = len(interaction_ind_pair_list1)
    interaction_count2 = len(interaction_ind_pair_list2)

    matched_pairs = []
    nucleotide_indices1 = []
    nucleotide_indices2 = []
    matched_positions1 = [0 for i in range(len(profile1.sequenceP))]
    matched_positions2 = [0 for i in range(len(profile2.sequenceP))]

    for i in range(len(estimated_clique)):
        ind_i = estimated_clique[i] // interaction_count2
        ind_j = estimated_clique[i] % interaction_count2

        # print(ind_i, ind_j)
        matched_pairs.append((ind_i, ind_j))

        if (interaction_ind_pair_list1[ind_i][2] == 'b' and interaction_ind_pair_list2[ind_j][2] == 'b'):
            evaluate_score += scores.weight_isosteric * match_isosteric_basepair_profile(interaction_ind_pair_list1[ind_i][3], interaction_ind_pair_list2[ind_j][3], basepair_index, scores, bp_scoring_data)
            # print('\tmatching bp')
            # print(evaluate_score)
        #  match of two non-adjacent stackings. hbond may present
        elif (interaction_ind_pair_list1[ind_i][2] == 's' and interaction_ind_pair_list2[ind_j][2] == 's'):
            evaluate_score += scores.weight_nonadjacent_stacking * match_stacking_profile(interaction_ind_pair_list1[ind_i][3], interaction_ind_pair_list2[ind_j][3], scores, stk_scoring_data)
            # print('\tmatching stk')
            # print(evaluate_score)

        evaluate_score += match_nucleotide_profile(profile1.sequenceP[interaction_ind_pair_list1[ind_i][0]], profile2.sequenceP[interaction_ind_pair_list2[ind_j][0]], nucl_scoring_data)
        evaluate_score += match_nucleotide_profile(profile1.sequenceP[interaction_ind_pair_list1[ind_i][1]], profile2.sequenceP[interaction_ind_pair_list2[ind_j][1]], nucl_scoring_data)

        # print('\tmatching nucleotides')
        # print(evaluate_score)

        nucleotide_indices1.append(interaction_ind_pair_list1[ind_i][0])
        nucleotide_indices1.append(interaction_ind_pair_list1[ind_i][1])
        nucleotide_indices2.append(interaction_ind_pair_list2[ind_j][0])
        nucleotide_indices2.append(interaction_ind_pair_list2[ind_j][1])

        if significant_match(interaction_ind_pair_list1[ind_i], interaction_ind_pair_list2[ind_j]):
            #  check for bonus score of triple interaction
            if matched_positions1[interaction_ind_pair_list1[ind_i][0]] > 0 and matched_positions2[interaction_ind_pair_list2[ind_j][0]] > 0:
                evaluate_score += scores.weight_isosteric * scores.bonuses.triple_interaction_bonus

            if matched_positions1[interaction_ind_pair_list1[ind_i][1]] > 0 and matched_positions2[interaction_ind_pair_list2[ind_j][1]] > 0:
                evaluate_score += scores.weight_isosteric * scores.bonuses.triple_interaction_bonus

            # print('\tsignificant_match')
            # print(evaluate_score)

            matched_positions1[interaction_ind_pair_list1[ind_i][0]] = 1
            matched_positions1[interaction_ind_pair_list1[ind_i][1]] = 1

            matched_positions2[interaction_ind_pair_list2[ind_j][0]] = 1
            matched_positions2[interaction_ind_pair_list2[ind_j][1]] = 1

    # Asymmetric Penalty (for segment size difference between bases of pairs)
    evaluate_score += asymmetric_penalty(matched_pairs, interaction_ind_pair_list1, interaction_ind_pair_list2, scores)

    # print('\tasymmetric_penalty')
    # print(evaluate_score)

    # print(evaluate_score)
    # sys.exit()
    # # evaluate_score += local_loop_align(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2)
    # evaluate_score += global_align_segments(profile1, profile2, scores, nucl_scoring_data)

    if len(nucleotide_indices1) != len(nucleotide_indices2):
        return 0
    elif len(nucleotide_indices1) == 0:

        score1 = maximal_seq_score(0, len(profile1.sequenceP), profile1, scores, nucl_scoring_data, stk_scoring_data)
        score2 = maximal_seq_score(0, len(profile2.sequenceP), profile2, scores, nucl_scoring_data, stk_scoring_data)

        return min(score1, score2)

    else:
        score1 = maximal_seq_score(0, len(profile1.sequenceP), profile1, scores, nucl_scoring_data, stk_scoring_data)
        score2 = maximal_seq_score(0, len(profile2.sequenceP), profile2, scores, nucl_scoring_data, stk_scoring_data)

        evaluate_score += min(score1, score2)

        evaluate_score += global_align_segments(profile1, profile2, scores, nucl_scoring_data)
        # print('\tbasic alignment')
        # print(evaluate_score)


    # print(evaluate_score)
    return evaluate_score

# def calculate_alignment_score_upper_bound(current_clique, candidates, profile1, profile2, interaction_compatibility_matrix, rank_dict, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, basepair_index, dBG):
def calculate_alignment_score_upper_bound(current_clique, candidates, profile1, profile2, interaction_compatibility_matrix, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, basepair_index, dBG):
    evaluate_score_t = 0
    evaluate_score_t += calculate_alignment_score_max(current_clique, profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2, basepair_index, dBG)


    interaction_count1 = len(interaction_ind_pair_list1)
    interaction_count2 = len(interaction_ind_pair_list2)

    num_edges_candidates = 0
    degrees_candidates = {}
    for i in range(len(candidates)):
        degrees_candidates[candidates[i]] = 0

    for i in range(len(candidates)-1):
        for j in range(i+1, len(candidates)):
            if interaction_compatibility_matrix[candidates[i]][candidates[j]] == 1:
                num_edges_candidates += 1
                degrees_candidates[candidates[i]] += 1
                degrees_candidates[candidates[j]] += 1

    min_degree = interaction_count1 * interaction_count2
    for item in degrees_candidates:
        if degrees_candidates[item] < min_degree:
            min_degree = degrees_candidates[item]

    P1_max_score = []
    P2_max_score = []
    for i in range(interaction_count1):
        P1_max_score.append(0)
    for i in range(interaction_count2):
        P2_max_score.append(0)

    for i in range(len(candidates)):
        ind_i = candidates[i] // interaction_count2
        ind_j = candidates[i] % interaction_count2

        if (interaction_ind_pair_list1[ind_i][2] == 'b' and interaction_ind_pair_list2[ind_j][2] == 'b'):
            evaluate_score = scores.weight_isosteric * match_isosteric_basepair_profile(interaction_ind_pair_list1[ind_i][3], interaction_ind_pair_list2[ind_j][3], basepair_index, scores, bp_scoring_data)
        #  match of two non-adjacent stackings. hbond may present
        elif (interaction_ind_pair_list1[ind_i][2] == 's' and interaction_ind_pair_list2[ind_j][2] == 's'):
            evaluate_score = scores.weight_nonadjacent_stacking * match_stacking_profile(interaction_ind_pair_list1[ind_i][3], interaction_ind_pair_list2[ind_j][3], scores, stk_scoring_data)

        if evaluate_score > P1_max_score[ind_i]:
            P1_max_score[ind_i] = evaluate_score

        if evaluate_score > P2_max_score[ind_j]:
            P2_max_score[ind_j] = evaluate_score


    expected_matches = (math.sqrt(1 + 8 * num_edges_candidates) + 1) / 2.0

    P1_score_bound = min(interaction_count1, expected_matches)
    P2_score_bound = min(interaction_count2, expected_matches)

    P1_max_score.sort(reverse=True)
    P2_max_score.sort(reverse=True)

    P1_sum_mScore = 0
    for i in range(int(P1_score_bound)):
        P1_sum_mScore += P1_max_score[i]

    P2_sum_mScore = 0
    for i in range(int(P2_score_bound)):
        P2_sum_mScore += P2_max_score[i]

    evaluate_score_t += min(P1_sum_mScore, P2_sum_mScore)

    return evaluate_score_t






#############################################################
# Code to improve efficiency end                            #
#############################################################
def partial_global_align_segments(seq1, seq2, scores, nucl_scoring_data):
    seq_len1 = len(seq1)
    seq_len2 = len(seq2)

    rows, cols = (seq_len1 + 1, seq_len2 + 1)
    DP_table_M = [[0 for i in range(cols)] for j in range(rows)]
    DP_table_Ix = [[0 for i in range(cols)] for j in range(rows)]
    DP_table_Iy = [[0 for i in range(cols)] for j in range(rows)]

    DP_table_M[0][0] = 0
    DP_table_Ix[0][0] = scores.penalties.gap_opening_penalty
    DP_table_Iy[0][0] = scores.penalties.gap_opening_penalty

    for i in range(1, rows):
        DP_table_M[i][0] = -Infinity
        DP_table_Ix[i][0] = scores.penalties.gap_opening_penalty + i * scores.penalties.gap_extension_penalty
        DP_table_Iy[i][0] = -Infinity

    for j in range(1, cols):
        DP_table_M[0][j] = -Infinity
        DP_table_Ix[0][j] = -Infinity
        DP_table_Iy[0][j] = scores.penalties.gap_opening_penalty + j * scores.penalties.gap_extension_penalty

    for i in range(1, rows):
        for j in range(1, cols):

            DP_table_M[i][j] = get_match_nucleotide_score(i-1, j-1, seq1, seq2, scores, nucl_scoring_data) + max(
                    DP_table_M[i-1][j-1],
                    DP_table_Ix[i-1][j-1],
                    DP_table_Iy[i-1][j-1]
            )

            DP_table_Ix[i][j] = max(
                    scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_M[i-1][j],
                    # scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_Ix[i-1][j],
                    scores.penalties.gap_extension_penalty + DP_table_Ix[i-1][j]
            )

            DP_table_Iy[i][j] = max(
                    scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_M[i][j-1],
                    scores.penalties.gap_extension_penalty + DP_table_Iy[i][j-1]
                    # scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_Iy[i][j-1]
            )

    # print('tracing')
    aln_P1, aln_P2 = trace_back_alignments(seq1, seq2, DP_table_M, DP_table_Ix, DP_table_Iy, scores, nucl_scoring_data)
    return aln_P1, aln_P2

def global_align_segments(profile1, profile2, scores, nucl_scoring_data):
    # dBG = True
    # print('global')
    opt = 0
    aln_P1 = []
    aln_P2 = []
    for b in range(len(profile1.break_points)):
        # print('aligning ' + str(b))
        seq_len1 = profile1.break_points[b]
        seq_len2 = profile2.break_points[b]

        seq1 = profile1.sequenceP[:seq_len1]
        seq2 = profile2.sequenceP[:seq_len2]

        # print('seq_len')
        # print(seq_len1, seq_len2)
        # print(seq1, seq2)

        if b > 0:
            seq_len1 -= profile1.break_points[b-1]
            seq_len2 -= profile2.break_points[b-1]

            seq1 = profile1.sequenceP[profile1.break_points[b-1]:profile1.break_points[b]]
            seq2 = profile2.sequenceP[profile2.break_points[b-1]:profile2.break_points[b]]

    # print(profile1.sequenceP)
    # print(profile2.sequenceP)

        if seq_len1 == 0 or seq_len2 == 0:
            return 0

        rows, cols = (seq_len1 + 1, seq_len2 + 1)
        DP_table_M = [[0 for i in range(cols)] for j in range(rows)]
        DP_table_Ix = [[0 for i in range(cols)] for j in range(rows)]
        DP_table_Iy = [[0 for i in range(cols)] for j in range(rows)]

        DP_table_M[0][0] = 0
        DP_table_Ix[0][0] = scores.penalties.gap_opening_penalty
        DP_table_Iy[0][0] = scores.penalties.gap_opening_penalty

        for i in range(1, rows):
            DP_table_M[i][0] = -Infinity
            DP_table_Ix[i][0] = scores.penalties.gap_opening_penalty + i * scores.penalties.gap_extension_penalty
            DP_table_Iy[i][0] = -Infinity

        for j in range(1, cols):
            DP_table_M[0][j] = -Infinity
            DP_table_Ix[0][j] = -Infinity
            DP_table_Iy[0][j] = scores.penalties.gap_opening_penalty + j * scores.penalties.gap_extension_penalty

        for i in range(1, rows):
            for j in range(1, cols):

                DP_table_M[i][j] = get_match_nucleotide_score(i-1, j-1, seq1, seq2, scores, nucl_scoring_data) + max(
                        DP_table_M[i-1][j-1],
                        DP_table_Ix[i-1][j-1],
                        DP_table_Iy[i-1][j-1]
                )

                DP_table_Ix[i][j] = max(
                        scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_M[i-1][j],
                        # scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_Ix[i-1][j],
                        scores.penalties.gap_extension_penalty + DP_table_Ix[i-1][j]
                )

                DP_table_Iy[i][j] = max(
                        scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_M[i][j-1],
                        scores.penalties.gap_extension_penalty + DP_table_Iy[i][j-1]
                        # scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_Iy[i][j-1]
                )

        # print('tracing')
        aln_P1_t, aln_P2_t = trace_back_alignments(seq1, seq2, DP_table_M, DP_table_Ix, DP_table_Iy, scores, nucl_scoring_data)
        aln_P1 += list(aln_P1_t)
        aln_P2 += list(aln_P2_t)

        # print(opt)
        # print(aln_P1_t)
        # print(aln_P2_t)

        opt += max(DP_table_M[seq_len1][seq_len2], DP_table_Ix[seq_len1][seq_len2], DP_table_Iy[seq_len1][seq_len2])

    # print(aln_P1)
    # print(aln_P2)
    # print_formatted_profile(aln_P1)
    # print_formatted_profile(aln_P2)
    # for i, nucl_dict in enumerate(aln_P1):
    #     if dBG:
    #         logger.debug(str(i) + '\t')
    #     else:
    #         print(i, end='\t')
    #     nucl_dict = dict(sorted(nucl_dict.items(), key=lambda item: item[1], reverse=True))
    #     for j, hashed_nucl in enumerate(nucl_dict):
    #         if j == len(nucl_dict) - 1:
    #             if dBG:
    #                 logger.debug(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]))
    #             else:
    #                 print(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]))
    #         else:
    #             if dBG:
    #                 logger.debug(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]) + ', ')
    #             else:
    #                 print(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]), end=', ')
    # for i, nucl_dict in enumerate(aln_P2):
    #     if dBG:
    #         logger.debug(str(i) + '\t')
    #     else:
    #         print(i, end='\t')
    #     nucl_dict = dict(sorted(nucl_dict.items(), key=lambda item: item[1], reverse=True))
    #     for j, hashed_nucl in enumerate(nucl_dict):
    #         if j == len(nucl_dict) - 1:
    #             if dBG:
    #                 logger.debug(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]))
    #             else:
    #                 print(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]))
    #         else:
    #             if dBG:
    #                 logger.debug(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]) + ', ')
    #             else:
    #                 print(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]), end=', ')
    # sys.exit()
    return opt

# def global_align_segments(profile1, profile2, scores, nucl_scoring_data):
    
#     seq_len1 = len(profile1.sequenceP)
#     seq_len2 = len(profile2.sequenceP)

#     # print(profile1.sequenceP)
#     # print(profile2.sequenceP)

#     if seq_len1 == 0 or seq_len2 == 0:
#         return 0

#     rows, cols = (seq_len1 + 1, seq_len2 + 1)
#     DP_table_M = [[0 for i in range(cols)] for j in range(rows)]
#     DP_table_Ix = [[0 for i in range(cols)] for j in range(rows)]
#     DP_table_Iy = [[0 for i in range(cols)] for j in range(rows)]

#     DP_table_M[0][0] = 0
#     DP_table_Ix[0][0] = scores.penalties.gap_opening_penalty
#     DP_table_Iy[0][0] = scores.penalties.gap_opening_penalty

#     for i in range(1, rows):
#         DP_table_M[i][0] = -Infinity
#         DP_table_Ix[i][0] = scores.penalties.gap_opening_penalty + i * scores.penalties.gap_extension_penalty
#         DP_table_Iy[i][0] = -Infinity

#     for j in range(1, cols):
#         DP_table_M[0][j] = -Infinity
#         DP_table_Ix[0][j] = -Infinity
#         DP_table_Iy[0][j] = scores.penalties.gap_opening_penalty + j * scores.penalties.gap_extension_penalty

#     for i in range(1, rows):
#         for j in range(1, cols):

#             DP_table_M[i][j] = get_match_nucleotide_score(i-1, j-1, profile1.sequenceP, profile2.sequenceP, scores, nucl_scoring_data) + max(
#                     DP_table_M[i-1][j-1],
#                     DP_table_Ix[i-1][j-1],
#                     DP_table_Iy[i-1][j-1]
#             )

#             DP_table_Ix[i][j] = max(
#                     scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_M[i-1][j],
#                     # scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_Ix[i-1][j],
#                     scores.penalties.gap_extension_penalty + DP_table_Ix[i-1][j]
#             )

#             DP_table_Iy[i][j] = max(
#                     scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_M[i][j-1],
#                     scores.penalties.gap_extension_penalty + DP_table_Iy[i][j-1]
#                     # scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_Iy[i][j-1]
#             )

#     # print('tracing')
#     aln_P1, aln_P2 = trace_back_alignments(profile1, profile2, DP_table_M, DP_table_Ix, DP_table_Iy, scores, nucl_scoring_data)
#     opt = max(DP_table_M[seq_len1][seq_len2], DP_table_Ix[seq_len1][seq_len2], DP_table_Iy[seq_len1][seq_len2])

#     # print(opt)
#     # print(aln_P1)
#     # print(aln_P2)
#     return opt

def trace_back_alignments(seq1, seq2, DP_table_M, DP_table_Ix, DP_table_Iy, scores, nucl_scoring_data):
    # print('DP_table_M')
    # print_2d_list(DP_table_M)
    # print('DP_table_Ix')
    # print_2d_list(DP_table_Ix)
    # print('DP_table_Iy')
    # print_2d_list(DP_table_Iy)

    
    # seq_len1 = len(profile1.sequenceP)
    # seq_len2 = len(profile2.sequenceP)
    seq_len1 = len(seq1)
    seq_len2 = len(seq2)

    aln_P1 = []
    aln_P2 = []

    i = seq_len1
    j = seq_len2
    # print(i, j)
    # sys.exit()
    while (i>0 and j>0):
        match_score = get_match_nucleotide_score(i-1, j-1, seq1, seq2, scores, nucl_scoring_data)
        if (DP_table_M[i][j] == DP_table_M[i-1][j-1] + match_score):
            aln_P1.append(seq1[i-1])
            aln_P2.append(seq2[j-1])
            i -= 1
            j -= 1

        elif (DP_table_M[i][j] == DP_table_Ix[i-1][j-1] + match_score):
            # aln_P1.append({'-':1})
            # aln_P2.append(profile2.sequenceP[j-1])
            # j -= 1
            aln_P1.append(seq1[i-1])
            aln_P2.append({'-':1})
            i -= 1

        elif (DP_table_M[i][j] == DP_table_Iy[i-1][j-1] + match_score):
            # aln_P1.append(profile1.sequenceP[i-1])
            # aln_P2.append({'-':1})
            # i -= 1
            aln_P1.append({'-':1})
            aln_P2.append(seq2[j-1])
            j -= 1

    while (i>0):
        aln_P1.append(seq1[i-1])
        aln_P2.append({'-':1})
        i -= 1
    while (j>0):
        aln_P1.append({'-':1})
        aln_P2.append(seq2[j-1])
        j -= 1

    # print(aln_P1)
    # print(aln_P2)
    aln_P1.reverse()
    aln_P2.reverse()
    # print('reversed')
    # print(aln_P1)
    # print(aln_P2)
    return aln_P1, aln_P2

# def trace_back_alignments(profile1, profile2, DP_table_M, DP_table_Ix, DP_table_Iy, scores, nucl_scoring_data):
#     # print('DP_table_M')
#     # print_2d_list(DP_table_M)
#     # print('DP_table_Ix')
#     # print_2d_list(DP_table_Ix)
#     # print('DP_table_Iy')
#     # print_2d_list(DP_table_Iy)

    
#     seq_len1 = len(profile1.sequenceP)
#     seq_len2 = len(profile2.sequenceP)

#     aln_P1 = []
#     aln_P2 = []

#     i = seq_len1
#     j = seq_len2
#     # print(i, j)
#     # sys.exit()
#     while (i>0 and j>0):
#         match_score = get_match_nucleotide_score(i-1, j-1, profile1.sequenceP, profile2.sequenceP, scores, nucl_scoring_data)
#         if (DP_table_M[i][j] == DP_table_M[i-1][j-1] + match_score):
#             aln_P1.append(profile1.sequenceP[i-1])
#             aln_P2.append(profile2.sequenceP[j-1])
#             i -= 1
#             j -= 1

#         elif (DP_table_M[i][j] == DP_table_Ix[i-1][j-1] + match_score):
#             # aln_P1.append({'-':1})
#             # aln_P2.append(profile2.sequenceP[j-1])
#             # j -= 1
#             aln_P1.append(profile1.sequenceP[i-1])
#             aln_P2.append({'-':1})
#             i -= 1

#         elif (DP_table_M[i][j] == DP_table_Iy[i-1][j-1] + match_score):
#             # aln_P1.append(profile1.sequenceP[i-1])
#             # aln_P2.append({'-':1})
#             # i -= 1
#             aln_P1.append({'-':1})
#             aln_P2.append(profile2.sequenceP[j-1])
#             j -= 1

#     while (i>0):
#         aln_P1.append(profile1.sequenceP[i-1])
#         aln_P2.append({'-':1})
#         i -= 1
#     while (j>0):
#         aln_P1.append({'-':1})
#         aln_P2.append(profile2.sequenceP[j-1])
#         j -= 1

#     # print(aln_P1)
#     # print(aln_P2)
#     aln_P1.reverse()
#     aln_P2.reverse()
#     # print('reversed')
#     # print(aln_P1)
#     # print(aln_P2)
#     return aln_P1, aln_P2


def local_loop_align(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, nucleotide_deletion_penalties1, nucleotide_deletion_penalties2):
    #  compute the semi-local alignment of the two sequences
    seq_len1 = len(profile1.sequenceP)
    seq_len2 = len(profile2.sequenceP)

    # print(profile1.sequenceP)
    # print(profile2.sequenceP)

    # sys.exit()

    #/********** Base Case **********/
    # If either or both segment is empty
    # For local alignment, alignment should be empty (Edited by Shahidul)
    if seq_len1 == 0 or seq_len2 == 0:
        return 0

    # #/********** Allocate space **********/
    rows, cols = (seq_len1, seq_len2)
    DP_table = [[0 for i in range(cols)] for j in range(rows)]
    # # Record the previous indices of (i, j) that lead to the best alignment
    record_i = [[0 for i in range(cols)] for j in range(rows)]
    record_j = [[0 for i in range(cols)] for j in range(rows)]
    # # Record the aligned region for (i, j) that lead to the best alignment
    align_region_i = [[[] for i in range(cols)] for j in range(rows)]
    align_region_j = [[[] for i in range(cols)] for j in range(rows)]
    

    #/********** Matrix initialization with boundary case **********/
    #  the boundary case for initial gaps are not required, as it is local alignment, and the first row/column would be zeros anyway
    for i in range(seq_len1):
        DP_table[i][0] = scores.weight_sequence * match_nucleotide_profile(profile1.sequenceP[i], profile2.sequenceP[0], nucl_scoring_data)
        record_i[i][0] = -1
        record_j[i][0] = -1
        add_a_region(align_region_i[i][0], i, i)
        add_a_region(align_region_j[i][0], 0, 0)
        # align_region_i[i][0].append((i, i))
        # align_region_j[i][0].append((0, 0))

    for j in range(1, seq_len2):
        DP_table[0][j] = scores.weight_sequence * match_nucleotide_profile(profile1.sequenceP[0], profile2.sequenceP[j], nucl_scoring_data)
        record_i[0][j] = -1
        record_j[0][j] = -1
        add_a_region(align_region_i[0][j], 0, 0)
        add_a_region(align_region_j[0][j], j, j)
        # align_region_i[0][j].append((0, 0))
        # align_region_j[0][j].append((j, j))

    print('Initital DP_table:')
    # print(DP_tasble)
    for i in range(seq_len1):
        for j in range(seq_len2):
            print(DP_table[i][j], end=', ')
        print('')

    # sys.exit()
    #/********** General Cases **********/
    for i in range(1, seq_len1):
        for j in range(1, seq_len2):

            max_score = -9999999
            temp_score = 0
            max_i_source = -1
            max_j_source = -1
            max_region_i = []
            max_region_j = []
            current_region_i = []
            current_region_j = []

            #  go over all entries in 0...i - 1, j - 1
            for k in range(i):
            # for(k = 0; k < i; ++ k)  {
                #  the previous alignment ends at k, j - 1
                temp_score = DP_table[k][j - 1]
                print(temp_score)
                temp_score += get_adjacent_stacking_score(k, i, j - 1, j, profile1.stksP, profile2.stksP, scores, stk_scoring_data)
                print(temp_score)
                temp_score += get_nucleotide_deletion_penalty(k, i, scores, nucleotide_deletion_penalties1, profile1.sequenceP)
                print(temp_score)
                temp_score += get_match_nucleotide_score(i, j, profile1.sequenceP, profile2.sequenceP, scores, nucl_scoring_data)
                print(temp_score)
                # add missing BP penalty for the index range (k + 1, i), (j, j) [Shahidul]
                temp_score += bp_and_stack_penalty_inclusion_exclusion(profile1.break_points, profile2.break_points, interaction_ind_pair_list1, interaction_ind_pair_list2, k, i, j - 1, j, 0, 0, current_region_i, current_region_j, align_region_i, align_region_j, scores)
                print(temp_score)

                print("temp score 1: " + str(temp_score))
                # if is_better_alignment(temp_score, max_score, current_region_i, current_region_j, max_region_i, max_region_j):
                if is_better_alignment(temp_score, max_score):
                    max_score = temp_score
                    max_i_source = k
                    max_j_source = j - 1
                    max_region_i = current_region_i
                    max_region_j = current_region_j
                
            
            #  go over all entries in i - 1, 0...j - 1
            for k in range(j):
            # for(k = 0; k < j; ++ k)  {
                #  the previous alignment ends at i - 1, k
                temp_score = DP_table[i - 1][k];
                temp_score += get_adjacent_stacking_score(i - 1, i, k, j, profile1.sequenceP, profile2.sequenceP, scores, stk_scoring_data)
                print(temp_score)
                temp_score += get_nucleotide_deletion_penalty(k, j, scores, nucleotide_deletion_penalties2, profile2.sequenceP)
                print(temp_score)
                temp_score += get_match_nucleotide_score(i, j, profile1.sequenceP, profile2.sequenceP, scores, nucl_scoring_data)
                print(temp_score)
                # add missing BP penalty for the index range (i, i), (k+1, j) [Shahidul]
                temp_score += bp_and_stack_penalty_inclusion_exclusion(profile1.break_points, profile2.break_points, interaction_ind_pair_list1, interaction_ind_pair_list2, i - 1, i, k, j, 0, 0, current_region_i, current_region_j, align_region_i, align_region_j, scores)
                print("temp score 2: " + str(temp_score))
                # if(is_better_alignment(temp_score, max_score, current_region_i, current_region_j, max_region_i, max_region_j))  {
                if is_better_alignment(temp_score, max_score):
                    max_score = temp_score;
                    max_i_source = i - 1;
                    max_j_source = k;
                    max_region_i = current_region_i;
                    max_region_j = current_region_j;
                # }
            # }


            initial_match_score = get_match_nucleotide_score(i, j, profile1.sequenceP, profile2.sequenceP, scores, nucl_scoring_data)
            initial_region_i = []
            initial_region_j = []
            add_a_region(initial_region_i, i, i)
            add_a_region(initial_region_j, j, j)


            # no missing BP penalty for the index range (i, i), (j, j), as only one index cannot have interaction [Shahidul]
            # if(is_better_alignment(initial_match_score, max_score, current_region_i, current_region_j, max_region_i, max_region_j))  {
            print("temp score 4: " + str(temp_score))
            if is_better_alignment(initial_match_score, max_score):
                max_score = initial_match_score
                max_region_i = initial_region_i
                max_region_j = initial_region_j
                max_i_source = -1
                max_j_source = -1
            
            DP_table[i][j] = max_score
            record_i[i][j] = max_i_source
            record_j[i][j] = max_j_source
            align_region_i[i][j] = max_region_i
            align_region_j[i][j] = max_region_j
        
    
    #  finds the maximum alignment score from the entire matrix
    # int maximum_align_score = 0, max_score_i = 0, max_score_j = 0;
    maximum_align_score = 0
    max_score_i = 0
    max_score_j = 0
    # for(i = 0; i < (int) M1_seq[seq_index].length(); ++ i)  {
    #     for(j = 0; j < (int) M2_seq[seq_index].length(); ++ j)  {
    for i in range(seq_len1):
        for j in range(seq_len2):
            if DP_table[i][j] > maximum_align_score:
                maximum_align_score = DP_table[i][j]
                max_score_i = i
                max_score_j = j
        
    

    # Trace code call
    # trace_fixed_end_open_beginning(0, 0, max_score_i, max_score_j, record_i, record_j, M1_alignment, M2_alignment, M1_nuc, M2_nuc, DP_table);

    print('Current DP_table:')
    # print(DP_tasble)
    for i in range(seq_len1):
        for j in range(seq_len2):
            print(DP_table[i][j], end=', ')
        print('')

    print(align_region_i)
    print(align_region_j)

    return maximum_align_score

def print_2d_list(lst):
    for i in range(len(lst)):
        for j in range(len(lst[i])):
            print(lst[i][j], end=', ')
        print('')

def get_non_gapped_aligned_regions(nucleotide_indices1, nucleotide_indices2, dBG):
    len1 = len(nucleotide_indices1)
    len2 = len(nucleotide_indices2)
    if len1 != len2:
        print('Aligned nucleotides length not equal. Please check.')
        sys.exit()

    aligned_regions1 = []
    aligned_regions2 = []

    if len1 > 0:
        aligned_regions1.append((nucleotide_indices1[0], nucleotide_indices1[0]))
    if len2 > 0:
        aligned_regions2.append((nucleotide_indices2[0], nucleotide_indices2[0]))

    for i in range(1, len1):
        if nucleotide_indices1[i] > nucleotide_indices1[i-1] and nucleotide_indices2[i] > nucleotide_indices2[i-1]:
            if aligned_regions1[-1][1] + 1 == nucleotide_indices1[i] and aligned_regions2[-1][1] + 1 == nucleotide_indices2[i]:
                s, e = aligned_regions1[-1]
                aligned_regions1[-1] = (s, nucleotide_indices1[i])
                s, e = aligned_regions2[-1]
                aligned_regions2[-1] = (s, nucleotide_indices2[i])
            else:
                aligned_regions1.append((nucleotide_indices1[i], nucleotide_indices1[i]))
                aligned_regions2.append((nucleotide_indices2[i], nucleotide_indices2[i]))
        else:
            if (nucleotide_indices1[i] > nucleotide_indices1[i-1] and nucleotide_indices2[i] == nucleotide_indices2[i-1]) or \
            (nucleotide_indices1[i] == nucleotide_indices1[i-1] and nucleotide_indices2[i] > nucleotide_indices2[i-1]):
                if dBG:
                    logger.debug('why both are not greater than the previous value? multiple nucleotide mapped to one?')
                    # print('why both are not greater than the previous value? multiple nucleotide mapped to one?')
                # sys.exit()

            # if nucleotide_indices1[i] > nucleotide_indices1[i-1]:
            #     if aligned_regions1[-1][1] + 1 == nucleotide_indices1[i]:
            #         s, e = aligned_regions1[-1]
            #         aligned_regions1[-1] = (s, nucleotide_indices1[i])
            #     else:
            #         aligned_regions1.append((nucleotide_indices1[i], nucleotide_indices1[i]))

            # if nucleotide_indices2[i] > nucleotide_indices2[i-1]:
            #     if aligned_regions2[-1][1] + 1 == nucleotide_indices2[i]:
            #         s, e = aligned_regions2[-1]
            #         aligned_regions2[-1] = (s, nucleotide_indices2[i])
            #     else:
            #         aligned_regions2.append((nucleotide_indices2[i], nucleotide_indices2[i]))

    return aligned_regions1, aligned_regions2

# def update_index_list_considering_breakpoints(ind1_list, ind2_list, profile1, profile2):
#     break_points = []
#     if len(profile1.break_points) != len(profile2.break_points):
#         print('Break points mismatch!')
#         sys.exit()
#     else:
#         for i in range(len(profile1.break_points)):
#             b_point1 = profile1.break_points[i]
#             b_point2 = profile2.break_points[i]

#             if (b_point1 < len(profile1.sequenceP) and b_point1 not in ind1_list) or \
#             (b_point2 < len(profile2.sequenceP) and b_point2 not in ind2_list):
#                 print('break point index not found in ind_lists')
#                 sys.exit()

#             elif b_point1 == len(profile1.sequenceP) and b_point2 == len(profile2.sequenceP):
#                 break_points.append(len(ind1_list))
#             else:
#                 b1 = ind1_list.index(b_point1)
#                 b2 = ind2_list.index(b_point2)
#                 print(b1)
#                 print(b2)
#                 if b1 != b2:
#                     print('break point index not matches!')
#                     print('prev ind_list')
#                     print(ind1_list)
#                     print(ind2_list)
#                     max_b = max(b1, b2)
#                     if max_b == b1:
#                         ind2_list = ind2_list[:b2+1] + ['-' for k in range(b2+1,b1+1)] + ind2_list[b2+1:]
#                     else:
#                         #max_b == b2
#                         ind1_list = ind1_list[:b1+1] + ['-' for k in range(b1+1,b2+1)] + ind1_list[b1+1:]

#                     break_points.append(max_b)
#                     print('updated ind_list')
#                     print(ind1_list)
#                     print(ind2_list)
#                     sys.exit()
#                 else:
#                     break_points.append(ind1_list.index(b_point1))

#     return ind1_list, ind2_list, break_points

def filter_aligned_regions_within_range(aligned_regions, seg_s, seg_e):
    filtered_regions = []
    for region in aligned_regions:
        s, e = region
        if (seg_s <= s and s <= seg_e) and (seg_s <= e and e <= seg_e):
            filtered_regions.append(region)

    return filtered_regions

def update_ind_list(ind_list, aln_P, p):
    k = 0
    for j in range(len(aln_P)):
        if '-' in aln_P[j] and aln_P[j]['-'] == 1: # new gap
            ind_list.append('-')
        else:
            ind_list.append(p + k)
            k += 1

    return ind_list

def get_gapped_index_mapping_v2(aligned_regions1, aligned_regions2, profile1, profile2, scores, nucl_scoring_data, dBG):
    len1 = len(aligned_regions1)
    len2 = len(aligned_regions2)

    p1_len = len(profile1.sequenceP)
    p2_len = len(profile2.sequenceP)

    if len1 != len2:
        print('Aligned regions are not equal, please check.')
        sys.exit()

    ind1_list = []
    ind2_list = []

    # print(profile1.break_points)
    # print(profile2.break_points)
    

    prev1 = 0
    prev2 = 0
    for break_point_i in range(len(profile1.break_points)):
        partial_ind1_list = []
        partial_ind2_list = []

        seg1_s = prev1
        seg1_e = profile1.break_points[break_point_i] - 1

        seg2_s = prev2
        seg2_e = profile2.break_points[break_point_i] - 1

        # print(seg1_s, seg1_e)
        # print(seg2_s, seg2_e)

        filtered_aligned_regions1 = filter_aligned_regions_within_range(aligned_regions1, seg1_s, seg1_e)
        filtered_aligned_regions2 = filter_aligned_regions_within_range(aligned_regions2, seg2_s, seg2_e)

        if len(filtered_aligned_regions1) != len(filtered_aligned_regions2):
            print('Filtered aligned regions are not equal, please check.')
            sys.exit()

        p1 = seg1_s
        p2 = seg2_s

        for filt_aln_reg_i in range(len(filtered_aligned_regions1)):
            region1_s, region1_e = filtered_aligned_regions1[filt_aln_reg_i]
            region2_s, region2_e = filtered_aligned_regions2[filt_aln_reg_i]

            # if filt_aln_reg_i == 0:
            if p1 < region1_s or p2 < region2_s:
                aln_P1a, aln_P2a = partial_global_align_segments(profile1.sequenceP[p1:region1_s], profile2.sequenceP[p2:region2_s], scores, nucl_scoring_data)
                # print(aln_P1a)
                # print(aln_P2a)
                ind1_list = update_ind_list(ind1_list, aln_P1a, p1)
                ind2_list = update_ind_list(ind2_list, aln_P2a, p2)

                if dBG:
                    print(ind1_list)
                    print(ind2_list)


            for i in range(region1_s, region1_e+1):
                # partial_ind1_list.append(i)
                ind1_list.append(i)
            for i in range(region2_s, region2_e+1):
                # partial_ind2_list.append(i)
                ind2_list.append(i)

            if dBG:
                print(ind1_list)
                print(ind2_list)


            p1 = region1_e + 1
            p2 = region2_e + 1

        # partial_ind1_list = sorted(partial_ind1_list)
        # partial_ind2_list = sorted(partial_ind2_list)

        # unaligned_regions1 = find_unaligned_regions(seg1_s, seg1_e, partial_ind1_list)
        # unaligned_regions2 = find_unaligned_regions(seg2_s, seg2_e, partial_ind2_list)

        # if len(unaligned_regions1) != len(unaligned_regions2):
        #     print('Unaligned regions are not equal, please check')
        #     sys.exit()
        if p1 <= seg1_e or p2 <= seg2_e:
            aln_P1a, aln_P2a = partial_global_align_segments(profile1.sequenceP[p1:seg1_e+1], profile2.sequenceP[p2:seg2_e+1], scores, nucl_scoring_data)
            # print(aln_P1a)
            # print(aln_P2a)
            ind1_list = update_ind_list(ind1_list, aln_P1a, p1)
            ind2_list = update_ind_list(ind2_list, aln_P2a, p2)

            if dBG:
                print(ind1_list)
                print(ind2_list)

        # ind1_list += partial_ind1_list
        # ind2_list += partial_ind2_list
        prev1 = profile1.break_points[break_point_i]
        prev2 = profile2.break_points[break_point_i]

    if dBG:
        print(ind1_list)
        print(ind2_list)

    # sys.exit()

    return ind1_list, ind2_list

def get_gapped_index_mapping(aligned_regions1, aligned_regions2, profile1, profile2, scores, nucl_scoring_data, dBG):
    len1 = len(aligned_regions1)
    len2 = len(aligned_regions2)

    p1_len = len(profile1.sequenceP)
    p2_len = len(profile2.sequenceP)

    if len1 != len2:
        print('Aligned regions are not equal, please check.')
        sys.exit()

    ind1_list = []
    ind2_list = []

    # if aligned_regions1[0][0] > 0:
    #     for i in range(aligned_regions1[0][0]):
    #         ind1_list.append('-')

    # if aligned_regions2[0][0] > 0:
    #     for i in range(aligned_regions2[0][0]):
    #         ind2_list.append('-')

    prev1 = -1
    prev2 = -1
    for i in range(len1):
        if aligned_regions1[i][0] > prev1 + 1 and aligned_regions2[i][0] > prev2 + 1:
            # both has missing index mapping
            # print('case 1')
            m1 = aligned_regions1[i][0] - prev1 - 1
            m2 = aligned_regions2[i][0] - prev2 - 1

            # print(m1, m2)
            # print(ind1_list)
            # print(ind2_list)

            # if m1 == m2:
            #     # good news
            #     for j in range(prev2 + 1, aligned_regions2[i][0]):
            #         ind1_list.append('-')
            #     for j in range(prev1 + 1, aligned_regions1[i][0]):
            #         ind1_list.append(j)

                
            #     for j in range(prev2 + 1, aligned_regions2[i][0]):
            #         ind2_list.append(j)
            #     for j in range(prev1 + 1, aligned_regions1[i][0]):
            #         ind2_list.append('-')

            # else:
            # print(m1, m2)
            aln_P1 = []
            aln_P2 = []
            # aln_P1, aln_P2 = partial_global_align_segments(profile1.sequenceP[prev1+1:aligned_regions1[i][0]], profile2.sequenceP[prev2+1:aligned_regions2[i][0]], scores, nucl_scoring_data)

            st1 = prev1+1
            st2 = prev2+1
            break_p_ind1 = get_break_point_ind(profile1.break_points, st1)
            break_p_ind2 = get_break_point_ind(profile2.break_points, st2)
            print('break_p_ind1')
            print(break_p_ind1)
            print('break_p_ind2')
            print(break_p_ind2)
            
            if dBG:
                logger.debug('Break points1: ' + str(profile1.break_points))
                print(profile1.break_points)
                logger.debug('Break points2: ' + str(profile2.break_points))
                print(profile2.break_points)

                logger.debug('Need to align from: ')
                # print('Need to align from ')
                logger.debug(str((prev1+1,aligned_regions1[i][0])))
                logger.debug(str((prev2+1,aligned_regions2[i][0])))
                # print(prev1+1,aligned_regions1[i][0])
                # print(prev2+1,aligned_regions2[i][0])

            if break_p_ind1 == break_p_ind2:
                print('here1')
                # print(break_p_ind1)
                en1 = min(profile1.break_points[break_p_ind1], aligned_regions1[i][0])
                en2 = min(profile2.break_points[break_p_ind2], aligned_regions2[i][0])
                print(st1, st2)
                print(en1, en2)
                aln_P1a, aln_P2a = partial_global_align_segments(profile1.sequenceP[st1:en1], profile2.sequenceP[st2:en2], scores, nucl_scoring_data)

                for x in aln_P1a:
                    aln_P1.append(x)
                for x in aln_P2a:
                    aln_P2.append(x)

                if en1 < aligned_regions1[i][0] or en2 < aligned_regions2[i][0]:
                    st1 = en1
                    st2 = en2
                    en1 = aligned_regions1[i][0]
                    en2 = aligned_regions2[i][0]
                    print(st1, st2)
                    print(en1, en2)
                    aln_P1b, aln_P2b = partial_global_align_segments(profile1.sequenceP[st1:en1], profile2.sequenceP[st2:en2], scores, nucl_scoring_data)

                    for x in aln_P1b:
                        aln_P1.append(x)
                    for x in aln_P2b:
                        aln_P2.append(x)
            else:
                print('here2')
                if break_p_ind1 < break_p_ind2:
                    en1 = min(profile1.break_points[break_p_ind1], aligned_regions1[i][0])
                    en2 = profile2.break_points[break_p_ind2-1]

                    aln_P1a, aln_P2a = partial_global_align_segments(profile1.sequenceP[st1:en1], profile2.sequenceP[st2:en2], scores, nucl_scoring_data)

                    for x in aln_P1a:
                        aln_P1.append(x)
                    for x in aln_P2a:
                        aln_P2.append(x)

                    st1 = en1
                    st2 = st2
                    en1 = aligned_regions1[i][0]
                    en2 = aligned_regions2[i][0]

                    aln_P1b, aln_P2b = partial_global_align_segments(profile1.sequenceP[st1:en1], profile2.sequenceP[st2:en2], scores, nucl_scoring_data)

                    for x in aln_P1b:
                        aln_P1.append(x)
                    for x in aln_P2b:
                        aln_P2.append(x)

                else:
                    en1 = profile1.break_points[break_p_ind1-1]
                    en2 = min(profile2.break_points[break_p_ind2], aligned_regions2[i][0])

                    aln_P1a, aln_P2a = partial_global_align_segments(profile1.sequenceP[st1:en1], profile2.sequenceP[st2:en2], scores, nucl_scoring_data)

                    for x in aln_P1a:
                        aln_P1.append(x)
                    for x in aln_P2a:
                        aln_P2.append(x)

                    st1 = st1
                    st2 = en2
                    en1 = aligned_regions1[i][0]
                    en2 = aligned_regions2[i][0]

                    aln_P1b, aln_P2b = partial_global_align_segments(profile1.sequenceP[st1:en1], profile2.sequenceP[st2:en2], scores, nucl_scoring_data)

                    for x in aln_P1b:
                        aln_P1.append(x)
                    for x in aln_P2b:
                        aln_P2.append(x)

                if dBG:
                    logger.debug('starting from two different segments, check')
                    # print('starting from two different segments, check')

            # print(aln_P1)
            # print(aln_P2)

            print('before')
            print(ind1_list)
            print(ind2_list)

            kk1 = 1
            kk2 = 1
            for j in range(len(aln_P1)):
                if '-' in aln_P1[j] and aln_P1[j]['-'] == 1: # new gap
                    ind1_list.append('-')
                else:
                    ind1_list.append(prev1 + kk1)
                    kk1 += 1

                if '-' in aln_P2[j] and aln_P2[j]['-'] == 1: # new gap
                    ind2_list.append('-')
                else:
                    ind2_list.append(prev2 + kk2)
                    kk2 += 1

            print('after')
            print(ind1_list)
            print(ind2_list)


            # print('Unfortunately the missing mapped portion is not equal, please implement partial alignment.')
            # sys.exit()

        elif aligned_regions1[i][0] > prev1 + 1:
            # print('case 2')
            for j in range(prev1 + 1, aligned_regions1[i][0]):
                ind1_list.append(j)
                ind2_list.append('-')

        elif aligned_regions2[i][0] > prev2 + 1:
            # print('case 3')
            for j in range(prev2 + 1, aligned_regions2[i][0]):
                ind1_list.append('-')
                ind2_list.append(j)
        # else:
        # print('case default')
        s, e = aligned_regions1[i]
        for j in range(s, e + 1):
            ind1_list.append(j)

        s, e = aligned_regions2[i]
        for j in range(s, e + 1):
            ind2_list.append(j)

        prev1 = aligned_regions1[i][1]
        prev2 = aligned_regions2[i][1]

        print('current status')
        print(ind1_list)
        print(ind2_list)


    
    if prev1 < p1_len - 1 and prev2 < p2_len - 1:
        aln_P1 = []
        aln_P2 = []
        st1 = prev1 + 1
        en1 = p1_len
        st2 = prev2 + 1
        en2 = p2_len
        aln_P1b, aln_P2b = partial_global_align_segments(profile1.sequenceP[st1:en1], profile2.sequenceP[st2:en2], scores, nucl_scoring_data)

        for x in aln_P1b:
            aln_P1.append(x)
        for x in aln_P2b:
            aln_P2.append(x)

        kk1 = 1
        kk2 = 1
        for j in range(len(aln_P1)):
            if '-' in aln_P1[j] and aln_P1[j]['-'] == 1: # new gap
                ind1_list.append('-')
            else:
                ind1_list.append(prev1 + kk1)
                kk1 += 1

            if '-' in aln_P2[j] and aln_P2[j]['-'] == 1: # new gap
                ind2_list.append('-')
            else:
                ind2_list.append(prev2 + kk2)
                kk2 += 1

    elif prev1 < p1_len - 1:
        for i in range(prev1 + 1, p1_len):
            ind1_list.append(i)
            ind2_list.append('-')

    elif prev2 < p2_len - 1:
        for i in range(prev2 + 1, p2_len):
            ind1_list.append('-')
            ind2_list.append(i)

    # ind1_list, ind2_list, break_points = update_index_list_considering_breakpoints(ind1_list, ind2_list, profile1, profile2)
    # return ind1_list, ind2_list, break_points
    return ind1_list, ind2_list

def generate_index_mapping(profile1, profile2, maximal_clique, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, nucl_scoring_data, dBG):
    # print('Profile1')
    # print_formatted_profile(profile1)
    # print('Profile2')
    # print_formatted_profile(profile2)
    # print(maximal_clique, interaction_ind_pair_list1, interaction_ind_pair_list2)
    nucleotide_indices1, nucleotide_indices2 = get_aligned_nucs(maximal_clique, interaction_ind_pair_list1, interaction_ind_pair_list2)
    # if len(maximal_clique) == 0:
    #     nucleotide_indices1, nucleotide_indices2 = get_aligned_nucs_from_global_align_segments(profile1, profile2, scores, nucl_scoring_data)
    
    # nucleotide_indices1 = [0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 5,  8, 8, 8, 9, 9, 10]
    # nucleotide_indices2 = [1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 6,  9, 9, 9, 10, 10, 11]

    if dBG:
        logger.debug('nucleotide_indices1: ' + str(nucleotide_indices1))
        logger.debug('nucleotide_indices2: ' + str(nucleotide_indices2))
        # print('nucleotide_indices')
        # print(nucleotide_indices1)
        # print(nucleotide_indices2)
        # sys.exit()

    s2t = {}
    s2t_r = []
    t2s = {}
    t2s_r = []
    for i, ind in enumerate(nucleotide_indices1):
        if ind not in s2t and nucleotide_indices2[i] not in s2t_r:
            s2t[ind] = nucleotide_indices2[i]
            s2t_r.append(nucleotide_indices2[i])
    for i, ind in enumerate(nucleotide_indices2):
        if ind not in t2s and nucleotide_indices1[i] not in t2s_r:
            t2s[ind] = nucleotide_indices1[i]
            t2s_r.append(nucleotide_indices1[i])


    if dBG:
        logger.debug('s2t: ' + str(s2t))
        logger.debug('t2s: ' + str(t2s))
        # print('s2t')
        # print(s2t)
        # print('t2s')
        # print(t2s)

    # s2t_indices = sorted(s2t)
    # t2s_indices = sorted(t2s)

    # ind1_list = []
    # if t2s_indices[0] > 0:
    #     for i in range(t2s_indices[0]):
    #         ind1_list.append('-')
    # prev = -1
    # for i in range(t2s_indices[0], t2s_indices[-1]+1):
    #     if i in t2s:
    #         if prev + 1 < t2s[i]:
    #             for k in range(prev + 1, t2s[i]):
    #                 ind1_list.append(k)
    #         ind1_list.append(t2s[i])
    #         prev = t2s[i]
    #     else:
    #         ind1_list.append('-')
    # if t2s[t2s_indices[-1]] < len(profile1.sequenceP) - 1:
    #     for i in range(t2s_indices[-1], len(profile1.sequenceP)):
    #         ind1_list.append('-')

    # ind2_list = []
    # if s2t_indices[0] > 0:
    #     for i in range(s2t_indices[0]):
    #         ind2_list.append('-')
    # prev = -1
    # for i in range(s2t_indices[0], s2t_indices[-1]+1):
    #     if i in s2t:
    #         if prev + 1 < s2t[i]:
    #             for k in range(prev + 1, s2t[i]):
    #                 ind2_list.append(k)
    #         ind2_list.append(s2t[i])
    #         prev = s2t[i]
    #     else:
    #         ind2_list.append('-')
    # if s2t[s2t_indices[-1]] < len(profile2.sequenceP) - 1:
    #     for i in range(s2t_indices[-1], len(profile2.sequenceP)):
    #         ind2_list.append('-')

    # print(ind1_list)
    # print(ind2_list)

    # aligned_regions1 = []
    # aligned_regions2 = []
    # get_aligned_regions(nucleotide_indices1, aligned_regions1, nucleotide_indices2, aligned_regions2)
    # print(aligned_regions1)
    # print(aligned_regions2)

    aligned_regions1, aligned_regions2 = get_non_gapped_aligned_regions(nucleotide_indices1, nucleotide_indices2, dBG)
    if dBG:
        # print('using new method')
        logger.debug('aligned_regions1: ' + str(aligned_regions1))
        logger.debug('aligned_regions2: ' + str(aligned_regions2))
        # print(aligned_regions1)
        # print(aligned_regions2)
    # get_gapped_index_mapping(aligned_regions1, aligned_regions2)
    # ind1_list, ind2_list, break_points = get_gapped_index_mapping(aligned_regions1, aligned_regions2, profile1, profile2)
    # ind1_list, ind2_list = get_gapped_index_mapping(aligned_regions1, aligned_regions2, profile1, profile2, scores, nucl_scoring_data, dBG)
    ind1_list, ind2_list = get_gapped_index_mapping_v2(aligned_regions1, aligned_regions2, profile1, profile2, scores, nucl_scoring_data, dBG)
    if dBG:
        logger.debug('ind1_list: ' + str(ind1_list))
        logger.debug('ind2_list: ' + str(ind2_list))
        # print(ind1_list)
        # print(ind2_list)

    if len(ind1_list) != len(ind2_list):
        logger.error('Index lengths should be of same length. Need to investigate.')
        # print('Index lengths should be of same length. Need to investigate.')
        sys.exit()

    # aligned_regions1 = []
    # aligned_regions2 = []
    # get_aligned_regions(nucleotide_indices1, aligned_regions1, nucleotide_indices2, aligned_regions2)
    # print(aligned_regions1)
    # print(aligned_regions2)
    
    # return s2t, t2s, ind1_list, ind2_list, break_points
    return s2t, t2s, ind1_list, ind2_list

# def get_nucleotides_in_current_location(ind1, ind2, profile1, profile2):
#     nucl_dict_in_cur_location = {}
#     if ind1 != '-':
#         for hashed_nucl in profile1.sequenceP[ind1]:
#             if hashed_nucl not in nucl_dict_in_cur_location:
#                 nucl_dict_in_cur_location[hashed_nucl] = 0
#             nucl_dict_in_cur_location[hashed_nucl] += profile1.sequenceP[ind1][hashed_nucl]

#     else:
#         hashed_nucl = hash_nucleotide('-')
#         if hashed_nucl not in nucl_dict_in_cur_location:
#             nucl_dict_in_cur_location[hashed_nucl] = 0
#         nucl_dict_in_cur_location[hashed_nucl] += profile1.loopCount

#     if ind2 != '-':
#         for hashed_nucl in profile2.sequenceP[ind2]:
#             if hashed_nucl not in nucl_dict_in_cur_location:
#                 nucl_dict_in_cur_location[hashed_nucl] = 0
#             nucl_dict_in_cur_location[hashed_nucl] += profile2.sequenceP[ind2][hashed_nucl]

#     else:
#         hashed_nucl = hash_nucleotide('-')
#         if hashed_nucl not in nucl_dict_in_cur_location:
#             nucl_dict_in_cur_location[hashed_nucl] = 0
#         nucl_dict_in_cur_location[hashed_nucl] += profile2.loopCount

#     return nucl_dict_in_cur_location

# def get_merged_interactionsP(ind1_list, ind2_list, profile1, profile2):
#     bpsP = {}
#     stksP = {}

#     for i in range(len(ind1_list)):
#         ind1 = ind1_list[i]
#         ind2 = ind2_list[i]

#         if ind1 != '-':


#         if ind2 != '-':


# def find_interactions_from_cur_index(ind, interaction_dict):
#     interactions_from_ind = {}
#     for (s, e) in interaction_dict:
#         if s == ind:
#             if (s, e) not in interactions_from_ind:
#                 interactions_from_ind[(s, e)] = {}
#             for hashed_interaction in interaction_dict[(s, e)]:
#                 if hashed_interaction not in interactions_from_ind[(s, e)]
#                     interactions_from_ind[(s, e)][hashed_interaction] = 0
#                 interactions_from_ind[(s, e)][hashed_interaction] += 

def copy_interactions_as_new_ind(new_interactionsP, interactionsP, ind_list):
    # print('received new_interactionsP')
    # print(new_interactionsP)
    # print('to add')
    # print(interactionsP)
    for (s, e) in interactionsP:
        if s not in ind_list or e not in ind_list:
            print('s or e not in the ind_list. check.')
            sys.exit()

        ns = ind_list.index(s)
        ne = ind_list.index(e)

        if (ns, ne) not in new_interactionsP:
            new_interactionsP[(ns, ne)] = {}
        for hashed_nucl in interactionsP[(s, e)]:
            if hashed_nucl not in new_interactionsP[(ns, ne)]:
                new_interactionsP[(ns, ne)][hashed_nucl] = 0
            new_interactionsP[(ns, ne)][hashed_nucl] += interactionsP[(s, e)][hashed_nucl]

    # print('returning')
    # print(new_interactionsP)

# def generate_merged_profile_based_on_alignment(profile1, profile2, ind1_list, ind2_list, break_points, maximal_clique, interaction_ind_pair_list1, interaction_ind_pair_list2):
def generate_merged_profile_based_on_alignment(profile1, profile2, ind1_list, ind2_list, maximal_clique, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, nucl_scoring_data, dBG):
    sequenceP = []
    bpsP = {}
    stksP = {}
    break_points = []

    # print('max clique size')
    # print(len(maximal_clique))
    # if len(maximal_clique) == 0:
    #     tt = global_align_segments(profile1, profile2, scores, nucl_scoring_data)
    #     print(tt)
    #     sys.exit()

    for i in range(len(ind1_list)):
        ind1 = ind1_list[i]
        ind2 = ind2_list[i]

        nucl_dict_in_cur_location = {}
        if ind1 != '-':
            for hashed_nucl in profile1.sequenceP[ind1]:
                if hashed_nucl not in nucl_dict_in_cur_location:
                    nucl_dict_in_cur_location[hashed_nucl] = 0
                nucl_dict_in_cur_location[hashed_nucl] += profile1.sequenceP[ind1][hashed_nucl]

            # bps = find_interactions_from_cur_index(ind1, profile1.bpsP)
            # stks = find_interactions_from_cur_index(ind1, profile1.stksP)

        else:
            hashed_nucl = hash_nucleotide('-')
            if hashed_nucl not in nucl_dict_in_cur_location:
                nucl_dict_in_cur_location[hashed_nucl] = 0
            nucl_dict_in_cur_location[hashed_nucl] += profile1.loopCount

        if ind2 != '-':
            for hashed_nucl in profile2.sequenceP[ind2]:
                if hashed_nucl not in nucl_dict_in_cur_location:
                    nucl_dict_in_cur_location[hashed_nucl] = 0
                nucl_dict_in_cur_location[hashed_nucl] += profile2.sequenceP[ind2][hashed_nucl]

        else:
            hashed_nucl = hash_nucleotide('-')
            if hashed_nucl not in nucl_dict_in_cur_location:
                nucl_dict_in_cur_location[hashed_nucl] = 0
            nucl_dict_in_cur_location[hashed_nucl] += profile2.loopCount

        # nucl_dict_in_cur_location = get_nucleotides_in_current_location(ind1, ind2, profile1, profile2)
        sequenceP.append(nucl_dict_in_cur_location)

    # bpsP, stksP = get_merged_interactionsP(ind1_list, ind2_list, profile1, profile2)
    # print('initial bpsP')
    # print(bpsP)
    copy_interactions_as_new_ind(bpsP, profile1.bpsP, ind1_list)
    # print('after adding bpsP1')
    # print(bpsP)

    # print('initial stksP')
    # print(stksP)
    copy_interactions_as_new_ind(stksP, profile1.stksP, ind1_list)
    # print('after adding stksP1')
    # print(stksP)

    copy_interactions_as_new_ind(bpsP, profile2.bpsP, ind2_list)
    # print('after adding bpsP2')
    # print(bpsP)
    copy_interactions_as_new_ind(stksP, profile2.stksP, ind2_list)
    # print('after adding stksP2')
    # print(stksP)

    # for (s, e) in profile1.bpsP:
    #     if s not in ind1_list or e not in ind1_list:
    #         print('s or e not in the ind1_list. check.')
    #         sys.exit()

    #     ns = ind1_list.index(s)
    #     ne = ind1_list.index(e)

    #     if (ns, ne) not in bpsP:
    #         bpsP[(ns, ne)] = {}
    #         for hashed_nucl in profile1.bpsP[(s, e)]:
    #             if hashed_nucl not in bpsP[(ns, ne)]:
    #                 bpsP[(ns, ne)][hash_nucl] = 0
    #             bpsP[(ns, ne)][hash_nucl] += profile1.bpsP[(s, e)][hash_nucl]

    if dBG:
        logger.debug('break_points1: ' + str(profile1.break_points))
        logger.debug('break_points2: ' + str(profile2.break_points))
        # print('break_points')
        # print(profile1.break_points)
        # print(profile2.break_points)
    # if len(profile1.break_points) != len(profile2.break_points):
    #     print('Break points mismatch!')
    #     sys.exit()
    # else:
    #     for i in range(len(profile1.break_points)):
    #         b_point1 = profile1.break_points[i]
    #         b_point2 = profile2.break_points[i]
    #         if (b_point1 < len(profile1.sequenceP) and b_point1 not in ind1_list) or \
    #         (b_point2 < len(profile2.sequenceP) and b_point2 not in ind2_list):
    #             print('break point index not found in ind_lists')
    #             sys.exit()
    #         else:
    #             print(ind1_list.index(b_point1))
    #             print(ind2_list.index(b_point2))
    #             if ind1_list.index(b_point1) != ind2_list.index(b_point2):
    #                 print('break point index not matches!')
    #                 sys.exit()
    #             else:
    #                 break_points.append(ind1_list.index(b_point1))



    if len(profile1.break_points) != len(profile2.break_points):
        print('Number of break points mismatch!')
        sys.exit()
    else:
        for i in range(len(profile1.break_points)):
            b_point1 = profile1.break_points[i]
            b_point2 = profile2.break_points[i]

            if (b_point1 < len(profile1.sequenceP) and b_point1 not in ind1_list) or \
            (b_point2 < len(profile2.sequenceP) and b_point2 not in ind2_list):
                print('break point index not found in ind_lists')
                sys.exit()

            elif b_point1 == len(profile1.sequenceP) and b_point2 == len(profile2.sequenceP):
                break_points.append(len(ind1_list))
            else:
                b1a = ind1_list.index(b_point1-1)
                b1b = ind1_list.index(b_point1)
                b2a = ind2_list.index(b_point2-1)
                b2b = ind2_list.index(b_point2)
                if dBG:
                    logger.debug('b1a, b1b: ' + str((b1a, b1b)))
                    logger.debug('b2a, b2b: ' + str((b2a, b2b)))
                    # print('b1a, b1b')
                    # print(b1a, b1b)
                    # print('b2a, b2b')
                    # print(b2a, b2b)
                if b1a == b2a or b1b == b2b:
                    # still okay
                    # print('okay')
                    break_points.append(ind1_list.index(b_point1))
                elif b2a <= b1b and b1b <= b2b:
                    break_points.append(b1b)
                elif b1a <= b2b and b2b <= b1b:
                    break_points.append(b2b)
                else:
                    # print('problem')
                    print('break point index not matches! please check')
                    sys.exit()

                # if b1 != b2:
                #     print('break point index not matches!')
                #     print('prev ind_list')
                #     print(ind1_list)
                #     print(ind2_list)
                #     max_b = max(b1, b2)
                #     if max_b == b1:
                #         ind2_list = ind2_list[:b2+1] + ['-' for k in range(b2+1,b1+1)] + ind2_list[b2+1:]
                #     else:
                #         #max_b == b2
                #         ind1_list = ind1_list[:b1+1] + ['-' for k in range(b1+1,b2+1)] + ind1_list[b1+1:]

                #     break_points.append(max_b)
                #     print('updated ind_list')
                #     print(ind1_list)
                #     print(ind2_list)
                #     sys.exit()
                # else:
                #     break_points.append(ind1_list.index(b_point1))



    # print(sequenceP)
    # print(bpsP)
    # print(stksP)
    if dBG:
        logger.debug('break_pointss: ' + str(break_points))
        # print(break_points)
    # print('adding ' + str(profile1.loopCount) + ' and ' + str(profile2.loopCount))
    return Profile(sequenceP, bpsP, stksP, break_points, profile1.loopCount + profile2.loopCount, "merged")

def csv_to_list(lines):
    list_of_lists = []
    for line in lines:
        pieces = line.strip().split(',')
        list_of_lists.append(list(map(lambda x: x.strip(), pieces)))

    return list_of_lists

def strToNode(loop_str):

    (chain, region) = loop_str.split(':')
    segments = region.split('_')
    node = Node(chain, segments)

    return node

def translate_interaction_indices(interactionsP, ind_map, interaction_type):
    new_interactionsP = {}
    for (s, e) in interactionsP:
        ns = ind_map[s]
        ne = ind_map[e]
        reverse_flag = False
        if ns > ne:
            reverse_flag = True
            t = ns
            ns = ne
            ne = t
        if (ns, ne) not in new_interactionsP:
            new_interactionsP[(ns, ne)] = {}
        for hashed_interaction in interactionsP[(s, e)]:
            new_hashed_interaction = None
            if reverse_flag == True:
                if interaction_type == 'bp':
                    interaction, ntd_pair = interpret_hash_value(hashed_interaction)
                    new_interaction = get_reversed_interaction(interaction, 'bp')
                    new_ntd_pair = ntd_pair[1] + ntd_pair[0]
                    new_hashed_interaction = hash_basepair(new_interaction, new_ntd_pair)
                elif interaction_type == 'stack':
                    interaction, ntd_pair = interpret_hash_value(hashed_interaction)
                    new_interaction = get_reversed_interaction(interaction, 'stack')
                    new_hashed_interaction = hash_stacking(new_interaction, '')
                else:
                    print('Invalid interaction type')
                    sys.exit()

            else:
                new_hashed_interaction = hashed_interaction

            if new_hashed_interaction != None and new_hashed_interaction not in new_interactionsP[(ns, ne)]:
                new_interactionsP[(ns, ne)][new_hashed_interaction] = 0
            new_interactionsP[(ns, ne)][new_hashed_interaction] += interactionsP[(s, e)][hashed_interaction]

    return new_interactionsP

# def translate_break_points(break_points, ind_map):
#     new_break_points = []
#     for break_point in break_points[:-1]:
#         new_break_points.append(ind_map[break_point])
#     new_break_points.append(break_points[-1])

#     return new_break_points

def translate_break_points(break_points, ind_map):
    new_break_points = []
    if ind_map[0] != 0:
        new_break_points.append(ind_map[0])
    for break_point in break_points[:-1]:
        if ind_map[break_point] != 0:
            new_break_points.append(ind_map[break_point])
    new_break_points.append(break_points[-1])

    return sorted(new_break_points)

def generate_all_profile_combinations(profile, dBG):
    if dBG:
        logger.debug('Processing profile:')
        # print('Processing profile:')
        print_formatted_profile(profile)

    profile_combinations = []
    # profile_combinations.append(profile)

    regions = []
    s = 0
    if dBG:
        logger.debug('break_points: ' + str(profile.break_points))
        # print(profile.break_points)
    for break_point in profile.break_points[:-1]:
        e = break_point - 1
        regions.append((s, e))
        s = break_point
    e = profile.break_points[-1] - 1
    regions.append((s, e))
    for i in range(len(regions)):
        region_i = regions[-i:] + regions[:-i]
        if dBG:
            logger.debug('region_i: ' + str(region_i))
            # print(region_i)
        seqP = []
        ind_map = {}
        ind = 0
        for r in region_i:
            s, e = r
            # seqP += profile.sequenceP[s, e+1]
            for k in range(s, e+1):
                seqP.append(profile.sequenceP[k])
                ind_map[k] = ind
                ind += 1

        if dBG:
            logger.debug('ind_map: ' + str(ind_map))
            # print(ind_map)
        # sys.exit()
        bpsP = translate_interaction_indices(profile.bpsP, ind_map, 'bp')
        stksP = translate_interaction_indices(profile.stksP, ind_map, 'stack')
        break_points = translate_break_points(profile.break_points, ind_map)
        profile_combinations.append(Profile(seqP, bpsP, stksP, break_points, profile.loopCount, profile.loopName))

    # sys.exit()

    # pdb_chain, regions = loop.strip().split(':')
    # regions = regions.strip().split('_')
    # loop = []
    # for region in regions:
    #     s, e = region.strip().split('-')
    #     loop.append((s, e))
    
    # for i in range(len(loop)):
    #     loop_i = rotate(loop, i)
    #     loop_combinations.append(pdb_chain + ':' + '_'.join(list(map(lambda x: '-'.join(x), loop_i))))

    if dBG:
        # print combinations
        for i, profile in enumerate(profile_combinations):
            logger.debug('Combination ' + str(i+1))
            # print('Combination ' + str(i+1))
            print_formatted_profile(profile)

    # sys.exit()

    return profile_combinations

def revise_profile_based_on_gap(profile, threshold=80.0):
    # print('revising')
    # print(profile.loopCount)
    candidates = []
    for i, nucl_dict in enumerate(profile.sequenceP):
        hashed_nucl, freq = list(sorted(nucl_dict.items(), key=lambda item: item[1], reverse=True))[0]
        # print(hashed_nucl, freq)
        if reverse_hash_nucleotide(hashed_nucl) == '-':
            p = freq * 100.0 / profile.loopCount
            # print(i, p)
            if p >= threshold:
                candidates.append(i)

    # print(candidates)
    for (s, e) in profile.bpsP:
        bpsP = profile.bpsP[(s, e)]
        max_percentage = 0.0
        for hashed_interaction in bpsP:
            percentage = bpsP[hashed_interaction] / profile.loopCount
            max_percentage = max(percentage, max_percentage)
        if s in candidates:
            if len(bpsP) > 1 or (len(bpsP) == 1 and max_percentage > 0.25):
                candidates.remove(s)
        if e in candidates:
            if len(bpsP) > 1 or (len(bpsP) == 1 and max_percentage > 0.25):
                candidates.remove(e)

    # print(candidates)
    if len(candidates) > 0:
        # print('profile before revising:')
        # print_formatted_profile(profile)

        ind_map = {}
        for i in range(len(profile.sequenceP)):
            if i in candidates:
                continue
            ind_map[i] = i - len(list(filter(lambda x: x < i, candidates)))

        ind_map[len(profile.sequenceP)] = len(profile.sequenceP) - len(candidates)

        # print(ind_map)

        sequenceP = []
        bpsP = {}
        stksP = {}
        for i, nucl_dict in enumerate(profile.sequenceP):
            if i not in candidates:
                sequenceP.append(nucl_dict)

        for (s, e) in profile.bpsP:
            if s not in candidates and e not in candidates:
                ns = ind_map[s]
                ne = ind_map[e]
                bpsP[(ns, ne)] = profile.bpsP[(s, e)]

        for (s, e) in profile.stksP:
            if s not in candidates and e not in candidates:
                ns = ind_map[s]
                ne = ind_map[e]
                stksP[(ns, ne)] = profile.stksP[(s, e)]

        break_points = []
        for i in profile.break_points:
            if i in candidates:
                break_points.append(i)
            else:
                break_points.append(ind_map[i])

        new_profile = Profile(sequenceP, bpsP, stksP, break_points, profile.loopCount, profile.loopName)
        # print('new profile after revising:')
        # print_formatted_profile(new_profile)
        # sys.exit()
        return new_profile

    return profile

def print_formatted_profile(profile, dBG=False):
    # sequenceP, bpsP, stksP, break_points, loopCount
    if dBG:
        logger.debug('\nSequence:\n-----------------------------------')
    else:
        print('\nSequence:\n-----------------------------------')
    # print(profile)
    for i, nucl_dict in enumerate(profile.sequenceP):
        if dBG:
            logger.debug(str(i) + '\t')
        else:
            print(i, end='\t')
        nucl_dict = dict(sorted(nucl_dict.items(), key=lambda item: item[1], reverse=True))
        for j, hashed_nucl in enumerate(nucl_dict):
            if j == len(nucl_dict) - 1:
                if dBG:
                    logger.debug(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]))
                else:
                    print(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]))
            else:
                if dBG:
                    logger.debug(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]) + ', ')
                else:
                    print(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]), end=', ')

    if dBG:
        logger.debug('\nBase-pairs:\n-----------------------------------')
    else:
        print('\nBase-pairs:\n-----------------------------------')
    for (s, e) in profile.bpsP:
        if dBG:
            logger.debug(str(s) + '-' + str(e) + '\t')
        else:
            print(str(s) + '-' + str(e), end='\t')

        profile.bpsP[(s, e)] = dict(sorted(profile.bpsP[(s, e)].items(), key=lambda item: item[1], reverse=True))
        for j, hashed_interaction in enumerate(profile.bpsP[(s, e)]):
            interaction, ntd_pair = interpret_hash_value(hashed_interaction)
            if j == len(profile.bpsP[(s, e)]) - 1:
                if dBG:
                    logger.debug(interaction + ' - ' + ntd_pair + ': ' + str(profile.bpsP[(s, e)][hashed_interaction]))
                else:
                    print(interaction + ' - ' + ntd_pair + ': ' + str(profile.bpsP[(s, e)][hashed_interaction]))
            else:
                if dBG:
                    logger.debug(interaction + ' - ' + ntd_pair + ': ' + str(profile.bpsP[(s, e)][hashed_interaction]) + ', ')
                else:
                    print(interaction + ' - ' + ntd_pair + ': ' + str(profile.bpsP[(s, e)][hashed_interaction]), end=', ')

    if dBG:
        logger.debug('\nBase-stackings:\n-----------------------------------')
    else:
        print('\nBase-stackings:\n-----------------------------------')
    for (s, e) in profile.stksP:
        if dBG:
            logger.debug(str(s) + '-' + str(e) + '\t')
        else:
            print(str(s) + '-' + str(e), end='\t')

        profile.stksP[(s, e)] = dict(sorted(profile.stksP[(s, e)].items(), key=lambda item: item[1], reverse=True))
        for j, hashed_interaction in enumerate(profile.stksP[(s, e)]):
            interaction, ntd_pair = interpret_hash_value(hashed_interaction)
            if j == len(profile.stksP[(s, e)]) - 1:
                if dBG:
                    logger.debug(interaction + ': ' + str(profile.stksP[(s, e)][hashed_interaction]))
                else:
                    print(interaction + ': ' + str(profile.stksP[(s, e)][hashed_interaction]))
            else:
                if dBG:
                    logger.debug(interaction + ': ' + str(profile.stksP[(s, e)][hashed_interaction]) + ', ')
                else:
                    print(interaction + ': ' + str(profile.stksP[(s, e)][hashed_interaction]), end=', ')

    if dBG:
        logger.debug('\nBreak points:\n-----------------------------------')
        logger.debug(', '.join(list(map(str, profile.break_points))))
        logger.debug('\nTotal loops: ' + str(profile.loopCount))
    else:
        print('\nBreak points:\n-----------------------------------')
        print(', '.join(list(map(str, profile.break_points))))
        print('\nTotal loops: ' + str(profile.loopCount))

def write_formatted_profile_to_file(profile, output_dir, output_fname):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fp = open(os.path.join(output_dir, output_fname), 'w')
    
    fp.write('#info=sequence\n')
    for i, nucl_dict in enumerate(profile.sequenceP):
        # fp.write(str(i) + '\t')
        # print(i, end='\t')
        nucl_dict = dict(sorted(nucl_dict.items(), key=lambda item: item[1], reverse=True))
        for j, hashed_nucl in enumerate(nucl_dict):
            if j == len(nucl_dict) - 1:
                fp.write(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]) + '\n')
            else:
                fp.write(reverse_hash_nucleotide(hashed_nucl) + ': ' + str(nucl_dict[hashed_nucl]) + ', ')

    fp.write('#info=basepair\n')
    for (s, e) in profile.bpsP:
        fp.write(str(s) + '-' + str(e) + '\t')

        profile.bpsP[(s, e)] = dict(sorted(profile.bpsP[(s, e)].items(), key=lambda item: item[1], reverse=True))
        for j, hashed_interaction in enumerate(profile.bpsP[(s, e)]):
            interaction, ntd_pair = interpret_hash_value(hashed_interaction)
            if j == len(profile.bpsP[(s, e)]) - 1:
                fp.write(interaction + ' - ' + ntd_pair + ': ' + str(profile.bpsP[(s, e)][hashed_interaction]) + '\n')
            else:
                fp.write(interaction + ' - ' + ntd_pair + ': ' + str(profile.bpsP[(s, e)][hashed_interaction]) + ', ')

    fp.write('#info=stacking\n')
    for (s, e) in profile.stksP:
        fp.write(str(s) + '-' + str(e) + '\t')

        profile.stksP[(s, e)] = dict(sorted(profile.stksP[(s, e)].items(), key=lambda item: item[1], reverse=True))
        for j, hashed_interaction in enumerate(profile.stksP[(s, e)]):
            interaction, ntd_pair = interpret_hash_value(hashed_interaction)
            if j == len(profile.stksP[(s, e)]) - 1:
                fp.write(interaction + ': ' + str(profile.stksP[(s, e)][hashed_interaction]) + '\n')
            else:
                fp.write(interaction + ': ' + str(profile.stksP[(s, e)][hashed_interaction]) + ', ')

    fp.write('#info=break_points\n')
    fp.write(', '.join(list(map(str, profile.break_points))) + '\n')

    fp.write('#info=total_loops\n')
    fp.write(str(profile.loopCount) + '\n')

    fp.close()

def rename_filename_if_exists(filename):
    if os.path.exists(filename):
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name, ext = os.path.splitext(filename)
        new_filename = base_name + '_bck_' + timestamp + ext
        os.rename(filename, new_filename)

def create_directory(dir_to_create):
    if not os.path.exists(dir_to_create):
        os.makedirs(dir_to_create)

def delete_directory(dir_to_delete):
    if os.path.exists(dir_to_delete):
        shutil.rmtree(dir_to_delete)

def remove_all_from_dir(mypath, exception_dirs=[]):
    for root, dirs, files in os.walk(mypath):
        # for file in files:  #files in all dirs
        #     if file not in exception_files:
        #         os.remove(os.path.join(root, file))
        for dir in dirs:
            if dir not in exception_dirs:
                shutil.rmtree(os.path.join(root, dir))

def create_required_directories(directories):

    create_directory(directories.pdbx_dir)
    create_directory(directories.fasta_dir)
    create_directory(directories.loop_dir)
    create_directory(directories.pdb_fasta_mapping_dir)
    create_directory(directories.annotation_dir)
    create_directory(directories.temp_dir)
    create_directory(directories.pickles_dir)
    create_directory(directories.partial_pdb_dir)

def wait_for_certain_time_according_to_wait_factor(n, wait_factor=0.1, max_wait_time=600):
    wait_time = n * wait_factor
    wait_time = min(wait_time, max_wait_time)
    time.sleep(wait_time)

# This was done to make the code compatible with any case-insensitive OS
def get_modified_chain_id_if_any_lowercase_letter(chain_id):
    if chain_id != chain_id.upper():
        chain_id = chain_id + '_'
    return chain_id

def get_way_counts(loops):
    way_counts = []
    for loop in loops:
        pdb_chain, regions = loop.strip().split(':')
        regions = regions.strip().split('_')
        way_counts.append(len(regions))

    return list(set(way_counts))

def get_zscore(x, m, std):
    if std == 0:
        return 0
    return (x-m)/std

def get_less_significant_filtered_interactions(interactionsP, total_value, threshold_value):
    new_interactionsP = {}
    for (s, e) in interactionsP:
        for hashed_interaction in interactionsP[(s, e)]:
            percentage = interactionsP[(s, e)][hashed_interaction] / total_value
            if percentage >= threshold_value:
                if (s, e) not in new_interactionsP:
                    new_interactionsP[(s, e)] = {}
                if hashed_interaction not in new_interactionsP[(s, e)]:
                    new_interactionsP[(s, e)][hashed_interaction] = 0
                new_interactionsP[(s, e)][hashed_interaction] = interactionsP[(s, e)][hashed_interaction]
    return new_interactionsP

def filter_less_significant_nucl_and_interactions_from_a_profile(profile, threshold_value = 0.1):
    new_sequenceP = []
    # print(profile.sequenceP)
    for i, nucl_dict in enumerate(profile.sequenceP):
        new_nucl_dict = {}
        for hashed_nucl in nucl_dict:
            percentage = nucl_dict[hashed_nucl] / profile.loopCount
            # print(nucl_dict[hashed_nucl], profile.loopCount, percentage)
            if percentage >= threshold_value:
                if hashed_nucl not in new_nucl_dict:
                    new_nucl_dict[hashed_nucl] = nucl_dict[hashed_nucl]
        new_sequenceP.append(new_nucl_dict)
        if len(new_nucl_dict) == 0:
            # print(i)
            logger.warning('Need to delete nucl which is unexpected. Please check for any errors.')
            sys.exit()

    new_bpsP = get_less_significant_filtered_interactions(profile.bpsP, profile.loopCount, threshold_value*0.75)
    new_stksP = get_less_significant_filtered_interactions(profile.stksP, profile.loopCount, threshold_value)
    return Profile(new_sequenceP, new_bpsP, new_stksP, profile.break_points, profile.loopCount, profile.loopName)

def get_loop_length(loop):
    pdb_chain, regions = loop.strip().split(':')
    regions = regions.strip().split('_')
    loop_len = 0
    junc_cnt = len(regions)
    for region in regions:
        s, e = region.strip().split('-')
        s = int(s)
        e = int(e)
        loop_len += (e - s + 1)

    return loop_len, junc_cnt

def print_loop_intera_stats(families, directories):
    bps_cnt = 0
    stks_cnt = 0
    loops_cnt = 0
    for family_id in families:
        loops = families[family_id]
        loops_cnt += len(loops)
        for loop in loops:
            if os.path.exists(os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf')):
                _, _, bp_cnt, _, stk_cnt, _ = load_loop_data(os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf'), False)
            else:
                loops_cnt -= 1
                continue
            bps_cnt += bp_cnt
            stks_cnt += stk_cnt
    avg_bps_cnt = bps_cnt / loops_cnt
    avg_stks_cnt = stks_cnt / loops_cnt

    print(avg_bps_cnt)
    print(avg_stks_cnt)

def print_input_stat(families):
    sum_sum_motif_len = 0
    max_max_motif_len = -1
    min_min_motif_len = 999
    total_loop_count = 0
    junc_dict = {}
    junc_wise_loops = {}
    for family_id in families:
        print(family_id, end='\t')
        loops = families[family_id]
        total_loop_count += len(loops)
        sum_motif_len = 0
        max_motif_len = -1
        min_motif_len = 999
        for loop in loops:
            motif_len, junc_cnt = get_loop_length(loop)
            if junc_cnt not in junc_dict:
                junc_dict[junc_cnt] = []
                junc_wise_loops[junc_cnt] = []
            junc_dict[junc_cnt].append(junc_cnt)
            junc_wise_loops[junc_cnt].append(loop)
            sum_motif_len += motif_len
            min_motif_len = min(min_motif_len, motif_len)
            max_motif_len = max(max_motif_len, motif_len)
        sum_sum_motif_len += sum_motif_len
        min_min_motif_len = min(min_min_motif_len, min_motif_len)
        max_max_motif_len = max(max_max_motif_len, max_motif_len)

        avg_motif_len = int(round(sum_motif_len / len(loops), 0))
        print(str(len(loops)) + '\t' + str(avg_motif_len) + ' (' + str(min_motif_len) + '-' + str(max_motif_len) + ')')
    avg_avg_motif_len = int(round(sum_sum_motif_len / total_loop_count, 0))
    print('Overall', end='\t')
    print(str(total_loop_count) + '\t' + str(avg_avg_motif_len) + ' (' + str(min_min_motif_len) + '-' + str(max_max_motif_len) + ')', end='\t')
    # print(str(junc_dict))
    for junc_cnt in junc_dict:
        print(str(junc_cnt) + ' - ' + str(len(junc_dict[junc_cnt])))

        loops = junc_wise_loops[junc_cnt]
        sum_motif_len = 0
        max_motif_len = -1
        min_motif_len = 999
        for loop in loops:
            motif_len, junc_cnt = get_loop_length(loop)
            # if junc_cnt not in junc_dict:
            #     junc_dict[junc_cnt] = []
            # junc_dict[junc_cnt].append(junc_cnt)
            sum_motif_len += motif_len
            min_motif_len = min(min_motif_len, motif_len)
            max_motif_len = max(max_motif_len, motif_len)
        # sum_sum_motif_len += sum_motif_len
        # min_min_motif_len = min(min_min_motif_len, min_motif_len)
        # max_max_motif_len = max(max_max_motif_len, max_motif_len)
        avg_motif_len = int(round(sum_motif_len / len(loops), 0))
        print(str(len(loops)) + '\t' + str(avg_motif_len) + ' (' + str(min_motif_len) + '-' + str(max_motif_len) + ')')