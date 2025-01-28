import os
import sys

from classes import *

def convert_a_cluster_from_FASTA_to_PDB(families, directories):
    families_pdb = {}
    for family_id in families:
        families_pdb[family_id] = []
        loops = families[family_id]
        for loop in loops:
            loop_pdb = convert_a_loop_from_FASTA_to_PDB(loop, directories)
            families_pdb[family_id].append(loop_pdb)
    return families_pdb

def convert_a_cluster_from_PDB_to_FASTA(families, directories):
    families_pdb = {}
    for family_id in families:
        families_pdb[family_id] = []
        loops = families[family_id]
        for loop in loops:
            loop_pdb = convert_a_loop_from_PDB_to_FASTA(loop, directories)
            families_pdb[family_id].append(loop_pdb)
    return families_pdb

def convert_a_loop_from_PDB_to_FASTA(loop, directories):
    # print(loop)
    pdb_chain, segments = loop.strip().split(':')
    pdb_id, chain_id = pdb_chain.strip().split('_')

    mapping_file_name = pdb_chain + '.rmsx.nch'

    # This was done to make the code compatible with any case-insensitive OS
    if chain_id != chain_id.upper():
        mapping_file_name = pdb_chain + '_.rmsx.nch'

    converter = PDB_FASTA_Index_Converter(directories.pdb_fasta_mapping_dir, mapping_file_name)

    segments = segments.strip().split('_')
    converted_segments = []
    for segment in segments:
        if segment.count('-') == 1:
            a, b = segment.strip().split('-')
        # elif segment.count('-') == 2:
        #     pcs = segment.strip().split('-')
        #     a, b = '-'.join(pcs[:1]), '-'.join(pcs[1:])
        else:
            pcs = segment.strip().split('-')
            a, b = '-'.join(pcs[:2]), '-'.join(pcs[2:])

        icode_a = ''
        icode_b = ''
        if '.' in a:
            a, icode_a = a.strip().split('.')
        if '.' in b:
            b, icode_b = b.strip().split('.')

        a_pdb = Chainindex(chain_id, int(a), icode_a)
        b_pdb = Chainindex(chain_id, int(b), icode_b)

        a_fasta = converter.convert_PDBindx_To_FASTAindx(a_pdb)
        b_fasta = converter.convert_PDBindx_To_FASTAindx(b_pdb)

        converted_segments.append(str(a_fasta) + '-' + str(b_fasta))

    converted_segments = '_'.join(converted_segments)
    converted_loop = pdb_chain + ':' + converted_segments

    return converted_loop

def convert_a_loop_from_FASTA_to_PDB(loop, directories):
    pdb_chain, segments = loop.strip().split(':')
    pdb_id, chain_id = pdb_chain.strip().split('_')

    mapping_file_name = pdb_chain + '.rmsx.nch'

    # This was done to make the code compatible with any case-insensitive OS
    if chain_id != chain_id.upper():
        mapping_file_name = pdb_chain + '_.rmsx.nch'

    converter = PDB_FASTA_Index_Converter(directories.pdb_fasta_mapping_dir, mapping_file_name)

    segments = segments.strip().split('_')
    converted_segments = []
    for segment in segments:
        a, b = segment.strip().split('-')

        a_fasta = a
        b_fasta = b

        a_pdb = converter.convert_FASTAindx_To_PDBindx(a_fasta)
        b_pdb = converter.convert_FASTAindx_To_PDBindx(b_fasta)

        if len(a_pdb.icode) == 0:
            a_pdb = str(a_pdb.seqnum)
        else:
            a_pdb = str(a_pdb.seqnum) + '.' + str(a_pdb.icode)

        if len(b_pdb.icode) == 0:
            b_pdb = str(b_pdb.seqnum)
        else:
            b_pdb = str(b_pdb.seqnum) + '.' + str(b_pdb.icode)

        converted_segments.append(a_pdb + '-' + b_pdb)
        
    converted_segments = '_'.join(converted_segments)
    converted_loop = pdb_chain + ':' + converted_segments

    return converted_loop
