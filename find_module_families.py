import sys
import os
import re

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

def load_pdb_organism_data(lib_dir):
	# lib_dir = 'my_lib'
	pdb_organism_details = read_pdb_chain_organism_details(os.path.join(lib_dir, 'PDB_Chain_Organism_Details.tsv'))
	pdb_organism_details_scrapped = read_pdb_chain_organism_details(os.path.join(lib_dir, 'PDB_Chain_Organism_Details_scrapped.tsv'))

	for pdb_chain in pdb_organism_details_scrapped:
		if pdb_chain not in pdb_organism_details:
			pdb_organism_details[pdb_chain] = pdb_organism_details_scrapped[pdb_chain]

	return pdb_organism_details

def main():

	lib_dir = 'my_lib'
	pdb_organism_details = load_pdb_organism_data(lib_dir)


	fp = open('input/known_families_HL.csv')
	lines = fp.readlines()
	fp.close()
	fp = open('input/known_families_IL.csv')
	lines += fp.readlines()
	fp.close()

	fp = open('input_rna_fam_strat.tsv', 'w')
	fp.write('')
	fp.close()

	for line in lines:
		pieces = line.strip().split(',')
		motif_fam_id = pieces[0]
		rna_families = {}
		unknown_pdb_chains = []
		for motif in pieces[1:]:
			pdb_chain = motif.strip().split(':')[0]
			rna_fam_name = 'N/A'
			if pdb_chain in pdb_organism_details:
				rna_fam_name = pdb_organism_details[pdb_chain][0]
			if rna_fam_name not in rna_families:
				rna_families[rna_fam_name] = 0
			rna_families[rna_fam_name] += 1

		formatted = ', '.join(f'"{k}": {v}' for k, v in rna_families.items() if k != 'N/A')
		if len(unknown_pdb_chains) > 0:
			formatted += ', unknown: ' + str(len(unknown_pdb_chains))
		print(motif_fam_id + '\t' + formatted)
		fp = open('input_rna_fam_strat.tsv', 'a')
		fp.write(motif_fam_id + '\t' + formatted + '\n')
		fp.close()

	# sys.exit()

	fp = open('output/identified_motif_modules.txt')
	lines = fp.readlines()
	fp.close()

	
	# print(pdb_organism_details)

	# sys.exit()
	fp = open('rna_family_stratification.tsv', 'w')
	fp.write('')
	fp.close()

	for line in lines[1:]:
		pieces = line.strip().split('\t')
		module = pieces[0]
		instances = pieces[2].strip()
		# print(instances)
		instances = [tuple(group.split(',')) for group in re.findall(r'\((.*?)\)', instances)]
		# print(instances)
		rna_families = {}
		unknown_pdb_chains = []
		for instance in instances:
			# instance.strip('()')
			pdb_chain = instance[0].split(':')[0]
			rna_fam_name = 'N/A'
			if pdb_chain in pdb_organism_details:
				rna_fam_name = pdb_organism_details[pdb_chain][0]
			if rna_fam_name not in rna_families:
				rna_families[rna_fam_name] = 0
			rna_families[rna_fam_name] += 1

			if rna_fam_name == 'N/A':
				unknown_pdb_chains.append(pdb_chain)

		formatted = ", ".join(f"{k}: {v}" for k, v in rna_families.items() if k != 'N/A')
		if len(unknown_pdb_chains) > 0:
			formatted += ', not found: ' + ','.join(unknown_pdb_chains)
		print(module + '\t' + formatted)
		fp = open('rna_family_stratification.tsv', 'a')
		fp.write(module + '\t' + formatted + '\n')
		fp.close()

if __name__ == '__main__':
	main()