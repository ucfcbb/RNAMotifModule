def parse_block(block):
	total_lines = len(block)
	if block[5].strip() == '***** No hits found *****':
		# print('No hit')
		pass
	elif block[4].strip().startswith('Sequences producing significant alignments:'):
		print('Hit with ' + block[0].strip().split('=')[1].strip())
		line_no = 19
		if block[line_no].startswith('Sbjct'):
			# print(block[line_no])
			pieces = block[line_no].strip().split()
			# print(pieces)
			end = pieces[1]
			start = -1
			while line_no < total_lines and block[line_no].startswith('Sbjct'):
				# print(block[line_no])
				pieces = block[line_no].strip().split()
				# print(pieces)
				start = pieces[3]
				line_no += 4

			print('FASTA locations: ' + start + ' - ' + end + '\n')

		else:
			print('Check block!')
			print('**************************************')
			for line in block:
				print(line)
			print('**************************************')


def main():
	# alignment_fname = 'CLIP-seq_test/test_z/GSE44959_RAW/dbr1_alignment_results_4V88_A6.txt'

	# alignment_fname = '../MotifModule_Explore3/CLIP-seq_test/repeat2/GSE44959_RAW/dbr1_rep1_alignment_results_4V88_A6.txt'
	# alignment_fname = '../MotifModule_Explore3/CLIP-seq_test/repeat2/GSE44959_RAW/dbr1_rep2_alignment_results_4V88_A6.txt'
	# alignment_fname = '../MotifModule_Explore3/CLIP-seq_test/repeat2/GSE44959_RAW/drn1_rep1_alignment_results_4V88_A6.txt'
	# alignment_fname = '../MotifModule_Explore3/CLIP-seq_test/repeat2/GSE44959_RAW/drn1_rep2_alignment_results_4V88_A6.txt'

	# alignment_fname = '../MotifModule_Explore3/CLIP-seq_test/repeat2/GSE44959_RAW/dbr1_rep1_alignment_results_5TBW_1.txt'
	# alignment_fname = '../MotifModule_Explore3/CLIP-seq_test/repeat2/GSE44959_RAW/dbr1_rep2_alignment_results_5TBW_1.txt'
	# alignment_fname = '../MotifModule_Explore3/CLIP-seq_test/repeat2/GSE44959_RAW/drn1_rep1_alignment_results_5TBW_1.txt'
	alignment_fname = '../MotifModule_Explore3/CLIP-seq_test/repeat2/GSE44959_RAW/drn1_rep2_alignment_results_5TBW_1.txt'

	fp = open(alignment_fname, 'r')
	lines = fp.readlines()
	fp.close()

	total_lines = len(lines)
	i = 0
	while i < total_lines:
		if lines[i].startswith('Query= chr'):
			block = []
			while not lines[i].startswith('Effective search space used:') and i < total_lines:
				block.append(lines[i].strip())
				i += 1
			block.append(lines[i].strip())
			i += 1
			parse_block(block)
		i += 1


if __name__ == '__main__':
	main()
