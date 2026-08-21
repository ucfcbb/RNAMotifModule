import os
import sys
import re
import random

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from classes import *
from config import *
from converter import *
from utils import *

def assign_peak_tracks(peak_regions):
    peak_regions = sorted(peak_regions, key=lambda x: x[0])
    tracks = []
    for start, end in peak_regions:
        placed = False
        for track in tracks:
            # check against the last peak placed in this row
            last_start, last_end = track[-1]
            # no overlap
            if start > last_end:
                track.append((start, end))
                placed = True
                break
        if not placed:
            tracks.append([(start, end)])
    return tracks

def convert_module_dict(module_dict):
    module_instances = []
    for module_name, instances in module_dict.items():
        for instance in instances:
            motifs = []
            for motif in instance:
                pdb_chain, regions = motif.split(":")
                segments = []
                for region in regions.split("_"):
                    start, end = region.split("-")
                    segments.append((int(start), int(end)))
                motifs.append(segments)
            module_instances.append((module_name, motifs))
    return module_instances

def get_module_shortcode(module_name):
	module_name = tuple(module_name.strip("()").split(","))
	fam_id_list = []
	for family_id in module_name:
		# print(family_id)
		family_id = '_'.join(family_id.strip().split('_')[1:])
		# print(family_id)
		fam_id_list.append(get_known_motif_shortcode(family_id))

	# print('-'.join(fam_id_list))
	# sys.exit()
	return '-'.join(fam_id_list)


def main():
	directories = get_base_dir_names()

	target_pdb_chain = '4V88_A6'
	RNA_LENGTH = 1800

	# target_pdb_chain = '5TBW_1'
	# RNA_LENGTH = 3149

	# loading peak locations
	fp = open(target_pdb_chain + '_peak_FASTA_locations.txt')
	lines = fp.readlines()
	fp.close()

	peak_regions = []
	for line in lines:
		if line.startswith('FASTA'):
			start, end = line.strip().split(':')[1].strip().split('-')
			start = int(start.strip())
			end = int(end.strip())
			peak_regions.append((start, end))

	# print(len(peak_regions))
	# sys.exit()
	fp = open('output/identified_motif_modules.txt')
	lines = fp.readlines()
	fp.close()


	target_instances = []
	target_instances_with_module_names = {}

	for line in lines[1:]:
		# print(line)
		pieces = line.strip().split('\t')
		module = pieces[0]
		module = get_module_shortcode(module)
		instances = pieces[2].strip()
		# print(instances)
		instances = [tuple(group.split(',')) for group in re.findall(r'\((.*?)\)', instances)]

		selected_instances = []
		for instance in instances:
			pdb_chain = instance[0].split(':')[0]
			if pdb_chain == target_pdb_chain:
				target_instances.append(instance)
				selected_instances.append(instance)

		if len(selected_instances) > 0:
			target_instances_with_module_names[module] = selected_instances

	target_instances = [tuple(map(lambda x: convert_a_loop_from_PDB_to_FASTA(x, directories), instance)) for instance in target_instances]
	module_instances = convert_module_dict(target_instances_with_module_names)

	peak_tracks = assign_peak_tracks(peak_regions)
	# fig_height = max(6, len(module_instances)*0.4 + len(peak_tracks)*0.4 + 2)
	fig_height = max(2.5, len(module_instances)*0.20 + len(peak_tracks)*0.15 + 0.5)

	fig, ax = plt.subplots(figsize=(10, fig_height))
	PEAK_TRACK_SPACING = 0.3
	MODULE_TRACK_SPACING = 0.3
	PEAK_MODULE_GAP = 0.5

	RNA_LINEWIDTH = 5
	PEAK_LINEWIDTH = 3
	MODULE_LINEWIDTH = 2

	current_y = 0

	ax.plot(
		[1, RNA_LENGTH],
		[current_y, current_y],
		color='black',
		linewidth=RNA_LINEWIDTH
	)

	ax.text(
		-40,
		current_y,
		"18S rRNA",
		ha='right',
		va='center',
		fontsize=10
	)

	peak_color = '#4C78A8'

	for track_index, track in enumerate(peak_tracks):
		y = -(track_index + 1) * PEAK_TRACK_SPACING

		ax.text(
			-40,
			y,
			f"Peak Track {track_index + 1}",
			ha='right',
			va='center',
			fontsize=9
		)

		for start, end in track:

			ax.plot(
				[start, end],
				[y, y],
				color=peak_color,
				linewidth=PEAK_LINEWIDTH,
				solid_capstyle='butt'
			)

	last_peak_y = -(len(peak_tracks)) * PEAK_TRACK_SPACING
	module_start_y = last_peak_y - PEAK_MODULE_GAP
	module_color = '#C0392B'

	for idx, (module_name, motifs) in enumerate(module_instances):
		# y = module_start_y - idx
		y = module_start_y - idx * MODULE_TRACK_SPACING

		ax.text(
			-40,
			y,
			module_name,
			ha='right',
			va='center',
			fontsize=9
		)

		for motif in motifs:
			for start, end in motif:
				ax.plot(
					[start, end],
					[y, y],
					color=module_color,
					linewidth=MODULE_LINEWIDTH,
					solid_capstyle='butt'
				)

	ax.set_xlim(0, RNA_LENGTH)
	# ax.set_ylim(
	# 	module_start_y - len(module_instances) - 1,
	# 	1
	# )
	lowest_module_y = module_start_y - (len(module_instances)-1) * MODULE_TRACK_SPACING
	ax.set_ylim(lowest_module_y - 0.3, 0.3)

	ax.set_xlabel("FASTA Position")
	ax.set_yticks([])

	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	legend_elements = [
		Line2D([0],[0],
			   color='black',
			   lw=8,
			   label='18S rRNA'),
		Line2D([0],[0],
			   color=peak_color,
			   lw=6,
			   label='CLIP-seq Peaks'),
		Line2D([0],[0],
			   color=module_color,
			   lw=4,
			   label='Motif Segments')
	]

	ax.legend(
		handles=legend_elements,
		loc='upper right'
	)

	plt.tight_layout()

	plt.savefig(
		"CLIP_Module_Overview.png",
		dpi=600,
		bbox_inches='tight'
	)

	plt.show()












if __name__ == '__main__':
	main()