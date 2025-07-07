## RNAMotifModule

**A method to identify RNA structural motifs appearing in closed spatial proximity**

* Md Mahfuzur Rahaman<sup>†</sup>, mahfuz at ucf dot edu
* Shaojie Zhang*<sup>†</sup>, shzhang at cs dot ucf dot edu

<sup>†</sup>Department Computer Science, University of Central Florida, Orlando, FL 32816 USA \
*To whom correspondence should be addressed.

---

RNAMotifModule is a method to identify motif modules from RNA structure data. It is developed and tested in a **64-bit Linux** machine. For the basic features, only python (3.x recommended) is required to be installed on the computer with **64-bit Linux** environment. The preliminary files to run this code is included here.

### 1. Installation

#### 1.1 Install prelimineries

All recent Linux systems normally come with python installed by default. If not, please install `python` and `pip` before proceeding to the next step.

#### 1.2: Install required Python libraries

It is required to install several python libraries to run RNAMotifModule. These libraries are included in the [requirements.txt](requirements.txt) file. To install all required python libraries, please navigate to the RNAMotifModule home directory in the terminal and execute the following command.

```
pip install -r requirements.txt
```

### 2. Input Specifications

RNAMotifModule takes input from one or more files. Default ones are located in the [RNAMotifModule/input/](input) directory. Each line in that file represents a motif family. Each line starts with the motif family name, followed by a comma-separated list of motifs (the indices for motifs are PDB index, but FASTA index can also be used by setting a parameter in the configuration file). Please check the default input files [known_families_IL.csv](input/known_families_IL.csv) and [known_families_HL.csv](input/known_families_HL.csv) provided in the [input](input) directory to look into the formatting of input data in detail.

### 3. Commands for usage

```
usage: python3 RNAMotifModule.py [-h] [-i I] [-d D] [-o O] [-s] [-c] [-x X] [-r] [-t] [-n] [-p] [-e E] [-k] [-a A]
I - <input_file_name>     [Space separated list of input files.                                                                   Default: input/known_families_IL.csv input/known_families_HL.csv]
D - <distance_threshold>  [Distance threshold to be identified as spatially close.                                                Default: 5.0 Å]
O - <output_subdir>       [Subdirectory inside the "output" directory to save the results.                                        Default: '']
X - <randomized_datasets> [Randomized dataset to generate to calculate expected frequency of the motif modules.                   Default: 100]
E - <extension>           [Residues to extend beyond loop boundary to generate the partial PDB (loop.cif) files.                  Default: 0]
A - <atom_set>            [Atom set to use in calculating distance between residues (0: all, 1: backbone and sugar, 2: backbone). Default: 0]
```

**Examples:**

To identify motif modules from [sample_IL.csv](input/sample_IL.csv), the following command can be used:

```
python3 RNAMotifModule.py -i input/sample_IL.csv
```

More than one input files can also be provided in the following format:

```
python3 RNAMotifModule.py -i input/sample_IL.csv input/sample_HL.csv
```

To adjust the distance threshold in order to define spatially close, -d parameter can be used. For example, to identify motif modules from the previously used input files with a distance threshold 1.0, the following command can be used:

```
python3 RNAMotifModule.py -i input/sample_IL.csv input/sample_HL.csv -d 1.0
```

In order to save time while testing the process, we have pre-generated nearest residue data and spatial proximity data and saved them in corresponding pickle files. These pickle files can be utilized for the default input files using the -k parameter. Any runs without the -k param will generate the data from scratch and will replace the previous pickle files.

To identify motif modules from the input data from [known_families_IL.csv](input/known_families_IL.csv) and [known_families_HL.csv](input/known_families_HL.csv), by utilizing pre-generated pickle files use the following command:

```
python3 RNAMotifModule.py -k
```

To include motif module statistics for the same dataset using pre-generated pickle files use the following command:

```
python3 RNAMotifModule.py -s -k
```

As the default number of randomized dataset is 100, it takes 3-4 hours to complete the process even while utilizing pre-generated pickle files. Using -x param will be helpful to increase or decrease the number of randomized dataset while testing the code. For example, the following command can be used to generate motif module statistics with only 2 randomized dataset while using the pre-generated pickle data:

```
python3 RNAMotifModule.py -s -k -x 2
```

We provided pre-generated data for 357 internal loop motifs from 10 families and 415 hairpin loop motifs from 3 families. For any new dataset, it will automatically download and/or generate required data files (e.g. *.cif, *.fasta, *.aln, etc.) which might take some time. Please make sure to provide valid (not obsolete) PDB number in the input data.

RNAMotifModule can generate the nearest residue data for each of the motif modules. It can also generate partial PDB files for each motif group from different motif modules. Using -r parameter will continue to generate the nearest residue information files, -t will generate the statistics of these nearest residues for each motif module, and -p will generate the partial PDB files which can be directly loaded in PyMOL tool to visualize the motif modules. These parts are a bit time consuming as it needs to read through a lot of data files. Following command can be used to generate nearest residue data:

```
python3 RNAMotifModule.py -k -r
```

To generate nearest residue data as well as the nearest residue statistics, following command can be used:

```
python3 RNAMotifModule.py -k -r -t
```

To generate partial PDB files for each motif module, following command can be used:

```
python3 RNAMotifModule.py -k -p
```

Use of any of the -t or -p parameter will enable generating nearest residue data as both of them reads through PDB files to get corresponding information.

RNAMotifModule can also generate nearest residue data and the statistics for the non-module loops from the provided input motif families. To do so, -n parameter can be utilized. Following is the example command:

```
python3 RNAMotifModule.py -k -r -t -p -n
```

While generating the partial PDB files, -e parameter can be utilized to extend the loops to visualize them better in PyMOL. The default value is set to 0. Also, -o parameter can be utilized to create a subdirectory for the output files inside the default output directory.




### 4. Output specifications

By default, all the outputs will be generated inside the ‘output’ directory. A subdirectory will be created inside the ‘output’ directory if any subdirectory name is provided. Following are the possible output files and directories that can be found in the output directory:

Files:
```
identified_motif_modules.txt            Motif modules identified by RNAMotifModule
motif_module_stats.txt                  Statistics of the identified motif modules
nearest_protein_stat.txt                Nearest protein information of the motif modules
nonmodule_loops_info.txt                Non-module loops information
nearest_protein_stat_nonmodule.txt      Nearest protein information of the non-module loops
```


Directories:
```
nearest_residue_data                    Nearest residue data for the identified motif modules
partial_pdb_for_modules                 Partial PDB files for identified motif modules
nearest_residue_data_non_module_loops   Nearest residue data for non-module loops

```


### ACKNOWLEDGEMENTS

RNAMotifModule is developed for an NIH funded project (R01GM102515).
  
### CONTACTS

The author of the source code is Md Mahfuzur Rahaman. For bug reports or comments please contact mahfuz@ucf.edu.