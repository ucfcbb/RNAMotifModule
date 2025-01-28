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

RNAMotifModule takes input from one or more files. Default ones are located in the [RNAMotifModule/input/](input) directory. Each line in that file represents a motif family. Each line starts with the motif family name, followed by a comma-separated list of motifs (the indices for motifs are PDB index, but FASTA index can also be used by setting a parameter in the configuration file). Please check the default input files ([known_families_IL.csv](input/known_families_IL.csv)) and ([known_families_HL.csv](input/known_families_HL.csv)) provided in the [input](input) directory to look into the formatting of input data in detail.

### 3. Commands for usage

```
usage: python3 RNAMotifModule.py [-h] [-i I] [-d D] [-o O] [-s] [-x X] [-r] [-t] [-n] [-p] [-e E] [-k]
I - <input_file_name> [Space separated list of input files.]
D - <distance_threshold> [Distance threshold to be identified as spatially close.]
O - <output_subdir> [Subdirectory inside the "output" directory to save the results.]
X - <randomized_datasets> [Randomized dataset to generate to calculate expected frequency of the motif modules.]
E - <extension> [Residues to extend beyond loop boundary to generate the partial PDB (loop.cif) files.]
```

**Examples:**

To identify motif module from the input data from ([known_families_IL.csv](input/known_families_IL.csv)) and ([known_families_HL.csv](input/known_families_HL.csv)), by utilizing pre-generated pickle files use the following command:

```
python3 RNAMotifModule.py -k
```

To include motif module statistics for the same dataset using pre-generated pickle files use the following command:

```
python3 RNAMotifModule.py -s -k
```

We provided pre-generated data for 358 internal loop motifs from 10 families and 415 hairpin loop motifs from 3 families. For any new dataset, it will automatically download and/or generate required data files (e.g. *.cif, *.fasta, *.aln, etc.) which might take some time. Please make sure to provide valid (not obsolete) PDB number in the input data.

### 4. Output specification

By default, all the outputs will be generated inside the ‘output’ directory. A subdirectory will be created inside the ‘output’ directory if any subdirectory name is provided.

### ACKNOWLEDGEMENTS

RNAMotifModule is developed for an NIH funded project (R01GM102515).
  
### CONTACTS

The author of the source code is Md Mahfuzur Rahaman. For bug reports or comments please contact mahfuz@ucf.edu.