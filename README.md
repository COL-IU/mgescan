MGEScan
===============================================================================

A probabilistic approach for *de novo* identification and classification of mobile genetic elements in eukaryotic genomes. 

Prerequisites
-------------------------------------------------------------------------------

* git v2.0.0
* python v2.7.13
* virtualenv v15.1.0
* perl v5.10.1
* HMMER v3.1b1
* EMBOSS v6.6.0
* Tandem Repeats Finder v4.07b
* NCBI BLAST+ v2.2.28

Optional Prerequisites
-------------------------------------------------------------------------------

* Galaxy v17.01 (If you intend to use Galaxy platform to run MGEScan)
* RepeatMasker v4.0.5 (If you indend to use RepeatMasker in the preprocessing step of MGEScan-LTR)


Installation
-------------------------------------------------------------------------------

#### Step 1: Exporting the necessary environment variables
```sh
export MGESCAN_HOME=$HOME/mgescan4
export MGESCAN_VENV=$HOME/venv/mgescan
```
Note: Instead of the paths specified above, you can provide your own paths where you would like to install MGEScan and the Python virtual environment for MGEScan. 

#### Step 2 (optional - if you wish to use MGEScan tool for Galaxy): Export GALAXY_HOME environment variable so that MGEScan knows where Galaxy is installed
```sh
export GALAXY_HOME=$HOME/galaxy
```
Note: Replace the above path with the path where you've installed Galaxy. 

#### Step 3: Create the Python virtual environment for MGEScan
```sh
virtualenv $MGESCAN_VENV
source $MGESCAN_VENV/bin/activate
```

#### Step 4: Install MGEScan
```sh
git clone https://github.com/COL-IU/mgescan.git
cd mgescan
python setup.py install
```

Command Line Tool (mgescan)
-------------------------------------------------------------------------------
**Note: Check the following before running MGEScan**
1. Make sure the environment variables (mentioned above) are exported every time before you run MGEScan. 
2. Make sure the Python virtual environment for MGEScan is activated every time before you run MGEScan. 

```sh
Usage:
    mgescan all <genome_dir> [--output=<data_dir>] [--mpi=<num>]
    mgescan ltr <genome_dir> [--output=<data_dir>] [--mpi=<num>]
    mgescan nonltr <genome_dir> [--output=<data_dir>] [--mpi=<num>]
    mgescan dna <genome_dir> [--output=<data_dir>] [--mpi=<num>]
    mgescan (-h | --help)
    mgescan --version
```

Citation
-------------------------------------------------------------------------------

How to cite MGEScan on Galaxy [here]

References
-------------------------------------------------------------------------------

1. H. Lee, M. Lee, W. M. Ismail, M. Rho, G. C. Fox, S. Oh, H. Tang (2016), MGEScan: 
   a Galaxy-based system for identifying retrotransposons in genomes, Bioinformatics, 
   32(16): 2502-2504. 

2. M. Rho, S. Schaack, X. Gao, S. Kim, M. Lynch and H. Tang (2010), LTR
   retroelements in the genome of Daphnia pulex, BMC Genomics, 11:425. Pubmed. 

3. M. Rho and H. Tang (2009), MGEScan-nonLTR: computational identification and
   classification of Non-LTR retrotransposons in eukaryotic genomes. Nucleic Acid
   Res, 37(21):e143. Free fulltext at NAR online 

4. M. Rho, J. H. Choi, S. Kim, M. Lynch and H. Tang (2007), De novo
   identification of LTR retrotransposons in eukaryotic genomes. BMC Genomics,
   8:90. Pubmed. 

Web Sites
-------------------------------------------------------------------------------

* [MGEScan-LTR](http://darwin.informatics.indiana.edu/cgi-bin/evolution/daphnia_ltr.pl)
* [MGEScan-nonLTR](http://darwin.informatics.indiana.edu/cgi-bin/evolution/nonltr/nonltr.pl)

License
-------------------------------------------------------------------------------

Copyright (C) 2015. See the LICENSE file for license rights and limitations
(GPL v3).

This program is part of MGEScan.

MGEScan is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
