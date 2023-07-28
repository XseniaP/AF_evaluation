# Evaluation of AlphaFold prediction accuracy

This script allows to evalute the prediction quality of the antibody heavy and light chains, to perform the calculations and combine all accuracy evaluation steps into a single pipeline. 

![Figure1_revised](https://github.com/XseniaP/AF_evaluation/assets/50076292/b61ee140-361e-4d7c-a37d-684aef8a440c)


The script is used to call the externally developed open-source products for alignment, superimposition, metrics calculation, and visualization, and then parses the results and moves the relevant information to the next step. In addition, the script provides on-demand calculations of the key metrics per domain and color coding, which are not covered by the externally developed packages. 

# Prerequisites

To run this script the following products and their prerequisites should be installed/compiled with the compilation files added to the main project folder:

(1) **BioPython** - it is used to align the sequences 
(2) **TM-align** - the TM-align.cpp can be downloaded from the official website (https://zhanggroup.org/TM-align/),following the compilation the executable file should be saved in the root directory of the project under the name "TMalign"
(3) **MM-align** - the MM-align.cpp can be downloaded from the official website (https://zhanggroup.org/MM-align/),following the compilation the executable file should be saved in the root directory of the project under the name "MMalign"
(4) **PyIR** should be installed according to the guidelines found of the official page https://github.com/crowelab/PyIR, "from crowelab_pyir import PyIR" is called in the script. Please note that PyIR requires a set of BLAST germline databases. Please follow the official instructions to set it up
(5) **json** - 
(4) 

# Running the script
