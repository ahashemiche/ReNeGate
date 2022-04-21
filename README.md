## ReNeGate 
Reaction Network Graph Theoretical tool for automated mechanistic studies in computational homogeneous catalysis


![Asset 1](https://user-images.githubusercontent.com/47638604/164565425-e93a3275-42fb-4534-81d9-169741a54b10.png)



RENEGATE (Reaction Network Graph Theoretical tool) is a free program for automated reaction network exploration and analysis. Conformer exploration, Reactive event identification and Reaction network analysis are the main steps taken for understanding the underlying mechanistic pathways in catalytic mixtures given the reaction mixture as the input. For conformer exploration configurational root mean squared deviations (RMSD) is used to bias metadynamics simulations to perform exhaustive configurational search on the potential energy surface of the given input catalytic mixture. The resulting metadynamics trajectories are then analysed with the developed basic graph theoretical algorithms in order to identify the conformations and transitions between them that occur along dynamics to form reaction networks. Reaction networks are then trimmed based on thermodynamic tresholds. Renegate is a series of codes written in C, bash and Python and can be run on Linux machines.

## Getting Started

---
Renegate program is free-ware and can be downloaded from [github]( https://github.com/ahashemiche/ReneGate). Hereafter the instructions of installations and requirements.

### Prerequisites

In order to run the Renegate program, you have to install the following software :

1. Crest(Conformer-Rotamer Ensemble Sampling Tool) : is used to conduct exhaustive RMSD biased explorations on the potential energy surface for the input (catalytic) mixture. [Download](https://github.com/grimme-lab/crest)

2. Gnuplot: is used to draw plots like the time evolution of hydrogen bonds. [Download](http://www.gnuplot.info/download.html)

3. Graphviz: is used to draw graphs as the 2D simplified graphs of conformations and the graph of transitions. [Download](https://www.graphviz.org/download/)

4. Nauty: is used to perform the isomorphism. [Download](http://pallini.di.uniroma1.it)

5. Python: is used for the trimming of chemical reaction networks after reactive event identification step [Download](https://www.python.org/)

### Installation and compilation

The code comes in a compressed file __Renegate.tar.gz__.
Use the following command to unpack the tar file :

```
tar -xvzf Renegate.tar.gz
```

After that, you will have the following arborescence:

* __Renegate/src__: containing the source and header files (_*.c_ , _*.h_ and *.sh),

* __Renegate/obj__: containing the object files (_*.o_),
 
* __Renegate/Makefile__: file for compilation,
 
* __Renegate/readMe__: file containing instructions of installations and requirements, and usage of the program. 
  
* __Renegate/examples__: containing some examples.

Then, in order to compile Renegate source code, use the make file provided :

``` 
 make -f Makefile  
```

Renegate contains three "Conformer Exploration", "Reactive event identification" and "Reaction network analysis" modules. 
If no error has occurred, this command-line "make" generates the following file:

* conformer_explorer: this script is used to do conformer exploration with the name of the input structure as the argument

``` 
 conformer_explorer input.xyz  
```

* singanalysis : this executable is used to analyse a single trajectory with a set of arguments (see usage section).

This file is located in the "bin" folder of the project.

You can run an example using the following command-line

```
  make example 
```

## Usage/ running code


---
Depending on the type of analysis, you can choose one module of Renegate by running one of the files generated by the "make" command-line. The usage of each module is described hereafter, with the list of required/optional arguments to be used.  

### Single analysis

In order to analyse a single trajectory file, the general syntax used is as follows :

```
 ./bin/singanalysis [-s] file [-i|-f|-c] [-d val] [-r] [-h car] [-a] [-l] [-t] [-g]
```

Where :

* __[-s]__ : means that you are asking for a single analysis, i.e., analysis of one trajectory.
 
* __file__ : the input file of _singanalysis module_. It contains the _xyz_ coordinates of the molecular system that you want to analyse. This input file should have a _.xyz_ or _.xmol_ extensions. User can add more extensions, this by editing the _"usage"_ function in the _"usage.c"_ file as described in the user-guide. The input file describes the molecular system geometry by giving the number of atoms with Cartesian coordinates of atoms along the dynamics. The formatting of this file format is as follows:

```
<number of atoms>
comment line
<atom_1> <X> <Y> <Z>
<atom_2> <X> <Y> <Z>
...
```

Note that for the current version the program takes into account 13 atoms types : Nitrogen (N), Hydrogen (H), Oxygen (O), Carbon (C) , Lithium (Li), Chlorine (Cl) , Sulfur (S), Manganese (Mn), Argon (Ar), Boron (B), Potassium (K), Brome (Br), and Silica (Si).

* __[-i|-f|-c]__ : this is a mandatory field where user specifies the type of system that he would analyse :
 
  * __-i__: MD trajectory for an isolated peptide in gas phase (by default: only the hydrogen bonds dynamics are analysed)
  * __-f__: CID (Collision Induced Dissociation) trajectory (by default: only the covalent bonds dynamics are analysed)
  * __-c__: MD trajectory for clusters in gas phase or catalysts (by default: the hydrogen bonds and the electrostatic interactions dynamics and the organometallic interactions are analysed). The cluster can contains either cation/anion as lithium, argon, potassium or chlorine, potassium or a metal as Manganese.

* __[-d val]__  : an optional field can be added to specify which bonds analyse. It is 4 bits number assigned (from left to right) to variables as follows:

|Bit|Variables |Value
|---|---|---
|1 |Hydrogen bonds| 1
|2 |Covalent bonds| 2
|3 |Electrostatic interactions | 4
|4 |Organometallic interactions | 8

Hence, to only analyse the hydrogen bonds dynamics, set __val__ to 1. To analyse the hydrogen bonds with the electrostatic interactions, set __val__ to 1+4=5. To analyse the hydrogen bonds, the covalent bonds and the electrostatic interactions, choose 1+2+4=7, and so on.

* __[-r]__ : an optional field to add a rotational motion analysis (identify the rotational axes). The rotational axes are identified based on the covalent bonds and the definition of isthmus on graphs (refer to the manual paper cited in the _References_ section for more details).

* __[-h car]__ : set the angle criterion for the hydrogen bonds analysis. __car__ is a binary value : set to 1 in order to add the angle criterion and 0 to remove it (a default value of 1 is used within the Renegate program). This is also an optional field.

* __[-a]__ : for some trajectories, energy value is saved at each step in the _.xyz_ input file (in the comment line at each step). Renegate can extract this value and compute for each conformation the energy average. User can get this information by setting the __-a__ option. This is also an optional field.

* __[-l]__ : an optional field that can be used in the analyses of catalysts. This field aims to perform a partial analysis of the molecular system. Only the fragments containing metals (_Mn_ in our case) are taking into account. One can make the same partial analysis with other type of atoms, and this by editing the function _"update\_atom\_in()"_ in the _"comfunct.c"_ file, as described in the user's manual.

* __[-t]__ : an optional field used to save the atoms involved in the dynamics of the analysed molecular system. With this option a second ".xyz" file is generated for each conformer identified along the trajectory. This ".xyz" will contain only the atoms involved in bonds that have been formed/broken at least once along the trajectory.

* __[-g]__ : an optional field used to analyze and identify the fragments that have been observed along the trajectory.

![Reaction Networks](https://user-images.githubusercontent.com/47638604/164558437-00eda2b3-d482-402f-acc0-09e1a19d7dd1.png)


## Code

---
Renegate consists of a series of sources files written in C, a header file, a make file for compilation, examples of input and output files stored in the Examples folder, a readMe file and a user's manual.

### Content

The code is organized within files according to the functions performed.
Files are organized into two types : (1) _*.c_ files which are the C source code files, and (2) _*.h_ files which are the header files including the function prototypes declaration, macro definitions, structures, constants, etc. Using such organization instead of a single file for the entire code let the latter more readable and easy to maintain. Moreover, the code is well commented in order to understand the different parts of the codes and the body of each function.
The Renegate code is composed of 4 main parts which are detailed in the user's manual:

* __Initialisation part__: there are two kinds of initialization, one for the input parameters  and one for the variables and structures used in the analyses.
* __Computational part__ : this part gathers all functions related to the conformational dynamics analyses.
* __Input/output part__ : it contains files used to save and display the different results.
* __Performance part__ : it contains the data structures, the function related to the memory management, check functions, etc.

### Output files/directories

Running the Renegate program will generate all results files in a directory having the same name of the input file/directory. Hereafter, the most important outputs provided by the Renegate program. Refer to the user's manual for the full description of the Renegate outputs.

* __Results files for a single analysis__ :
  * _"FileName\_conf"_ subdirectory : this subdirectory contains the 2D simplified graphs of the conformations explored along the trajectory.
  * _"FileName\_trans"_ subdirectory : this subdirectory contains the 2D simplified graphs of the conformations that are considered as intermediate states.
  * _"FileName\_xyz"_ subdirectory : for each conformation explored along the trajectory, Renegate saves one snapshot of this conformation into "coordinates\_confi.xyz" file. The snapshot saved is the first one where the conformation has been appeared.
  * _"FileName\_2\_xyz"_ subdirectory : same as _"FileName\_xyz"_, but the _.xyz_ files contain an extra column which shows the atoms involved in the bonds/interaction change.  
  * _"FileName\_frag"_ subdirectory : this subdirectory contains the 2D simplified graphs of the fragments explored along the trajectory.
  * _"FileName\_graphTrans.pdf"_ : this contains the graph of transitions labelled with the changes that have been occurred along the dynamics.
  * _"FileName\_perdFin.pdf"_ : this file contains the time evolution of the conformations explored along the trajectory.
  * _"FileName\_HbondsFin.pdf"_ : this file contains the time evolution of the hydrogen bonds identified along the trajectory.

Note that in our results files, atoms are labelled according to their chemical type and their order of appearance in the input file. For example if we have the list (_O_, _H_, _H_, _O_, _N_), the corresponding labels are respectively : (_O1_, _H1_, _H2_, _O2_, _N1_).

### Threshold values

A file __<constants.h>__ in the __src__ folder contains all threshold values used in the analyses done by the program, with the default values. These values include covalent radii, reference distances and angles for forming bonds/interactions, PBC values, etc. User can modify the values according to his needs.

## References

---
Refer to the following papers for more details

1. Bougueroua, S.; Barth, D.; Quessette, F.; Spezia, R.; Vial, S.; and Gaigeot, M.-P. Graph theory for automatic structural recognition in molecular dynamics simulations. J. Chem. Phys. (2018). DOI: 10.1063/1.5045818
2. manual\_paper / in progress
3. material\_paper / in progress

## Support

---
If you find a bug, please report it to the main developer (sana.bougueroua@univ-evry.fr).
Comments and suggestions are also welcome.

## Main developers

---

* _2014-2022_ Dr. Sana Bougueroua (LAMBE-UMR8587, University Evry Val d'Essonne, Université Paris Saclay, France)
* _2020-2022_ Dr. Ali Hashemi (TU Delft, Faculty of Applied Sciences, Netherlands)

## Acknowledgements

---
All co-authors acknowledge funding from the LABEX CHARMMMAT
(LABoratoire d'EXcellence CHARMMMAT, CHimie des ARchitectures Moléculaires Multifonctionnelles et des MATériaux).
Sana Bougueroua especially benefited from a PhD funding from LABEX CHARMMMAT. Pr. Marie-Pierre Gaigeot, Pr. Dominique Barth, Dr. Franck Quessette, Dr. Riccardo Spezia, Dr. Alvaro Cimas, Dr. Simone Pezzotti, and Prof. Evgeny Pidko are acknowledged for the very helpful and fruitful discussions.
Ali Hashemi acknowledges the financial support from the European Research Council (ERC) under the European Union's Horizon 2020 Research and Innovation Programme (grant agreement no. 725686). This work was sponsored by NWO Domain Science for the use of the national computer facilities. Ali Hashemi acknowledges Prof. Evgeny Pidko for very helpful and fruitful discussions.

## License

---
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
For more details see <http://www.gnu.org/licenses/>.

## Project status

---
The project is under development.
