/*****************************************************************************
    ReNeGaTe

    Copyright (c) 2020-2022 Ali Hashemi
    		  2014-2022 Sana Bougueroua
                  
    Please cite:  J. Chem. Phys. 2018, 149 (18), 184102.         (DOI 10.1063/1.5045818 )
    		  J. Chem. Theory Comput. 2022, 18, 12, 7470–7482 (DOI 10.1021/acs.jctc.2c00404)
		  J. Chem. Inf. Model. 2023, XXXX, XXX, XXX-XXX (DOI 10.1021/acs.jcim.3c00660)

    ---------------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

/**
 * \file single_analysis.c 
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief The main program file to analyse a single trajectory.
 * \details This file contains the main steps used to analyse one trajectory.
 * The program read the coordinate file snapshot by snapshot, and for each snapshot compute the bonds based on geometric criteria.
 * If at least one bond has been changed, an isomorphism test is applied in order to identify the conformation already explored.
 * The coordinate file should have an \e .xyz or \e .mol as extension, and it is formatted as follows:
 *
 * <number of atoms>
 \n  comment line
 \n <atom_1> <X_1> <Y_1> <Z_1>
 \n <atom_2> <X_2> <Y_2> <Z_2>
 \n ... 
 \n<atom_n> <X_n> <Y_n> <Z_n>
 \n  <number of atoms>
 \n  comment line
 \n <atom_1> <X_1> <Y_1> <Z_1>
 \n <atom_2> <X_2> <Y_2> <Z_2>
 \n ... 
 \n<atom_m> <X_n> <Y_m> <Z_m>
 \n ... 
 * 
 * This file describes the molecule geometry by giving the number of atoms with Cartesian coordinates of atoms along the dynamics (for each snapshot). \n
 * Since an isomorphism test is applied, the number of atoms as well as the order of these atoms may change from one snapshot to another. 
 * 
 * In order to analyse a single trajectory after compilation through the makefile use the following line: 
 *
 * ./singanalysis -s {coordinate_file} [-i|-f|-c] [-d {bonds_to_analyze}] [-r] [-h 0/1] [-a]
 *
 * - -s : single analysis is performed
 * - file : trajectory input file for single analysis. It must have the extension : .xyz or .mol 
 * - [-i|-f|-c] : the type of the molecular system of the input trajectory:
 * - -i: MD of an isolated peptide in gas phase (by default: only hydrogen bonds dynamics are analysed)
 * - -f: Follow a collision leading to the fragmentation of a peptide (by default: only covalent bonds dynamics are analysed)
 * - -c: MD of an isolated cluster in gas phase (by default:  intermolecular electrostatic interactions and hydrogen bonds dynamics are analysed)
 * - [-d val] : specify the bonds types to take into account in the conformational analyses. It is 3 bits number that let you choose which kinds of bonds do you would take into account for the conformational dynamics analysis. These bits 1 to 3 are assigned to variables as follows: 
 * 
 * | Bit  | Variables                   | Value 
 * --------------------------------------------
 * | 1    | H-Bonds                     | 1 
 * --------------------------------------------
 * | 2    | Covalent bonds              | 2 
 * --------------------------------------------
 * | 3    | Electrostatic interactions  | 4 
 * --------------------------------------------
 * To only analyse H-bonds, choose 1. To analyse H-bonds with electrostatic interactions, choose 1+4=5. To analyse H-bonds, covalent bonds and electrostatic interactions, choose 1+2+4=7, etc.
 * 
 * - [-r] : rotational motion analysis is added
 * - [-h 0/1] : 1 to take into account angle in hydrogen bonds, 0 else (by default 1)
 * - [-a] : save the real energy average computed from the simulation for each conformation identified. This information is extracted from the comment line at each snapshot in the input file
 */

/*========================Libraries=========================*/ 
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "struct.h"
#include "init.h"
#include "conformation.h"
#include "io.h"
#include "memory.h"
#include "comfunct.h"
#include "fragment.h"
#include "usage.h"

/*========================Main program======================*/

/**
 * \fn int main(int argc, char *argv[])
 * \brief  Main steps to analyse a single molecular dynamics trajectory
 *
 * \param argc  number of parameters introduced by the user.
 * \param argv  table containing all parameters introduced by the user.
 * \return  \e 0 if the program has ended without problem.
 */


int main(int argc, char *argv[]){
	struct MyModel molecule; // Global variable used to contain all parameters and results of the analyses
  clock_t t= clock(); // Variable to calculate the running time 

//Read the input (arguments) and update parameters of variable \a molecule  (cut-off values, analysis type, bonds to analyse, etc.)
  struct MyOpt opt;
  opt = usage(argc,argv);
  get_anal_opt(&molecule,opt);
 //Initialization of the model, creation of the result directories and initialization of the adjacency matrices 
  init_model(&molecule);
  printf("read traj\n");
//Open the coordinate file (.xyz or .mol) 
  char inputF[512]; 
  sprintf(inputF,"%s/%s",molecule.inputDir,molecule.inputFN);
  FILE* traj;
  if((traj = fopen(inputF,"r")) == NULL){ printf("Can\'not open file for trajectory %s\n", inputF); exit(3);}
  struct confList* pred=NULL;
//Analyse the trajectory
  while(molecule.num_img< molecule.nbMolecule){
  //Get the conformational dynamics (isomorphism) for the current snapshot
    if(treat_snapshot(traj,&molecule)) {
      //At least one bond has been changed and so check the isomorphism 
      molecule.conformations = add_transit(&molecule,&pred,molecule.num_img); 
    } 
  //Go to next snapshot
    molecule.num_img++;
  //Tracking 
    printf("%d\n",molecule.num_img );
    if(molecule.num_img%1000==0){printf("\r %5d / %d ", molecule.num_img,molecule.nbMolecule); fflush(stdout);} 
  }  //End while 

//Close the coordinate file 
  fclose(traj);
//Update the last state (conformation) explored
  pred->imgList[2]=molecule.nbMolecule-1;
//Save the last period of the conformer 
  if((traj = fopen(molecule.periodFile,"a")) == NULL){ printf("Can\'not open file of periods\n"); exit(3);}
  fprintf(traj,"%d \t %d \t %d \n",pred->name,pred->imgList[1],pred->imgList[2]);
  pred->imgList[0]+=pred->imgList[2]-pred->imgList[1]+1;
  fclose(traj);
//Save the time evolution of hydrogen bonds if they have been taken into account in the trajectory analysis  
  if(bit_1(molecule.level,0) && molecule.nbHbond>0){
    //Save the HB dynamics
    save_Hbond_dynamics(&molecule);
  }
//Save the conformational dynamics of the analysed trajectory (the 2D graphs, .xyz files, graph of transitions, etc.) 
  printf("save conformations\n");
  save_conformations(&molecule);
// Analyse the fragments that have been observed along the trajectory
  if(opt.frg_on){
    printf("Fragments analysis\n");
    get_frag_dynamics(&molecule);
    // save_fragments(&molecule);
  }

//Extract the energy of conformations from the input file (should be present in the comment line at each snapshot : E= value)
  if(opt.eng_anal){
    save_real_energy(&molecule);
    save_diff_with_ref(&molecule);
  }

// Release memory allocated to variable \a molecule
  free_mymodel(&molecule);

//Display the running time
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
  printf("\n Running time = %f\n", time_taken);

//End the program
  return 0 ;
}       
