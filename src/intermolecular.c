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
 * \file intermolecular.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief intermolecular electrostatic interactions analyses.
 * \details 
 * 				 						
 * We use the cut-off distances and angles to identify intermolecular interactions. 
 * According to the ion's atoms (cation, anion), threshold  distances are used. 
 * For this version we use the following distances :
 * - Li-Ar, Li-Atom, where the type of atom is not an Argon
 * - Ar-Ar , Ar-Atom
 * - NA-N , NA-O
 * - Cl-Cl , Cl-O
 * - Br-O
 * - I-O 
 * - K-O , K-atom
 * - F-H	
 * Note: user can add more atoms / distances. He has to define them in the following functions and in the constant.h file. 
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
#include "constants.h"
#include "struct.h"
#include "comfunct.h"
#include "io.h"
#include "check.h"
#include "intermolecular.h"

/*========================Functions=========================*/

/*========
Identify the intermolecular interactions and check if there has been a change with the previous snapshot
========*/
bool get_ion(struct MyModel* molecule){  
	int j=0;
	int i=0;
	bool change=false;

	//Browse all atoms and find if the ion approaches some of them.
	for(i=0;i<molecule->nbIon;i++){
		for(j=0;j<molecule->nbAtom;j++){
			if(verif_ion_bond(molecule,molecule->ionBond[i][molecule->nbAtom],j)){
				// printf("%d : %d : %d : %c \n", molecule->num_img, molecule->ionBond[i][molecule->nbAtom],j, molecule->img[molecule->ionBond[i][molecule->nbAtom]].atomType);
				if(molecule->ionBond[i][j]==-1){ //New intermolecular interaction
					molecule->ionBond[i][j]=1;
					molecule->ionBond[i][molecule->nbAtom+1]++; //CN increases (coordinate number)
					change = true;
					if(molecule->covBond[j][j] != -5 || j>molecule->ionBond[i][molecule->nbAtom])
						save_event(molecule,molecule->ionBond[i][molecule->nbAtom],j,'4','a');

				}
				//else : put else if we consider the atoms in position of a bond
			}
			else{ //Intermolecular interaction broken
				if(molecule->ionBond[i][j]!=-1 && j!=molecule->ionBond[i][molecule->nbAtom]){
					molecule->ionBond[i][j]=-1;
					molecule->ionBond[i][molecule->nbAtom+1]--; //CN decreases
					change = true;
					if(molecule->covBond[j][j] != -5 && j>molecule->ionBond[i][molecule->nbAtom])				
						save_event(molecule,molecule->ionBond[i][molecule->nbAtom],j,'5','a');
				}			
			}			
		}	
	}

	return change;
}

/*========
Verify if atoms of indexes index1 and index2 form an intermolecular interaction
========*/
bool verif_ion_bond(struct MyModel* molecule, int index1, int index2){ //check the parameter 	
	//FILE *traj_ang= fopen("traj_ang.txt","w");
	if(index1==index2) return false;
	if(molecule->covBond[index1][index2]>0 || molecule->covBond[index2][index1]>0) return false; // Atoms covalently bonded 
	double dist=0.00;
	if(molecule->sysType[0]=='w')
		dist = distance_pbc(molecule->img[index1], molecule->img[index2],molecule->tshVal.boxX,molecule->tshVal.boxY,molecule->tshVal.boxZ);
	else
		dist = distance(molecule->img[index1].x,molecule->img[index1].y,molecule->img[index1].z,molecule->img[index2].x,molecule->img[index2].y,molecule->img[index2].z);
	
	if(strcmp(molecule->img[index1].atomName,"Li")==0 || strcmp(molecule->img[index2].atomName,"Li")==0){ //One atom is Li
		if(strcmp(molecule->img[index1].atomName,"Ar")!=0 && strcmp(molecule->img[index2].atomName,"Ar")!=0){
			if(dist<=molecule->tshVal.distLiAt)
				return true;
			else
				return false;
		}
		else{//Second atom is Ar
			if(dist<=molecule->tshVal.distLiAr)
				return true;
			else
				return false;			
		}
	} //atoms != Li
	else{
		// NA-N , NA-O
		if(strcmp(molecule->img[index1].atomName,"Na")==0 ){
			if(strcmp(molecule->img[index2].atomName,"N")==0){
				if(dist<=molecule->tshVal.distNaN)
					return true;
				else
					return false;			
			}
			if(strcmp(molecule->img[index2].atomName,"O")==0){
				if(dist<=molecule->tshVal.distNaO)
					return true;
				else
					return false;			
			}
			return false;
		}
		// Cl-Cl , Cl-O
		if(strcmp(molecule->img[index1].atomName,"Cl")==0 ){
			if(strcmp(molecule->img[index2].atomName,"Cl")==0){
				if(dist<=molecule->tshVal.distClCl)
					return true;
				else
					return false;			
			}
			if(strcmp(molecule->img[index2].atomName,"O")==0){
				if(dist<=molecule->tshVal.distClO)
					return true;
				else
					return false;			
			}
			return false;
		}
		// Br-O
		if(strcmp(molecule->img[index1].atomName,"Br")==0 ){
			if(strcmp(molecule->img[index2].atomName,"O")==0){
				if(dist<=molecule->tshVal.distBrO)
					return true;
				else
					return false;	
			}		
			return false;
		}
		// I-O 
		if(strcmp(molecule->img[index1].atomName,"I")==0 ){
			if(strcmp(molecule->img[index2].atomName,"O")==0){
				if(dist<=molecule->tshVal.distIO)
					return true;
				else
					return false;			
			}
			return false;
		}
		// K-O , K-atom
		if(strcmp(molecule->img[index1].atomName,"K")==0 || strcmp(molecule->img[index2].atomName,"K")==0){
			if(strcmp(molecule->img[index2].atomName,"C")==0) return false; //ignore the carbone for the moment
			if(strcmp(molecule->img[index2].atomName,"O")==0){
				if(dist<=molecule->tshVal.distKO){
					return true;
				}
				else{
					if(dist<=molecule->tshVal.distKAt && molecule->sysType[0]!='w')
						return true;
					else
						return false;			
				}
			}
		}
		else{ 
			// F-H	
			if(strcmp(molecule->img[index1].atomName,"F")==0){
				if(strcmp(molecule->img[index2].atomName,"H")==0){
					if(dist<=molecule->tshVal.distFH)
						return true;
					else
						return false;								
				}
				return false;
			} 
			else{ //Ar-atom 
				if(strcmp(molecule->img[index1].atomName,"Ar")==0 || strcmp(molecule->img[index2].atomName,"Ar")==0){
					if(dist<=molecule->tshVal.distArAt){
						// if(strcmp(molecule->img[index1].atomName,"Ar")==0 && strcmp(molecule->img[index2].atomName,"Ar")!=0){
						// 	double ang = (angle(molecule->img[get_index(molecule,"Li", 1)].x,molecule->img[get_index(molecule,"Li", 1)].y,molecule->img[get_index(molecule,"Li", 1)].z,molecule->img[index1].x,molecule->img[index1].y,molecule->img[index1].z,molecule->img[index2].x,molecule->img[index2].y,molecule->img[index2].z))*180/PI;
						// 	// fprintf(traj_ang,"file=%s, numImg=%d, ang(%d,%d,%d)=%lf \n", molecule->inputFN, molecule->num_img, index1,get_index(molecule,"Li", 1),index2 , ang);
						// 	// fclose(traj_ang);
						// 	if(ang>=79.5 && ang<=99.5)
						// 		return true;
						// 	else 
						// 		return false;
						// }
						// else{
						// 	if(strcmp(molecule->img[index1].atomName,"Ar")!=0 && strcmp(molecule->img[index2].atomName,"Ar")==0){
						// 		double ang = (angle(molecule->img[get_index(molecule,"Li", 1)].x,molecule->img[get_index(molecule,"Li", 1)].y,molecule->img[get_index(molecule,"Li", 1)].z,molecule->img[index2].x,molecule->img[index2].y,molecule->img[index2].z,molecule->img[index1].x,molecule->img[index1].y,molecule->img[index1].z))*180/PI;
						// 		// fprintf(traj_ang,"file=%s, numImg=%d, ang2(%d,%d,%d)=%lf \n", molecule->inputFN, molecule->num_img, index1,get_index(molecule,"Li", 1),index2 , ang);
						// 		// fclose(traj_ang);
						// 		if(ang>=84.5 && ang<=94.5)
						// 			return true;
						// 		else 
						// 			return false;					
						// 	}
						// 	else
						// 		return true;
						// }
						return true;
					}
					else
						return false;					

				}
			}

		}	
	}
	return false;
}
