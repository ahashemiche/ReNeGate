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
 * \file organometallic.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief organometallic interactions analyses.
 * \details 
 * 				 						
 * We use the cut-off distances to identify organometallic interactions. 
 * According to the metal atom, threshold  distances are used. 
 * For this version we use the following distances :
 * - Metal-Atom
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
#include "organometallic.h"

/*========================Functions=========================*/

/*========
Identify the organometallic interactions and check if there has been a change with the previous snapshot
========*/
bool get_metal_dynamics(struct MyModel* molecule){  
	int j=0;
	int i=0;
	bool change=false;
	//Browse all atoms and find if the metal approaches some of them.
	for(i=0;i<molecule->nbMetal;i++){
		for(j=0;j<molecule->nbAtom;j++){
			if(verif_metal_bond(molecule,molecule->metalBond[i][molecule->nbAtom],j)){
				if(molecule->metalBond[i][j]==-1){ //New organometallic interaction
					molecule->metalBond[i][j]=1;
					molecule->metalBond[i][molecule->nbAtom+1]++; //CN increases (coordinate number)
					change = true;
					// if( j>molecule->metalBond[i][molecule->nbAtom])//molecule->covBond[j][j] != -5 ||
						save_event(molecule,molecule->metalBond[i][molecule->nbAtom],j,'6','a');

				}
				//else : put else if we consider the atoms in position of a bond
			}
			else{ //organometallic interaction broken
				if(molecule->num_img !=0){
					if(molecule->metalBond[i][j]==1){ // && j!=molecule->metalBond[i][molecule->nbAtom]
						molecule->metalBond[i][j]=-1;
						molecule->metalBond[i][molecule->nbAtom+1]--; //CN decreases
						change = true;
						// if(j>molecule->metalBond[i][molecule->nbAtom] )	//molecule->covBond[j][j] != -5 && 			
							save_event(molecule,molecule->metalBond[i][molecule->nbAtom],j,'7','a');
					}			
				}
			}			
		}	
	}
	return change;
}

/*========
Verify if atoms of indexes index1 and index2 form an organometallic interaction
========*/
bool verif_metal_bond(struct MyModel* molecule, int index1, int index2){ //check the parameter 	
	double dist = distance(molecule->img[index1].x,molecule->img[index1].y,molecule->img[index1].z,molecule->img[index2].x,molecule->img[index2].y,molecule->img[index2].z);
	// if(dist<2.50 && (strcmp(molecule->img[index2].atomName,"O")==0 || strcmp(molecule->img[index2].atomName,"H")==0))printf("distance %s is %lf \n", molecule->img[index2].atomName,  dist);
	if(index1==index2 ) return false; //The metal atom
	if(molecule->covBond[index1][index2]>0 || molecule->covBond[index2][index1]>0) return false; //atoms covalently bonded

	if(strcmp(molecule->img[index2].atomName,"C")==0) return false; //Carbon not allowed
	if(molecule->covBond[index1][index2]>0 || molecule->covBond[index2][index1]>0 ) return false; // The metal is covalently bonded to atom of index  index2
	
	if(strcmp(molecule->img[index1].atomName,"Mn")==0 && strcmp(molecule->img[index2].atomName,"O")==0){ //check 
		if(dist<=molecule->tshVal.distMnO)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Mn")==0 && strcmp(molecule->img[index2].atomName,"H")==0){
		if(dist<=molecule->tshVal.distMnH)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Mn")==0 && strcmp(molecule->img[index2].atomName,"Br")==0){
		if(dist<=molecule->tshVal.distMnBr)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Mn")==0 && strcmp(molecule->img[index2].atomName,"F")==0){
		if(dist<=molecule->tshVal.distMnF)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Mn")==0 && strcmp(molecule->img[index2].atomName,"S")==0){
		if(dist<=molecule->tshVal.distMnS)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Mn")==0 && strcmp(molecule->img[index2].atomName,"P")==0){
		if(dist<=molecule->tshVal.distMnP)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Ru")==0 && strcmp(molecule->img[index2].atomName,"F")==0){
		if(dist<=molecule->tshVal.distRuF)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Ru")==0 && strcmp(molecule->img[index2].atomName,"Cl")==0){
		if(dist<=molecule->tshVal.distRuCl)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Ru")==0 && strcmp(molecule->img[index2].atomName,"P")==0){
		if(dist<=molecule->tshVal.distRuP)
			return true;
		else
			return false;		
	}
	if(strcmp(molecule->img[index1].atomName,"Au")==0 && strcmp(molecule->img[index2].atomName,"Au")==0){
		printf("dist=%.2lf\n",dist);
		if(dist<=molecule->tshVal.distAuAu)
			return true;
		else
			return false;		
	}

	return false;
}
