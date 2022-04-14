/*****************************************************************************
    GATEwAY - GrAph ThEory for conformAtional dYnamics 

    Copyright (c) 2014-2022 Sana Bougueroua
                  2020-2022 Ali Hashemi
    Please cite:  J. Chem. Phys. 2018, 149 (18), 184102.         (DOI 10.1063/1.5045818 )
		   	
    This file written by Sana Bougueroua.

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
 * \file covalent.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief  covalent bonds treatment.
 * \details 
 *
 * This file contains all functions used to compute covalent bonds. We use the "covalent radius" (covR) of atoms to identify
 * the covalent bonds. This calculation is done on two steps:
 *
 * Step (I): put a covalent bond between two atoms (i,j) 	
 * each at a distance less < (covRi + covRj)*1,3 (margin 30%).
 *
 * Step (II): handle missing bonds and overflow.				
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
#include "memory.h"
#include "io.h"
#include "check.h"
#include "covalent.h"

/*========================Functions=========================*/

/*========
Compute covalent bonds for the molecule at instant molecule->num_img
========*/
bool get_cov(struct MyModel* molecule){
	int i=0;
	int j=0 ;
	int **copy_cov;
//Step (0) : save the previous covalent bonds
	if(!molecule->nbAtomChange){
		copy_cov= allocate_matrix(molecule->nbAtom,molecule->nbAtom,0,"copy_cov"); // copy of covalent bonds
		copy_matrix(molecule->covBond,copy_cov,molecule->nbAtom);
		for (i = 0; i < molecule->nbAtom; i++) for (j = 0; j < molecule->nbAtom; j++) molecule->covBond[i][j]=0;
	}

//Step (I) : create a simple covalent bond between atoms at a distance less < (covRi + covRj)*1,4
	//Complexity : O(nbAtom x nbAtom/2)
	for (i = 0; i < molecule->nbAtom; i++){
		if(strcmp(molecule->img[i].atomName,"Ar")==0 || strcmp(molecule->img[i].atomName,"Na")==0 || strcmp(molecule->img[i].atomName,"Li")==0 || strcmp(molecule->img[i].atomName,"I")==0 || strcmp(molecule->img[i].atomName,"K")==0 ||  strcmp(molecule->img[i].atomName,"Br")==0  ){
			molecule->covBond[i][i] = -5; //ignore cation / anion
		} 
		else{
			molecule->covBond[i][i] = 0;
		}
		for (j = i+1; j < molecule->nbAtom; j++){
			if (strcmp(molecule->img[i].atomName,"H")!=0 && strcmp(molecule->img[i].atomName,"H")!=0){
				//We ignore Hydrogen atoms at the first step
				if( molecule->covBond[j][i] !=-1 && molecule->covBond[j][i] !=1){ 
			   		//if the atom is not implicated in a proton-transfer & H-bond we can make change
					//add a covalent bond
					if (verif_cov_bond(molecule,i,j)){  
						molecule->covBond[i][j]=1;
						molecule->img[j].bondMax--;
						molecule->img[i].bondMax--;
					}
					else
						molecule->covBond[i][j]=0;
				}
				else 
					if(molecule->img[i].atomType=='H') {//Hydrogen
						if((molecule->covBond[j][i] ==1 && sens_Hbond(molecule,i)==0) || (molecule->covBond[j][i] ==-1 && sens_Hbond(molecule,i)==1)){
							molecule->img[j].bondMax--; 
							molecule->img[i].bondMax--; 									
						}
					}
					else {
						if((molecule->covBond[j][i] ==1 && sens_Hbond(molecule,j)==0) || (molecule->covBond[j][i] ==-1 && sens_Hbond(molecule,j)==1)){
							molecule->img[j].bondMax--; 
							molecule->img[i].bondMax--; 									
						}
					}
			}
		}
	}
//Step (II) : handle missing bonds and overflow
	//Complexity : O(nbAtom x nbAtom/2)
  	check_cov_bond(molecule);
//Display covalent bonds
	#ifdef DISP_COV
		display_cov_bonds(molecule);
	#endif  	

//Check if there is a dynamics on covalent bonds
	if(!molecule->nbAtomChange){
  		if(verif_change_cov(molecule,molecule->covBond,copy_cov,molecule->nbAtom) || molecule->num_img==0){
	  		free_matrix(copy_cov, molecule->nbAtom);
	  		return true;
	  	}
	  	else{
		   	free_matrix(copy_cov, molecule->nbAtom);
	  		return false; 		
	  	}   	
	}	
	else{
		return true;
	}
}

/*========
Check the covalent bonds found in the get_cov function at instant "molecule->num_img".
========*/
void check_cov_bond(struct MyModel* molecule){
	//Delete all triangle structures 
	check_triangle_bond(molecule);

	int i=0; 
	bool change=true;
	while(change){ 
		while(i<molecule->nbAtom && molecule->img[i].bondMax==0)  i++;
		if(i==molecule->nbAtom) change =false;
		else{
			//we let metal and the phosphorus at the end 
			if(strcmp(molecule->img[i].atomName,"Mn")==0 || strcmp(molecule->img[i].atomName,"Ru")==0 || strcmp(molecule->img[i].atomName,"P")==0){
				i++;
			}
			else{
				if(molecule->img[i].bondMax>0){ //bonds missed || strcmp(molecule->img[i].atomName,"Mn")==0
					int j=check_nearb(molecule,i,'-');//check if there is an adjacent atom with a missed bond
					if (j!=-1){
							if (verif_cov_bond(molecule,i,j)){   
								add_bond(molecule,i,j);//add a bond between i and j
								i=0;								
							}
							else
								i++;
					} 
					else{
						if(molecule->img[i].atomType=='S' && molecule->img[i].bondMax==1){
							molecule->img[i].bondMax=0;
						}
						if(strcmp(molecule->img[i].atomName,"Mn")==0){
							molecule->img[i].bondMax--;
						}
						if(strcmp(molecule->img[i].atomName,"Ru")==0){
							molecule->img[i].bondMax--;
						}
						if(strcmp(molecule->img[i].atomName,"P")==0){
							molecule->img[i].bondMax--;
						}
						i++;					
					}
				}
				else{ //Overflow 
					if((strcmp(molecule->img[i].atomName,"N")==0 ||strcmp(molecule->img[i].atomName,"B")==0 ) && molecule->img[i].bondMax==-1 ){
						molecule->img[i].bondMax++;
					}
					else{
						int j=check_nearb(molecule,i,'+'); //check if there is an adjacent atom with an overflow				
					    // printf("i=%d ; j=%d ; bondMax=%d ; type=%s\n",i,j, molecule->img[i].bondMax, molecule->img[i].atomName);

						if (j!=-1){
							del_bond(molecule,i,j); //delete the bond between i and j 
		 					i=0;
						}
						else
							i++;											
					}
				}
			}
		}
	}  
	//Get the number of covalent bonds
	int nbCov =0;
	int j=0;
	for (i = 0; i < molecule->nbAtom; i++) for (j = i+1; j < molecule->nbAtom; j++)if(molecule->covBond[i][j] >0) nbCov++;
	molecule->nbCov = Max(nbCov,molecule->nbCov);
}
/*========
Check if there are triangles and delete the longest bond for each triangle
========*/
void check_triangle_bond(struct MyModel* molecule){
	int i;
	int j;
	int k;
	for(i=0;i<molecule->nbAtom;i++){
		for(j=i+1;j<molecule->nbAtom;j++){
			if(molecule->covBond[i][j]!=0){
				for(k=i+1;k<molecule->nbAtom;k++){
					if(k!=i && k!=j && molecule->covBond[i][k]!=0 && molecule->covBond[j][k]!=0){
						//Triangle found 
						double dist1 = distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z);
						double dist2 = distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[k].x,molecule->img[k].y,molecule->img[k].z);
						double dist3 = distance(molecule->img[k].x,molecule->img[k].y,molecule->img[k].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z);
						double minDist = Max(dist1,dist2);
						minDist = Max(minDist,dist3);
						//Delete the longest bond
						if(minDist==dist1) del_bond(molecule,i,j);
						if(minDist==dist2) del_bond(molecule,i,k);
						if(minDist==dist3) del_bond(molecule,j,k);
					}
				}				
			}			
		}
	}
}
/*========
Provide the internal covalent bonds NI-NI of molecule.
========*/
void get_inter_bonds(struct MyModel* molecule){
	int i=0;
	int j=0;
	for (i = 0; i < molecule->nbAtom; i++){
		int cpt=0;
		for (j = i+1; j <molecule->nbAtom; j++){
			if(check_inter_node(molecule,i,j) && check_inter_node(molecule,j,i)){
				molecule->interBond[i][j]=1;
				cpt++;
			}
			else
				molecule->interBond[i][j]=0;
			molecule->interBond[j][i]=0;			
		}
		molecule->interBond[i][i]=cpt; //# of internal bonds of the atom i
	}
}

/*========
Check if a-b is an internal covalent bond
========*/
bool check_inter_node(struct MyModel* molecule, int a, int b){
	
	if(molecule->img[a].atomType=='H'|| molecule->img[b].atomType=='H')
		return false; //hydrogen atoms are leaves
	if((a<b  && molecule->covBond[a][b]==0) || (a>b && molecule->covBond[b][a]==0))  
		return false;
	int j=0;
	//research by lines
	while (j<a){
		if(molecule->covBond[j][a]!=0 && j!=b && molecule->img[j].atomType!='H')
			return true;
		else
			j++;
	}
	j++;
	//research by columns
	while (j<molecule->nbAtom){
		if(molecule->covBond[a][j]!=0 && j!=b && molecule->img[j].atomType!='H')
			return true;
		else
			j++;
	}
	return false;
}

/*========
Verify changes in the covalent bonds
========*/
bool verif_change_cov(struct MyModel *molecule,int **covBond,int **covCopy,int size){
	int i =0 ;
	int j=0 ;
	bool covChang=false;
	for (i = 0; i < size; i++)
		for (j = i+1; j < size; j++)
			if(covCopy[i][j]!=covBond[i][j] && covBond[j][i]!=-1 && molecule->num_img !=0){
				//Save the change 
				if(covBond[i][j]==0 || covCopy[i][j]==0){ //we don't take into account change from double/triple bond to single/double
					if(molecule->img[i].atomType=='H') molecule->covBond[i][i]=11; //make change
					if(molecule->img[j].atomType=='H') molecule->covBond[j][j]=11;				
					covChang=true;
					if(covBond[i][j]> covCopy[i][j])
						save_event(molecule,i,j,'2','a');
					else
						save_event(molecule,i,j,'3','a');				

				}
			}
	return covChang;
}

/*========
Check if there is a bond between atoms i and j
========*/
bool verif_cov_bond(struct MyModel* molecule, int i, int j) {
	// if(molecule->img[i].atomType=='H' && molecule->img[j].atomType=='H') 
	// 	return false ;
	
	if(strcmp(molecule->img[i].atomName,"I")==0 || strcmp(molecule->img[i].atomName,"Au")==0 || strcmp(molecule->img[i].atomName,"Na")==0 || strcmp(molecule->img[i].atomName,"Li")==0 || strcmp(molecule->img[i].atomName,"Ar")==0 || strcmp(molecule->img[i].atomName,"K")==0 || strcmp(molecule->img[i].atomName,"Br")==0 )
		return false;
	if(strcmp(molecule->img[j].atomName,"I")==0 || strcmp(molecule->img[j].atomName,"Au")==0 || strcmp(molecule->img[j].atomName,"Na")==0 || strcmp(molecule->img[j].atomName,"Li")==0 || strcmp(molecule->img[j].atomName,"Ar")==0 || strcmp(molecule->img[j].atomName,"K")==0 || strcmp(molecule->img[j].atomName,"Br")==0 )
		return false;
	if((strcmp(molecule->img[i].atomName,"Cl")==0 || strcmp(molecule->img[j].atomName,"Cl")==0) && molecule->sysType[0]!='w')
		return false;
	if((strcmp(molecule->img[i].atomName,"Mn")==0 && strcmp(molecule->img[j].atomName,"O")==0) || (strcmp(molecule->img[i].atomName,"O")==0 && strcmp(molecule->img[j].atomName,"Mn")==0) )
		return false;
	
	//Compute distances	
	double dist;
	//We can use either the simple euclidean distance or distance with respect to the Periodic Boundary Conditions (PBC)
	if(molecule->sysType[0]=='w')
		dist = distance_pbc(molecule->img[i], molecule->img[j],molecule->tshVal.boxX,molecule->tshVal.boxY,molecule->tshVal.boxZ);
	else
		dist = distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z);
	
	double covr1 = get_covr_value(molecule,molecule->img[i].atomName);
	double covr2 = get_covr_value(molecule,molecule->img[j].atomName);
	//Condition of covalent bond using the covalent radius and margin of 40%
	if(i!=j && dist<= (covr1+covr2)*molecule->tshVal.covMarg) return true ; else return false;
	
}

/*========
Get the covalent radius of the atom with type atomName
========*/
double get_covr_value(struct MyModel *molecule, char* atomName){
	
	//replace switch by if(strcmp(molecule->img[j].atomName,"Mn")
	if(strcmp(atomName,"C")==0) return molecule->tshVal.covrC ; 
	if(strcmp(atomName,"H")==0) return molecule->tshVal.covrH ; 
	if(strcmp(atomName,"O")==0) return molecule->tshVal.covrO ; 
	if(strcmp(atomName,"N")==0) return molecule->tshVal.covrN ; 
	if(strcmp(atomName,"S")==0) return molecule->tshVal.covrS ; 
	if(strcmp(atomName,"Mn")==0) return molecule->tshVal.covrMn ; 
	if(strcmp(atomName,"B")==0) return molecule->tshVal.covrB ; 
	if(strcmp(atomName,"F")==0) return molecule->tshVal.covrF ; 
	if(strcmp(atomName,"P")==0) return molecule->tshVal.covrP ; 
	if(strcmp(atomName,"Ru")==0) return molecule->tshVal.covrRu ; 
	if(strcmp(atomName,"Cl")==0) return molecule->tshVal.covrCl ;
	if(strcmp(atomName,"Si")==0) return molecule->tshVal.covrSi ;
	if(strcmp(atomName,"Zn")==0) return molecule->tshVal.covrZn ; 

	return 0.00;
}

/*========
Verify the sens of H-bonds (proton transfer or not)
========*/
int sens_Hbond(struct MyModel* molecule, int index){
	int i=0 ; 
	while(i<molecule->nbHAtom && molecule->bondHList[i].numH != index) i++;
	if(i==molecule->nbHAtom) return -1 ; else return molecule->bondHList[i].sens-'0';
}

/*========
Add a covalent bond between two atoms with indexes i and j
========*/
void add_bond(struct MyModel *molecule, int i, int j){
	if(i<j)
		molecule->covBond[i][j]++;
	else
		molecule->covBond[j][i]++;
	molecule->img[i].bondMax--;
	molecule->img[j].bondMax--;

}

/*========
Delete a covalent bond between two atoms with indexes i and j 
========*/
void del_bond(struct MyModel *molecule, int i, int j){
	if(i<j)
		molecule->covBond[i][j]--;
	else
		molecule->covBond[j][i]--;
	molecule->img[i].bondMax++;
	molecule->img[j].bondMax++;
}

/*========
Check if i has a covalent bond with an atom under some conditions (take the nearest for the '-' and the farest for the '+')
========*/
int check_nearb(struct MyModel *molecule, int i,int type){
	int j=0;
	switch(type){
		case '-':{ //bonds missed
			int k=-1;
			//research by lines
			while (j<i){
				if(molecule->covBond[j][i]!=0 || molecule->img[i].atomType=='H')
					if(molecule->img[j].bondMax>0){
						if(k==-1)
							k=j;
						else
							if(distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z)<distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[k].x,molecule->img[k].y,molecule->img[k].z))
								k=j;
						j++;
					}
					else
						j++;
				else
					j++;
			}
			j++;
			//research by clowns
			while (j<molecule->nbAtom){
				if(molecule->covBond[i][j]!=0 || molecule->img[i].atomType=='H')
					if(molecule->img[j].bondMax>0){
						if(k==-1)
							k=j;
						else
							if(distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z)<distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[k].x,molecule->img[k].y,molecule->img[k].z))
								k=j;
						j++;
					}
					else
						j++;
				else
					j++;
			}
			return k;
		}
		break;
		case '+':{ //overflow
			int k=-1;
			// research by lines
			while (j<i){
				if(molecule->covBond[j][i]!=0)
					if(molecule->img[j].bondMax<=0 || strcmp(molecule->img[j].atomName,"Mn")==0 || strcmp(molecule->img[j].atomName,"Ru")==0){
						if(k==-1)
							k=j;
						else
							if(distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z)>distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[k].x,molecule->img[k].y,molecule->img[k].z))
								k=j;
						j++;
					}
					else
						j++;
				else
					j++;
			}
			//research by colowns
			j++;
			while (j<molecule->nbAtom){
				if(molecule->covBond[i][j]!=0)
					if(molecule->img[j].bondMax<=0 || strcmp(molecule->img[j].atomName,"Mn")==0){
						if(k==-1)
							k=j;
						else
							if(distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z)>distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[k].x,molecule->img[k].y,molecule->img[k].z))
								k=j;
						j++;
					}
					else
						j++;
				else
					j++;
			}
			return k;
		}
		break;
		default : printf("ERROR on the index type\n");
	}
	return -1;
}

