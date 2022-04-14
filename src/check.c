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
 * \file check.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief Functions to display the content of variables.
 * \details 
 *
 * This file contains all functions used to display the content of variables.
 */

/*========================Libraries=========================*/ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "struct.h"
#include "comfunct.h"
#include "check.h"

/*========================Functions=========================*/

/*========
Display the content of the variable molecule->inputFNList : a list of files names
========*/
void display_fileList(struct MyModel *molecule){
	int i=0;

	printf("dir=%s , nb=%d\n",molecule->inputDir, molecule->nbInputF);
	for(i=0;i<molecule->nbInputF;i++)
		printf("file (%d) %s\n",i, molecule->inputFNList[i]);
}

/*========
Display the covalent bonds based on  molecule->covBond
========*/
void display_cov_bonds(struct MyModel *molecule){
	int i= 1 ;  int k=1;
	printf("Bond \t AtomA-atomB \t Bond type \t distance\n");
	for ( i = 0; i < molecule->nbAtom; i++){ 	
		int j = 0 ;
		for (j = 0; j < molecule->nbAtom; j++){
			if (molecule->covBond[i][j]>0 ||  molecule->covBond[j][i]>0){
				printf("%3d \t %c%2d-%c%2d \t %1d \t %.2lf \n",k, molecule->img[i].atomType , num_occ(molecule->img, i, molecule->nbAtom), molecule->img[j].atomType,num_occ(molecule->img, j, molecule->nbAtom) , molecule->covBond[i][j] , distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z));
				k++;
			}
		}
	}
	printf("\n");
}

/*========
Display the content of a two dimensional matrix of integers
========*/
void display_matrix(int **matrix, int sizeL,int sizeC, char matrix_name[256] ){
	int i;
	int j;

	for(i=0;i<sizeL;i++) for(j=i+1; j<sizeC; j++) if(matrix[i][j]!=-1) printf("%s[%d][%d]=%d\n",matrix_name,i,j,matrix[i][j]);
}

/*========
 Display the content of the vector img. It contains a list of atoms with their chemical types and their (x,y,z) positions
========*/
void display_img(struct atom* img, int size, char img_name[256]){
	int i;

	for(i=0;i<size;i++) printf("%s[%d]:%c, %s, (%.2lf,%.2lf,%.2lf)\n",img_name,i,img[i].atomType,img[i].atomName, img[i].x,img[i].y,img[i].z);
}

/*========
Display the content of the set of characters conf and count the number of '1' in it. 
The variable conf suppose to be a series of 1 and 0 which corresponds to hydrogen bonds. 
At the i-th character, 1 means that the i-th HB is formed in the conformation conf and 0 else.
========*/
void display_conf(char *conf, int size){
	int i=0;
	int cpt=0;
	for(i=0;i<size;i++){
		if(conf[i]=='1') cpt++;
		printf("%c",conf[i] );			
	}
	printf("\t%d\n",cpt);
}

/*========
Display the content of the hash table which contains list of conformation hashed
========*/
void display_tabHash(struct nodHash* tabHash){
	int i; 
	for(i=0;i<(int)pow(2,MAX_PACK);i++){
		printf("conformation of index %d\n",i );
		struct conf_elt *head=tabHash[i].confList;
		while(head){
			printf("%s \t %d \t %d \t %d \n",head->conf, head->index, head->eng, head->pen );
			head = head->suiv;
		}		
	} 
}

/*========
Check and display if there is an atom which has been more than 2 times donor or more than 2 times acceptor. 
One atom should form at most 2 H-bonds as donor, and at most 2 H-bonds as acceptor. 
========*/
void check_hbond_conf(struct confList* conf){
	int i;
	int j;
	int cptA[conf->nbAtom];
	int cptD[conf->nbAtom];

	init_tab(cptD,conf->nbAtom,0);
	init_tab(cptA,conf->nbAtom,0);

	for(i=0;i<conf->nbAtom;i++){
		for(j=i+1;j<conf->nbAtom;j++){ //Update
			if(conf->HB[i][j]==5 && conf->img[i].atomType!='H' && conf->img[j].atomType!='H'){
				//check if there is proton transfer or not
				int D=-1;
				int A=-1;
				check_ptr(conf,i,j,conf->nbAtom,&D);							
				if(D==i){
					A=j;
				}
				if(D==j){
					A=i;
				}
				if(D==-1){
					printf("Error on the donor\n");
					exit(10);
				}	
				cptD[D]++; //Donor
				cptA[A]++; //acceptor
			}
		}					
	}
	printf("check for donors\n");
	for(i=0;i<conf->nbAtom;i++) if(cptD[i]>2) printf("atom (%d) %c%d is %d time as donor\n",i,conf->img[i].atomType, num_occ(conf->img,i , conf->nbAtom), cptD[i] );

	printf("check for acceptors\n");
	for(i=0;i<conf->nbAtom;i++) if(cptA[i]>2) printf("atom (%d)  %c%d is %d time as acceptor\n",i,conf->img[i].atomType, num_occ(conf->img,i , conf->nbAtom), cptA[i] );
}
