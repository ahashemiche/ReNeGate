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
 * \file memory.c 
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief Memory management
 * \details 
 *
 * This file contains functions used for memory management.
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
#include "check.h"
#include "memory.h"

/*========================Functions=========================*/

/*========
Allocate memory for a 2 dimensional matrix of integer numbers with sizeL rows and sizeC columns.
========*/
int**  allocate_matrix( int sizeL,int sizeC, int init_val, char matrix_name[256]){
	//Allocation of rows
	int **matrix = malloc(sizeL*sizeof(*matrix));
	if (!matrix) {fprintf(stderr,"cannot malloc filed for matrix %s \n",matrix_name); exit(1);}
	int i;
	//Columns
	for (i=0 ; i<sizeL ; i++) {
		matrix[i] = malloc(sizeC*sizeof(int));
		if (!matrix[i]) {fprintf(stderr,"cannot malloc filed for line %d of matrix %s\n",i, matrix_name); exit(1);}
	}
	//Initialization
	int j; 
	for(i=0;i<sizeL;i++) for(j=0;j<sizeC; j++) matrix[i][j]=init_val;

	return matrix;
}

/*========  
Release the memory allocated to a matrix of integer numbers.
========*/
void free_matrix(int **matrix, int size){
	int i;
	for (i=0 ; i<size ; i++) free(matrix[i]);
	free(matrix);
}

/*======== 
Allocate memory for a 2 dimensional matrix of char with sizeL rows and sizeC columns.
========*/
char**  allocate_char_matrix(int sizeL,int sizeC, char matrix_name[256]){
	//Allocation 
	char **matrix = malloc(sizeL*sizeof(*matrix));
	if (!matrix) {fprintf(stderr,"cannot malloc filed for matrix %s \n",matrix_name); exit(1);}
	int i;
	for (i=0 ; i<sizeL ; i++) {
		matrix[i] = malloc(sizeC*sizeof(char));
		if (!matrix[i]) {fprintf(stderr,"cannot malloc filed for line %d of matrix %s\n",i, matrix_name); exit(1);}
	}
	return matrix ;
}

/*======== 
Release the memory allocated to a matrix of char.
========*/
void free_char_matrix(char **matrix, int size){
	int i;
	for (i=0 ; i<size ; i++) free(matrix[i]);
	free(matrix);
}

/*======== 
Allocate memory for a vector of type int 
========*/
int* allocate_int_tab(int size, int init_val){
  int *tab=malloc((size)*sizeof(int));
  if (!tab) {fprintf(stderr,"cannot malloc filed for a table\n"); exit(-4);}     
  for (int i = 0; i < size; i++){
    tab[i] = init_val;
  }
  return tab;
}

/*======== 
Release the memory allocated to a vector of type int
========*/
void free_int_tab(int *tab){
  free(tab);
}

/*========
Allocate memory for a vector of struct atom (it represents one snapshot).
========*/
struct atom* allocate_img(int size){
	struct atom* img = (struct atom*)malloc((size)*sizeof(struct atom)); // atoms (Type,Name,(x,y,z),bondMax)
	if (!img) {fprintf(stderr,"cannot malloc filed for  image  \n"); exit(1);}	
	return img;
}

/*========
Release the memory allocated to vector of struct atom.
========*/
void free_img(struct atom* img){
	free(img);
}

/*========
Allocate memory for a struct confList.
========*/
struct confList* allocate_confList( int sizeL, int sizeM, int sizeC,char* confList_name){
	struct confList* newConf =malloc(sizeof(struct confList));
	if (!newConf) {fprintf(stderr,"cannot malloc filed for  conformation node %s  \n",confList_name); exit(1);}	
	//Snapshot
	newConf->img=allocate_img(sizeC);
	//Covalent bonds
	newConf->CB=allocate_matrix(sizeC,sizeC,0,"CB");
	//H-bonds 
	newConf->HB=allocate_matrix(sizeC,sizeC,-1,"HB");
	//Intermolecular bonds
	newConf->IB=allocate_matrix(sizeL,sizeC+2,0,"IB");
	//Organometallic bonds
	newConf->MB=allocate_matrix(sizeM,sizeC+2,0,"MB");

	return newConf;
}

/*========
Release memory allocated for a struct confList.
========*/
void free_confList(struct confList* conf, int sizeL, int sizeM, int sizeC){
	//Snapshot
	free_img(conf->img);
	//Covalent bonds
	free_matrix(conf->CB,sizeC);
	//H-bonds
	free_matrix(conf->HB,sizeC);
	//Intermolecular bonds
	free_matrix(conf->IB,sizeL);
	//Organometallic bonds
	free_matrix(conf->MB,sizeM);
	
	free(conf);
}

/*========
Allocate memory for a struct confList.
========*/
struct frgList* allocate_frgList( int sizeL, int sizeM, int sizeC,char* frgList_name){
	struct frgList* newfrg =malloc(sizeof(struct frgList));
	if (!newfrg) {fprintf(stderr,"cannot malloc filed for  conformation node %s  \n",frgList_name); exit(1);}	
	//Snapshot
	newfrg->img=allocate_img(sizeC);
	//Bonds
	newfrg->B=allocate_matrix(sizeL,sizeC,-1,"B");
	//Conformers
	newfrg->conftab=allocate_int_tab(sizeM,-1);

	return newfrg;
}

/*========
Release memory allocated for a struct frgList.
========*/
void free_frgList(struct frgList* frg, int sizeL){
	//Snapshot
	free_img(frg->img);
	//Bonds
	free_matrix(frg->B,sizeL);
	//Conformers
	free(frg->conftab);	

	free(frg);
}
/*========
Allocate memory for a struct MyModel (global struct that contains all variables and matrices needed in the program).
========*/
void allocate_matrix_model(struct MyModel *molecule,int size){
  	//Allocate memory for atom's coordinates 
    molecule->img = allocate_img(size); // atoms (Type,(x,y,z),bondMax)    	
	//Covalent bonds
  	molecule->covBond=allocate_matrix(size,size,0,"covBond");
  	//Internal covalent bonds
  	molecule->interBond=allocate_matrix(size,size,0,"interBond");
  	//H-bonds
  	molecule->donHacc=allocate_matrix(size,size,-1,"donHacc");
  	//Intermolecular bonds
  	molecule->ionBond=allocate_matrix(size,size+2,-1,"ionBond");
  	//Organometallic bonds
  	molecule->metalBond=allocate_matrix(size,size+2,-1,"metalBond");
  	
  	//Checks
  	#ifdef VERIF_VAL
  		display_matrix(molecule->covBond,size,size, "covBond" );
  		display_matrix(molecule->covBond,size,size, "interBond" );
  		display_matrix(molecule->covBond,size,size, "donHacc" );
  		display_matrix(molecule->covBond,size,size, "ionBond" );
  		display_matrix(molecule->covBond,size,size, "metalBond" );
  	#endif
}

/*========
Allocate memory for a struct tabhash , hash table that stores conformation using hash function
========*/
struct nodHash*  allocate_tabHash(){
	struct nodHash* tabHash = (struct nodHash*)malloc((pow(2,MAX_PACK))*sizeof(struct nodHash)); ////Store index of conformations (tabConf) using hash function
	#ifdef VERIF_MEM_ACCESS
		if(!tabHash){fprintf(stderr,"Cannot malloc filed for tabHash\n"); exit(1);}
	#endif	
	int i; for(i=0;i<(int)pow(2,MAX_PACK);i++) tabHash[i].confList=NULL;	
	return tabHash;	
}

/*========
Release memory for a struct MyModel (global struct that contains all variables and matrices needed in the program).
========*/
void free_matrix_model(struct MyModel *molecule,int size){
  free_matrix(molecule->covBond, size);
  free_matrix(molecule->interBond, size);
  free_matrix(molecule->donHacc, size);
  free_matrix(molecule->ionBond, size);	
  free_matrix(molecule->metalBond, size);	
  free_img(molecule->img);  
}

/*========
Release memory allocated to a struct MyModel.
========*/

void free_mymodel(struct MyModel *molecule){	
  free_img(molecule->imgRef);   
  free_matrix_model(molecule, molecule->nbAtom);
  free_matrix(molecule->bondDyn, molecule->nbAtomDyn);	
  if(molecule->nbInputF>1) free_char_matrix(molecule->inputFNList,molecule->nbInputF);
  struct confList* headConf=molecule->conformations;
  while( headConf!=NULL){
  	struct confList* confTmp = headConf->suiv ;
  	free_confList(headConf,headConf->nbIon, headConf->nbMetal, headConf->nbAtom);
  	headConf = confTmp ;
  }
  struct frgList* headFrg=molecule->fragments;
  while( headFrg!=NULL){
  	struct frgList* frgTmp = headFrg->suiv ;
  	free_frgList(headFrg,headFrg->nbAtom);
  	headFrg = frgTmp ;
  }
  //free(molecule.conformations);
}

/*========
Release memory allocated to a struct MyGraph.
========*/
void free_graph(struct MyGraph *graph){
	free_matrix(graph->adjMatrix,graph->nbdynm+1);
}	

/*========
Release memory allocated to a struct MyTrajGraph.
========*/
void free_graphList(struct MyTrajGraph *graphTraj){
	free_char_matrix(graphTraj->files, graphTraj->nbInputF);
	free(graphTraj->filegraph);
	free(graphTraj->fState);
	free(graphTraj->lState);
	struct MyGraphList* head=graphTraj->graphList;
	while(head){
		free_graph(head->stateG);
		head = head->suiv;
	}
}

/*========
Get the first element of the LIFO
========*/
int get_elt_LIFO(struct path* F){
	return F->head->index;
}

/*========
Delete the first element of the LIFO
========*/
struct path* delet_elt_LIFO(struct path* F){
	struct path_elt* elt = F->head;
	F->head = elt->suiv ; 
	//if(F->head ==NULL) F->queue = NULL;
	free(elt);
	return F;
}

/*=======
Release memory allocated for a LIFO
========*/
void free_LIFO(struct path *F){
	while(F->head !=NULL){
		get_elt_LIFO(F);
	}
	free(F);
}

/*=======
Get the first element of the queue (FIFO)
=======*/
int get_elt_FIFO(struct path* F){
	struct path_elt* elt = F->head;
	int val = elt->index ;
	F->head = elt->suiv ; 
	if(F->head ==NULL) F->queue = NULL;
	free(elt);
	return val;
}

/*=======
Release memory allocated for a FIFO
========*/
void free_FIFO(struct path *F){

	while(F->head !=NULL){
		get_elt_FIFO(F);
		// #ifdef VERIF_VAL
		// 	printf("%d\n",x );
		// #endif
	}
	free(F);
}

/*========
Release the memory allocated to the list of conformations in "level".
========*/
void free_level(struct conf_elt* level){
	struct conf_elt* head=level;
	while(head){
		struct conf_elt* tmp=head->suiv;
		free(head->conf);
		free(head);
		head=tmp;
	}
}

/*========
Release the memory allocated to the hash table.
========*/
void free_tabHash(struct nodHash* tabHash, int size){
	int i=0;
	for(i=0;i<size;i++){
		free_level(tabHash[i].confList);
	}
	free(tabHash);
}

/*========
Release the allocated memory space for tabConf
========*/
void free_tabConf(struct nodConf* tabConf, int size){
	int i=0;
	for(i=0;i<size;i++)
		free(tabConf[i].conf);
	free(tabConf);
}

/*========
Release the allocated memory space for tabCC
========*/
void free_tabCC(struct nodCC* tabCC, int size){
	int i=0;
	for(i=0;i<=size;i++){
	struct confIndex* conf = tabCC[i].confH;
		while(conf){
			struct confIndex* tmp = conf->suiv;
			free(conf); 
			conf = tmp;
		}		
	}
	free(tabCC);
}

/*========
Release the allocated memory space for listEng
========*/
void free_listEng(struct nodEng *listEng ){
	struct nodEng *head= listEng;
	while(head){
		struct nodEng *tmp=head->suiv;
		struct conf_elt* headConf = head->confList;
		while(headConf){ //Browse all conformations of head 
			struct conf_elt* tmp2=headConf->suiv;
			free(headConf->conf);
			free(headConf);
			headConf =tmp2;
		}
		free(head);
		head=tmp;
	}

}

/*========
Release the allocated memory space for graphPoss
========*/
void free_graphPoss(struct GraphModel* graphPoss){
	free_listEng(graphPoss->listEng);
	free_tabConf(graphPoss->tabConf,graphPoss->nbConfT);
	free_tabCC(graphPoss->tabCC,graphPoss->ccLabel);
	free_tabHash(graphPoss->tabHash,(int)pow(2,MAX_PACK));	
}

/*========
Release the memory allocated for tabhash and list of neighbouring conformations and conformations visited 
========*/
void free_path(struct DijkModel path){
	free_tabHBAtoms(path.tabHBAtoms,path.nbHBonds);
	free_nodDijk(path.neighConfList);
	free_nodDijk(path.visConfList);
	free_tabHash(path.tabHash,(int)pow(2,MAX_PACK));
}

/*========
Release the memory allocated for tabHBAtoms table 
========*/
void free_tabHBAtoms(struct nodHB* tabHBAtoms,int size){
	int i=0;
	for(i=0;i<size;i++){
		//Release memory allocated for atom list
		free(tabHBAtoms[i].tabAtoms);
	}
	free(tabHBAtoms);
}

/*========
Release the memory allocated for list of struct nodDijk
========*/
void free_nodDijk(struct nodDijk *list){
	struct nodDijk *head =list;
	while(head){
		struct nodDijk* tmp=head->suiv;
		if(head->tabPen) free(head->tabPen);
		head->confP->d=NULL;
		free(head);
		head=tmp;
	}
}

