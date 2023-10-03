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
 * \file conformation.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief This file contains all functions used in the algorithm in order to analyse the conformational dynamics of molecules
 * \details  
 *
 * This file contains treatment of function used to analyse the conformational dynamics of molecules.
 * For this analysis, a graph theoretic method is used to check a change in conformations. 
 * In this file, function related to the construction of mixed graphs and graphs of transitions are described.
 * For conformational change, we use two levels of analyses:
 *
 *	- A adjacency matrices comparisons: this consists to take the adjacency matrix of each type of bonds and compare element by element.
 *	- Isomorphism test: the Mackay algorithm (nauty) is used in order to decide if two conformations are already isomorphic or not. 
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
#include "check.h"
#include "conformation.h"

/*========================Functions=========================*/
/*========
Analyse the current conformation (numImg) according to the conformations list (add it to the list, if it is not identified yet)
========*/
struct confList* add_transit(struct MyModel* molecule,struct confList** pred, int numImg){
	struct confList* prev=*pred;
	//Create the current conformation : based on molecule->covBond , molecule->bondHList, molecule->ionBond
	struct confList* courState = create_state(molecule);
	//Compare to the previous one
	int transf=0;
	FILE *outputP = fopen(molecule->periodFile,"a");
	if(!verif_isom(prev,courState,molecule->level,molecule->nbAtom,&transf,molecule->bondDyn, molecule->nbAtomDyn)){ //Gi and Gi-1 are not isomorphic 
		if(bit_1(transf,0) || bit_1(transf,1) || bit_1(transf,2)) molecule->changType= add_1(molecule->changType,0); // H-bond
		if(bit_1(transf,3) || bit_1(transf,4)) molecule->changType=add_1(molecule->changType,1); //Covalent 
		if(bit_1(transf,5) || bit_1(transf,6)) molecule->changType=add_1(molecule->changType,2); // Intermolecular 
		if(bit_1(transf,7) || bit_1(transf,8)) molecule->changType=add_1(molecule->changType,3); // Organometallic 
		
		struct confList* state= verif_state(molecule->conformations,courState,molecule->nbAtom,molecule->level,&transf,molecule->bondDyn, molecule->nbAtomDyn); //Verify if the conformation already exists
		if(state != NULL){ //The conformation already exists
			int j=0;//state->imgList[0];
			if(j*2==NB_PERIOD){
				printf("Maximum number of periods allowed (NB_PERIOD=%d) is reached at snapshot %d \n",NB_PERIOD,molecule->num_img);
				exit(-1);
			}
			state->imgList[j*2+1] =numImg; //New period -> think to change the method for the periods (take whole evolution)
			free_confList(courState,molecule->nbIon,molecule->nbMetal,molecule->nbAtom);
		}
		else{ //New conformation
				state = add_state(molecule,courState,numImg);
		}
		if(prev!=NULL){
			int j=1;//prev->imgList[0]+1;
			if(j*2==NB_PERIOD){
				printf("Maximum number of periods allowed (NB_PERIOD=%d) is reached at snapshot %d \n",NB_PERIOD,molecule->num_img);
				exit(-1);
			}
			prev->imgList[j*2] =numImg-1; //End of period
			fprintf(outputP,"%d \t %d \t %d \n",prev->name,prev->imgList[1],prev->imgList[2]);
			prev->imgList[0]+=prev->imgList[2]-prev->imgList[1]+1; //j;
			//Add the successor
			struct successor* succState = verif_succ(*pred,state->name);
			if(succState==NULL){
				verif_isom(*pred,state,molecule->level,molecule->nbAtom,&transf,molecule->bondDyn, molecule->nbAtomDyn);
				succState=new_succ(*pred,state,transf);
			}
			succState->freq++;
		}
        *pred=state;
    }
	else{
	//Gi and Gi-1 are isomorphic 
		// printf("Gi and Gi-1 isomorphic\n");
		free_confList(courState,molecule->nbIon,molecule->nbMetal,molecule->nbAtom);
	}
	fclose(outputP);
	return molecule->conformations;
}

/*========
Update the list of conformations. Check if each conformation is stable or an intermediate state
========*/
void update_conf_states(struct MyModel* molecule){
//Put the type of conformations found
	struct confList* head = molecule->conformations;
	//Update according to the frequencies of appearance 
	while(head){	
		head->totalPerc= ((double)(head->imgList[0])/(double)(molecule->nbMolecule) )*100.00;
		if(head->totalPerc>molecule->tshVal.pourtMinResT){
		//if(verif_timeRes_conf(head,molecule->tshVal.pourtMinResT,molecule->nbMolecule)){
			head->state ='C'; //Real conformation
		}
		else{
			head->state ='T'; //Transitional (intermediate) state
			molecule->nbTransConf++;
		}
		head = head->suiv;
	}
}

/*========
Verify if the conformation courState already exists in the conformations list (head)
========*/
struct confList* verif_state(struct confList* head, struct confList* courState, int nbAtom, int level, int *transf,int **bondDyn, int nbAtomDyn){
	while( head ) {
		if(verif_isom(head,courState,level,nbAtom,transf, bondDyn, nbAtomDyn)) return head;
		head = head->suiv;
	}
	return head;
}
/*========
Verify if confA and confB are isomorphic or not. Comparison of adjacency matrices and isomorphism tests are performed.
========*/
bool verif_isom(struct confList* confA,struct confList* confB,int level, int nbAtom,int *transf, int **bondDyn, int nbAtomDyn){
	*transf =0;
	if(confA == NULL || confB== NULL)
		return false;
	//check if confA and confB are strongly isomorphic (check the sets of atoms)
	if(confA->nbAtom != confB->nbAtom)
		return false;

	//check if confA and confB are weakly isomorphic (check the sets of H-bonds)
	if(bit_1(level,0)){ //check 1/0 with -1
		int i=0;// add comparaison of donor/acceptor
		int j=0;
		for (i = 0; i < nbAtom; i++){
			if(confA->img[i].atomIn && confB->img[i].atomIn){	
				//Proton transfer
				if(confB->HB[i][i]!=confA->HB[i][i]){
					*transf=add_1(*transf,2);
				}			
				for (j = i+1; j < nbAtom; j++){
					if(confA->img[j].atomIn && confB->img[j].atomIn){	
						// printf("j=%d %d %d\n",j, confB->HB[i][j],confA->HB[i][j] );
						if(confA->HB[i][j]!=confB->HB[i][j]) {
							//New Hydrogen bond
							if(confA->HB[i][j]==-1 && confB->HB[i][j]==-5 && confA->img[i].atomType!='H' &&  confA->img[j].atomType!='H'){
								if(nbAtomDyn==nbAtom) add_bondDyn(bondDyn,i,j,2);
								*transf=add_1(*transf,1);
							}
							//H-bond broken
							if(confA->HB[i][j]==-5 && confB->HB[i][j]==-1){
								if(nbAtomDyn==nbAtom) add_bondDyn(bondDyn,i,j,2);						
								*transf=add_1(*transf,0);
							}
						}
					}			
				}
			}
		}
	}
	//check if confA and confB are strongly isomorphic (check the sets of covalent bonds)
	if(bit_1(level,1)){ //check 1/0 with -1
		int i=0;
		int j=0;
		
		for(i=0;i<nbAtom;i++)
			for(j=i+1;j<nbAtom;j++){
				if(confA->img[i].atomIn && confB->img[i].atomIn && confA->img[j].atomIn && confB->img[j].atomIn){				
					if(confA->CB[i][j]!=confB->CB[i][j] && (confA->CB[i][j]==0 || confB->CB[i][j]==0)){
						if(nbAtomDyn==nbAtom) add_bondDyn(bondDyn,i,j,1);
						if(confA->CB[i][j]> confB->CB[i][j]) *transf=add_1(*transf,4); else *transf=add_1(*transf,3);
					}			
				}
			}
	}
	//check if confA and confB are strongly isomorphic (check the sets of intermolecular bonds)
	if(bit_1(level,2)){ //check 1/0 with -1
		int i=0;
		int j=0;
		// for(i=0; i<confA->nbIon; i++){
		// 	if(confA->IB[i][nbAtom+1]!=confB->IB[i][nbAtom+1]){ //Difference of CN
		// 		if(confA->IB[i][nbAtom+1]> confB->IB[i][nbAtom+1]) *transf=add_1(*transf,6); else	*transf=add_1(*transf,5);
		// 	}				
		// }	
		for(i=0; i<confA->nbIon; i++){
			for(j=0;j<nbAtom;j++){
				if(confA->img[confA->IB[i][nbAtom]].atomIn && confB->img[confB->IB[i][nbAtom]].atomIn && confA->img[j].atomIn && confB->img[j].atomIn){
					if(confA->IB[i][j]!=confB->IB[i][j] || (confA->IB[i][j]>0 || confB->IB[i][j]>0)){ //Repenser à comment presenter l'info !! si on a mm CN mais pas la mm distribution 
						if(nbAtomDyn==nbAtom) add_bondDyn(bondDyn,confA->IB[i][nbAtom],j,3);
						if(confA->IB[i][j]> confB->IB[i][j]) *transf=add_1(*transf,6); else	*transf=add_1(*transf,5);
					}												
				}
			}
		}
	}
	//check if confA and confB are strongly isomorphic (check the sets of organometallic bonds)
	if(bit_1(level,3)){ //check 1/0 with -1
		int i=0;
		int j=0;
		// for(i=0; i<confA->nbMetal; i++){
		// 	if(confA->MB[i][nbAtom+1]!=confB->MB[i][nbAtom+1]){ //Difference of CN
		// 		if(confA->MB[i][nbAtom+1]> confB->MB[i][nbAtom+1]) *transf=add_1(*transf,8); else	*transf=add_1(*transf,7);
		// 	}				
		// }	
		for(i=0; i<confA->nbMetal; i++){
			for(j=0;j<nbAtom;j++){
				if(confA->img[confA->MB[i][nbAtom]].atomIn && confB->img[confB->MB[i][nbAtom]].atomIn && confA->img[j].atomIn && confB->img[j].atomIn){
					if(confA->MB[i][j]!=confB->MB[i][j] && (confA->MB[i][j]>0 || confB->MB[i][j]>0) ){ //Repenser à comment presenter l'info !! si on a mm CN mais pas la mm distribution 
						if(nbAtomDyn==nbAtom) add_bondDyn(bondDyn,confA->MB[i][nbAtom],j,4);
						if(confA->MB[i][j]> confB->MB[i][j]) *transf=add_1(*transf,8); else	*transf=add_1(*transf,7);
					}							
				}
			}
		}
	}

	//Check isomorphism
	if(*transf==0){
		return true;
	}
	else{
		construct_dreadnaut_file(confA,confB, level,"graph.dreadnaut");
		return verif_dreadnaut_isom( "graph.dreadnaut");
		// return false;
	}

	//if(*transf==0) return true; else return false;
}
/*========
Add the bond between atoms index1 and index2 as a dynamic bond
========*/
void add_bondDyn(int **bondDyn, int index1, int index2, int type){
	bondDyn[index1][index2]=type;
	bondDyn[index2][index1]=type;
}

/*========
construct a dreadnaut file for the two conformation graphs, in order to apply the nauty program for the isomorphism.
========*/
void construct_dreadnaut_file(struct confList* confA,struct confList* confB, int level, char inputFN[]){
	FILE *outputF =fopen(inputFN,"w");
	int i;
	int j;

	//Two graphs , vertices start with 1
	fprintf(outputF,"$=%d\n",1 );
	//Save first graph for nauty test
	fprintf(outputF, "n=%d g \n",confA->nbAtom);
	int **CCA=allocate_matrix(confA->nbAtom,confA->nbAtom,0,"adjacencyMat"); //Adjacency matrix
	init_adj_conf_matrix(CCA, confA->nbAtom, confA,level,'i');//Initialize adjacency matrix depending on the level (type of bonds taken into account)		
	//Save the graph from the adjacency matrix CCA	  
	for(i=0;i<confA->nbAtom;i++){
		fprintf(outputF, "%d: ",i+1 );
		//Save neighboors
		for(j=0;j<confA->nbAtom;j++){
			if(CCA[i][j]==1)
				fprintf(outputF, "%d ",j+1);
		}
		//Go to the next vertex (atom)
		if(i==confA->nbAtom-1)
			fprintf(outputF, ".\n");
		else
			fprintf(outputF, ";\n");
	}
	//Save the partitions for the first graph
	save_partition(outputF,confA->img, confA->nbAtom);

	//Get and save the canonical graph h (and copy it to h') for the first graph
	fprintf(outputF, ">%s.%s\n",inputFN,"tmp" );
	fprintf(outputF, "%s\n","c x @" );

	//Save the canonical form of the first graph
	fprintf(outputF, ">%s.%s\n",inputFN,"canonical1" );
	fprintf(outputF, "%s\n", "b");

	//Save the second graph graph for nauty test
	fprintf(outputF, "n=%d g \n",confB->nbAtom);
	int **CCB=allocate_matrix(confB->nbAtom,confB->nbAtom,0,"adjacencyMat"); //Adjacency matrix
	init_adj_conf_matrix(CCB, confB->nbAtom, confB,level,'i');//Initialize adjacency matrix from covalent bonds matrix 		

	//Save the graph from the adjacency matrix CCA	  
	for(i=0;i<confB->nbAtom;i++){
		fprintf(outputF, "%d: ",i+1 );
		//Save neighboors
		for(j=0;j<confB->nbAtom;j++){
			if(CCB[i][j]==1)
				fprintf(outputF, "%d ",j+1);
		}
		//Go to the next vertex (atom)
		if(i==confB->nbAtom-1)
			fprintf(outputF, ".\n");
		else
			fprintf(outputF, ";\n");
	}

	//Save the partitions for the second graph
	save_partition(outputF,confB->img, confB->nbAtom);

	//Get and save the canonical graph h for the second graph
	fprintf(outputF, ">%s.%s\n",inputFN,"tmp");
	fprintf(outputF, "%s\n","x" );

	//Save the canonical form of the second graph
	fprintf(outputF, ">%s.%s\n",inputFN,"canonical2");
	fprintf(outputF, "%s\n", "b");

	//Compare the two graphs 
	fprintf(outputF, ">%s.%s\n",inputFN,"isom" );
	fprintf(outputF, "%s\n", "##");

	//Release memory
	free_matrix(CCA,confA->nbAtom);
	free_matrix(CCB,confB->nbAtom);

	//Close file
	fclose(outputF);
}

/*========
Create a state for the current conformation. Based upon molecule->covBond , molecule->bondHList, molecule->ionBond
========*/
struct confList* create_state(struct MyModel* molecule){ 
	struct confList *newConf =allocate_confList(molecule->nbIon,molecule->nbMetal,molecule->nbAtom,"newConf"); 
	int i=0;
	int j=0;
	//Save the number of atoms 
	newConf->nbAtom=molecule->nbAtom;
	newConf->state='C'; //We consider the conformation as stable conformation
	newConf->engMin=1000.00;
	newConf->snapMin=molecule->num_img;
	//Save the Cartesian coordinates	
	for (i = 0; i < molecule->nbAtom; i++){
		// printf("atom %d\n",i);
		newConf->img[i].atomType=molecule->img[i].atomType;
		newConf->img[i].atomName[0]=molecule->img[i].atomName[0];		
		newConf->img[i].atomName[1]=molecule->img[i].atomName[1];		
		newConf->img[i].atomName[2]=molecule->img[i].atomName[2];		
		newConf->img[i].atomIn=false;
		newConf->img[i].x= molecule->img[i].x;
		newConf->img[i].y=molecule->img[i].y;
		newConf->img[i].z=molecule->img[i].z;

	}
	#ifdef VERIF_CONF
		display_img(newConf->img,molecule->nbAtom,"newconf-img");
	#endif
	//Check if there is a parial conformational analysis (analysis of fragments containing metal or cation or anion atoms, only)
	if(molecule->partAnal){ 
		//browse the fragments and update the atomIn field in img structure
		get_atom_in_conf(molecule,newConf->img);
	}
	else{
		for(i=0;i<molecule->nbAtom;i++) newConf->img[i].atomIn=true;
	}
	//Put the covalent bonds
	for(i=0;i<molecule->nbAtom;i++){
		if(newConf->img[i].atomIn){
			newConf->CB[i][i]= 0;//molecule->covBond[i][i];
			for(j=i+1;j<molecule->nbAtom;j++){
				newConf->CB[i][j]=molecule->covBond[i][j]; //Get the covalent bond
				newConf->CB[j][i]=0; //we suppose that we haven't an Isthmus (bridge) 
			}
		}
	}
	#ifdef VERIF_CONF
		display_matrix(newConf->CB,molecule->nbAtom,molecule->nbAtom,"newconf-CB");
	#endif

	//Put the H-Bonds
	if(bit_1(molecule->level,0)){
		for(i=0;i<molecule->nbAtom; i++){
			if(newConf->img[i].atomIn){
				newConf->HB[i][i]=molecule->donHacc[i][i];
				for(j=i+1;j<molecule->nbAtom;j++){
					//Add only the donor (heavy atom) and acceptor (D,A) 
					if(molecule->donHacc[i][j]==-5 && molecule->img[i].atomType!='H' && molecule->img[j].atomType!='H' ){
						// printf("1:HB %d , i=%s%d, j=%s%d \n", molecule->num_img, newConf->img[i].atomName, num_occ(newConf->img, i, molecule->nbAtom), newConf->img[j].atomName, num_occ(newConf->img, j, molecule->nbAtom));	
						//Search D than put all HD at 2 
						int D=get_HHB(molecule->donHacc, i,j,molecule->nbAtom);	
						int A ; if(D==i) A=j ; else A=i;
						if(D!=-1){
							int k=0;
							for(k=0;k<molecule->nbAtom;k++){//Put for all hydrogen bond the index of donor
								// printf("D=%d , k=%d , cov=%d , %d\n",D,k,molecule->covBond[k][D],molecule->covBond[D][k] );
								//Proton transfer
								if(k!=D && molecule->img[k].atomType=='H' && (molecule->donHacc[k][k]==D || molecule->donHacc[k][A]==1 || molecule->donHacc[k][A]==0) ){
									if(molecule->donHacc[k][k]==D){//Proton transfer
										newConf->HB[k][k]=D;	
										newConf->HB[i][j]=molecule->donHacc[i][j];
										newConf->HB[j][i]=molecule->donHacc[j][i];	
									}
									else{
										if(molecule->covBond[k][D]==1 || molecule->covBond[D][k]==1){										
												// printf("2:HB %d , i=%s%d, j=%s%d , %.2lf , %.2lf\n", molecule->num_img, newConf->img[i].atomName, num_occ(newConf->img, i, molecule->nbAtom), newConf->img[j].atomName, num_occ(newConf->img, j, molecule->nbAtom),DA, ang);	
											// Simple HB
											newConf->HB[k][D]=D;
											newConf->HB[D][k]=D;								
											if(molecule->donHacc[k][A]==1 || molecule->donHacc[k][A]==0){
												// printf("3:HB %d , i=%s%d, j=%s%d \n", molecule->num_img, newConf->img[i].atomName, num_occ(newConf->img, i, molecule->nbAtom), newConf->img[j].atomName, num_occ(newConf->img, j, molecule->nbAtom));	
												newConf->HB[A][k]=-2;								
												newConf->HB[k][A]=-2;	
												newConf->HB[i][j]=molecule->donHacc[i][j];
												newConf->HB[j][i]=molecule->donHacc[j][i];	
											}
										}							
									}
								}
							}
						}	
						else{
							printf("problemn img %d , i=%s%d, j=%s%d  %d , %d\n", molecule->num_img, newConf->img[i].atomName, num_occ(newConf->img, i, molecule->nbAtom), newConf->img[j].atomName, num_occ(newConf->img, j, molecule->nbAtom), i,j);	
							exit(10);
						}
					}
				}
			}
		}
		#ifdef VERIF_CONF
			display_matrix(newConf->HB,molecule->nbAtom,molecule->nbAtom,"newconf-HB");
		#endif
	} 

	//Intermolecular bonds
	newConf->nbIon=molecule->nbIon;
	if(bit_1(molecule->level,2)){
		for(i=0;i<molecule->nbIon; i++){
			if(newConf->img[molecule->ionBond[i][molecule->nbAtom]].atomIn){ //Check how to put atomIn
				for(j=0;j<molecule->nbAtom+2;j++){
					newConf->IB[i][j]=molecule->ionBond[i][j];
				}
			}
			else{
				for(j=0;j<molecule->nbAtom;j++){
					newConf->IB[i][j]=-1;
				}
				newConf->IB[i][molecule->nbAtom]=molecule->ionBond[i][molecule->nbAtom];
				newConf->IB[i][molecule->nbAtom+1]=0;
				newConf->IB[i][molecule->ionBond[i][molecule->nbAtom]]=-2;
			}
		}
		#ifdef VERIF_CONF
			display_matrix(newConf->IB,molecule->nbIon,molecule->nbAtom+2,"newconf-IB");
		#endif
	}
	//Organometallic bonds
	newConf->nbMetal=molecule->nbMetal;
	if(bit_1(molecule->level,3)){
		for(i=0;i<molecule->nbMetal; i++){
			if(newConf->img[molecule->metalBond[i][molecule->nbAtom]].atomIn){ //Check how to put atomIn		
				for(j=0;j<molecule->nbAtom+2;j++){
					newConf->MB[i][j]=molecule->metalBond[i][j];
				}
			}
			else{
				for(j=0;j<molecule->nbAtom;j++){
					newConf->MB[i][j]=-1;
				}
				newConf->MB[i][molecule->nbAtom]=molecule->metalBond[i][molecule->nbAtom];
				newConf->MB[i][molecule->nbAtom+1]=0;
				newConf->MB[i][molecule->metalBond[i][molecule->nbAtom]]=-2;
			}
		}
		#ifdef VERIF_CONF
			display_matrix(newConf->MB,molecule->nbMetal,molecule->nbAtom+2,"newconf-MB");
		#endif
	}
	newConf->succ=NULL; 			
	newConf->suiv=NULL;	
	if(molecule->sysType[0]=='w')get_atom_in_interface(newConf);
	return newConf;
}

/*========
Add the new conformation newConf at the "end" of the conformations list (molecule->conformations)
========*/
struct confList* add_state(struct MyModel* molecule,struct confList *newConf,int numImg){
	molecule->nbConf++;
	newConf->name = molecule->nbConf;
	newConf->imgList[0]=0;
	newConf->imgList[1]=numImg;

	//Get rotational parts of molecule (covalent bonds dynamics except leaves)
	newConf->nbFrg = 0;//molecule->nbFrg; //Get the number of fragments
	if(bit_1(molecule->level,4))
		get_bondcovdyn(molecule,&newConf);
	
	if(molecule->conformations==NULL) molecule->conformations=newConf; //first conformation
	else molecule->lastConf->suiv =newConf ;
	molecule->lastConf=newConf;

	return molecule->lastConf;
}

/*========
Verify if a node head has confName as successor
========*/
struct successor* verif_succ(struct confList* head, int confName){
	struct successor* q = head->succ;
	while(q){
		if(q->conf->name==confName)return q;
		q=q->suiv;
	}
	return q;
}

/*========
Create a new successor at the beginning of head->succ list of head.
========*/
struct successor* new_succ(struct confList* head,struct confList* conf,int transf){
	struct successor *newSucc =malloc(sizeof(struct successor));
	#ifdef VERIF_ACCES
		if(!newSucc){fprintf(stderr,"Cannot malloc filed for successor of conformation %d\n",head->name); exit(1);}
	#endif
	newSucc->conf = conf;
	newSucc->freq=0; 
	newSucc->transf=transf;
	newSucc->suiv = head->succ;
	head->succ=newSucc;
	return head->succ;
}
/*========
Get rotational axes of the molecule in conformation conf (cov-bonds dynamics except leaves) 
========*/
void get_bondcovdyn(struct MyModel* molecule,struct confList** conf){
	struct confList* cf=*conf;
//Initialization of covalent bonds
	int i=0;
	int j=0;
	int** CC = allocate_matrix(molecule->nbAtom,molecule->nbAtom,0,"get_bondcovdyn");//Adjacency matrix

//Construct the global Adjacency matrix
	for(i=0;i<molecule->nbAtom;i++){
		CC[i][i]=cf->CB[i][i];
		for(j=i+1;j<molecule->nbAtom;j++){
			if(cf->CB[i][j]>0){ //Internal covalent bond
				CC[i][j]=1;
				CC[j][i]=1; 
			}
			else{ 
				if(cf->HB[i][j]==-5 || cf->HB[i][j]==-5){
					CC[i][j]=2; 
					CC[j][i]=2;  					
				}
				else{//leaves
					CC[i][j]=0; //We don't consider the leaves
					CC[j][i]=0;  
				}
			}
		}
	}
//Get isthmuses (bridges)
	for(i=0;i<molecule->nbAtom;i++){
		for(j=i+1;j<molecule->nbAtom;j++){
			if(molecule->interBond[i][j]==1){ //Internal covalent bond
				CC[i][j]=0; //Delete the edge
				CC[j][i]=0;
				int nbFrg = get_cc(CC,molecule->nbAtom,false,NULL);
				if(nbFrg != cf->nbFrg){
					cf->CB[j][i]=1; //We have an isthmus (bridge) 
					add_isthm(molecule,i,j);
				}
				else{
					cf->CB[j][i]=0; //We don't have an isthmus (bridge)
				}
				CC[i][j]=1; //Reset the edge
				CC[j][i]=1;				
			}		
		}
	}
	free_matrix(CC,molecule->nbAtom);
}

/*========
Add an isthmus to the isthmus list
========*/
void add_isthm(struct MyModel *molecule, int index1, int index2){
	if(verif_exist_isthm(molecule->isthList, index1,index2)==-1){
		int i = molecule->isthList[0];
		if(i==NB_MaxIsthm){
			printf("Maximum number of rotaional axes allowed (NB_MaxIsthm=%d) is reached at snapshot %d \n",NB_MaxIsthm,molecule->num_img);
			exit(-1);
		}							
		molecule->isthList[i*3+1]=index1;
		molecule->isthList[i*3+2]=index2;
		molecule->isthList[i*3+3]=1; //we suppose first that's a simple rotation
		molecule->isthList[0]++;
	}
}

/*========
Verify if the isthmus already exists
========*/
int verif_exist_isthm(int isthList[NB_MaxIsthm], int index1, int index2){
	int i=0;
	while(i<isthList[0])
		if(isthList[i*3+1]==index1 && isthList[i*3+2]==index2)
			return isthList[i*3+3];
		else
			i++;
	return -1;
}

/*========
Update isthmus list to identify the conformational rotation
========*/
bool update_isthList(struct MyModel *molecule){ 
	int i=0;
	bool confRot =false;	
	for(i=0;i<molecule->isthList[0];i++){
		struct confList* head = molecule->conformations;
		while(head){ //For each conformation check if the isthmus (i) already exists, else it means there is a conformational rotation
			if(head->CB[molecule->isthList[i*3+2]][molecule->isthList[i*3+1]]==0){ //isthmus doesn't exist in this conformation
				molecule->isthList[i*3+3]=2; //conformational rotation 
				confRot = true;
				head = NULL;
			}
			else{
				head = head->suiv;				
			} 
		}
	}
	return confRot;
}

/*========
Delete the H-bond that has been never formed (always in position of HB)
Note: this is not used yet. The idea is to optimize the number of conformations according to the H-bonds list.
========*/
void update_HbondList(struct MyModel *molecule, int index){
	//Update conformations
	int D= molecule->bondHdyn[index].bondH[0];
	int H= molecule->bondHdyn[index].bondH[1];
	int A= molecule->bondHdyn[index].bondH[2];
	struct confList* head= molecule->conformations;
	while(head){
		head->HB[A][H]=-1;		
		head->HB[H][A]=-1;		
		head->HB[D][H]=-1;
		head->HB[H][D]=-1;
		head->HB[H][H]=-1;
		// printf("name=%d, HB: %c%d-%c%d-%c%d\n",head->name, molecule->imgRef[D].atomType,num_occ(molecule->imgRef, D, molecule->nbAtom),molecule->imgRef[H].atomType, num_occ(molecule->imgRef, H, molecule->nbAtom),molecule->imgRef[A].atomType, num_occ(molecule->imgRef, A, molecule->nbAtom));
		// int k;int t;
		// for(k=0;k<molecule->nbAtom;k++)
		// 	for(t=k;t<molecule->nbAtom;t++)
		// 		if(head->HB[k][t]!=-1)
		// 			printf("HB[%d][%d]=%d, HB[%d][%d]=%d  \n", k,t,head->HB[k][t], t,k,head->HB[t][k]);
	
		
		if(get_Hbond_sens(head->HB,molecule->nbAtom,A,D)==-1){//There is no H-bond with the same Donor-acceptor
			// printf("on supprime %d %d name=%d, HB: %c%d-%c%d-%c%d\n",D,A,head->name, molecule->imgRef[D].atomType,num_occ(molecule->imgRef, D, molecule->nbAtom),molecule->imgRef[H].atomType, num_occ(molecule->imgRef, H, molecule->nbAtom),molecule->imgRef[A].atomType, num_occ(molecule->imgRef, A, molecule->nbAtom));
			head->HB[D][A]=-1;
			head->HB[A][D]=-1;
		}
		head = head->suiv;
	}

	//Update bondHdyn vectors
	// int i=0;
	// int j=0;
	// for(i=index;i<molecule->nbHbond-1;i++){
	// 	molecule->bondHdyn[i].nChange=molecule->bondHdyn[i+1].nChange;
	// 	molecule->bondHdyn[i].bondH[0]=molecule->bondHdyn[i+1].bondH[0];
	// 	molecule->bondHdyn[i].bondH[1]=molecule->bondHdyn[i+1].bondH[1];
	// 	molecule->bondHdyn[i].bondH[2]=molecule->bondHdyn[i+1].bondH[2];
	// 	for(j=0;j<molecule->bondHdyn[i].nChange;j++){
	// 		molecule->bondHdyn[i].imgList[j*3]=molecule->bondHdyn[i+1].imgList[j*3];
	// 		molecule->bondHdyn[i].imgList[(j*3)+1]=molecule->bondHdyn[i+1].imgList[(j*3)+1];
	// 		molecule->bondHdyn[i].imgList[(j*3)+2]=molecule->bondHdyn[i+1].imgList[(j*3)+2];
	// 	}
	// }
	// molecule->nbHbond--;
}

/*========
Update the list of conformations if there are similar one, once the H-bond were updated (some H-bonds may be deleted)
Note: this is not used yet. The idea is to optimize the number of conformations according to the H-bonds list.
========*/
void update_ConfList(struct MyModel *molecule){
	struct confList* head = molecule->conformations;

	//Update frequencies 
	while(head){	
		struct successor* succHead = head->succ;
		while(succHead){
			succHead->freq = 0;
			succHead = succHead->suiv;
		}
		head = head->suiv;
	}
	//Update links between conformations
	head = molecule->conformations;
	int i=0;
	int b = head->imgList[2];	
	while(b < molecule->nbMolecule-1){ //While not the end of the trajectory
		int j=-1;
		struct confList* suivHead = get_suiv(molecule,b+1,&j);							
		struct successor* succHead= verif_succ(head,suivHead->name);
		int transf=0;
		if(head->name !=suivHead->name && verif_isom(head,suivHead,molecule->level,molecule->nbAtom,&transf,molecule->bondDyn, molecule->nbAtomDyn)){
			if(succHead !=NULL) succHead->transf=0;  
			if(suivHead->name == molecule->lastConf->name){
				molecule->lastConf = head;
			}
			head->imgList[i*2+2]=suivHead->imgList[j*2+2];
			suivHead = shift_periods(suivHead,j); // Delete [c,d] from suivHead
			if(suivHead->imgList[0]==0){
				suivHead->name = -1;
			}
			else{ 
				if(head->name!=suivHead->name){ //Transfer periods 
					merge_period(head,suivHead);
					suivHead->name = -1;
				}
			}
			// if(i< head->imgList[0]-1 && head->imgList[i*2+2]+1 == head->imgList[i*2+3]){ //[a,b] next [b+1,c] => [a,c]
			// 	head->imgList[i*2+2]=head->imgList[i*2+4];
			// 	head = shift_periods(head,i+1); 
			// } 
		}
		else{
			struct confList* state= verif_state(molecule->conformations,suivHead,molecule->nbAtom,molecule->level,&transf, molecule->bondDyn, molecule->nbAtomDyn); //Verify if the conformation already exists			
			if(state->name==suivHead->name || state->name==-1){
				if(succHead==NULL){
					verif_isom(head,suivHead,molecule->level,molecule->nbAtom,&transf, molecule->bondDyn, molecule->nbAtomDyn);
					succHead=new_succ(head,suivHead,transf);
					succHead->freq=1;//+=succA->freq;			
				}
				else{
					succHead->freq++;
				}
			}
			else{ //There is a conformation isomorphic to suivhead
				if(state->name != suivHead->name && state->name != -1){ //update the conformation list
					merge_period(state,suivHead);	
					if(suivHead->name == molecule->lastConf->name){
						molecule->lastConf = state;
					}					
					suivHead->name = -1;
					struct successor* succHead2= verif_succ(head,state->name);
					if(succHead2==NULL){
						verif_isom(head,state,molecule->level,molecule->nbAtom,&transf,molecule->bondDyn, molecule->nbAtomDyn);
						succHead2=new_succ(head,state,transf);
						succHead2->freq =1; //+= succHead->freq; //1;//+=succA->freq;	
					}
					else{
						succHead2->freq ++; //= succHead->freq;//++;
					}
					j=get_snap_index(state,b+1);
					suivHead =state;
				}
			}
			head = suivHead;
			i=j ;			
		}
		b = head->imgList[i*2+2];
	} 

	head=molecule->conformations; //Deallocate memory !!!
	while(head){
		struct successor* succHead = head->succ;
		struct successor* prev = NULL;
		while(succHead){
			if(succHead->conf->name==-1 || succHead->transf ==0){
				struct successor* emptySucc = succHead->suiv;
				if(prev==NULL)
					head->succ = succHead->suiv;
				else
				 	prev->suiv=succHead->suiv;
				free(succHead);
				succHead=emptySucc;
			}
			else{
				int transf=0;
				verif_isom(head,succHead->conf,molecule->level,molecule->nbAtom,&transf,molecule->bondDyn, molecule->nbAtomDyn);
				succHead->transf=transf;									
				prev = succHead;
				succHead=succHead->suiv;
			}
		}
		head=head->suiv;
	}
	//Update conformation name 
	head = molecule->conformations;
	struct confList* prevHead = NULL;
	i=1;
	while(head){
		if(head->name != -1 && head->imgList[0]!=0){
			int j=0;
			// printf("conformation=%d\n",head->name );
			while(j<head->imgList[0]){
		  		//printf("period %d = %d \n",i,head->imgList[j*2+2]-head->imgList[j*2+1]+1);
				j++;
			}	
			head->name = i;
			// if(head->name==6 || head->name==7){
			// 	printf("conf=%d\n",head->name );
			// 	int k;int t;
			// 	for(k=0;k<molecule->nbAtom;k++)
			// 		for(t=k;t<molecule->nbAtom;t++)
			// 			if(head->HB[k][t]!=-1)
			// 				printf("HB[%d][%d]=%d, HB[%d][%d]=%d \n", k,t,head->HB[k][t], t,k,head->HB[t][k]);
			// }
			i++;
			prevHead = head;
			head = head->suiv;						
		}
		else{ //head->name =-1 , we delete link
			struct confList* emptyHead= head->suiv;
			if(prevHead==NULL)
				molecule->conformations = head->suiv;
			else
			 	prevHead->suiv=head->suiv;
			//free(head);
			free_confList(head,molecule->nbIon,molecule->nbMetal,molecule->nbAtom);
			head = emptyHead;			
		}
	}
	molecule->nbConf = i-1;
}

/*========
Merge all periods of confA in confB
========*/
void merge_period(struct confList* confA,struct confList* confB){
	int i =0;
	int j=0;
	for(i=0; i<confB->imgList[0];i++){
		while(j<confA->imgList[0] && confB->imgList[i*2+2]>confA->imgList[j*2+1]) j++;
		add_period(confA,j);
		confA->imgList[j*2+1]=confB->imgList[i*2+1];
		confA->imgList[j*2+2]=confB->imgList[i*2+2];
	}
	confB->imgList[0]=0;
}

/*========
Delete succ from the list of successors of Head
========*/
void shift_succ( struct confList* head, struct successor* succ){
	struct successor* succHead = head->succ;
	struct successor* prev=NULL;
	while(succHead && succHead->conf->name != succ->conf->name){
		prev = succHead;
		succHead = succHead->suiv;
	}
	if(succHead==NULL){
		printf("ERROR\n");
	}
	else{
		if(prev==NULL)
			head->succ = succHead->suiv;
		else
		 	prev->suiv=succHead->suiv;
	}
}

/*========
Find the conformation that appears at snapshot b+1 
========*/
struct confList* get_suiv(struct MyModel *molecule,int numImg,int *index){
 struct confList* head = molecule->conformations;
 while(head){
 	int i=0;
 	while(i< head->imgList[0] && head->imgList[i*2+1]!=numImg) i++;
	 	if(i<head->imgList[0] && head->imgList[i*2+1]==numImg){
	 		*index = i ;
	 		return head;
	 	}
	 	else
	 		head = head->suiv;
 }
 return NULL;
}

/*========
Find in imgList of conformation conf the index of snapshot b+1
========*/
int get_snap_index(struct confList *conf,int numImg){
 	int i=0;
 	while(i< conf->imgList[0] ){
 		if(conf->imgList[i*2+1]==numImg)
 			return i;
 		else
 			i++;
 	} 
 	return -1;
}

/*========
Delete the period of index in imgList of conformation conf.
========*/
struct confList*  shift_periods(struct confList *conf,int index){
	int i=index;
	for(i=index;i<conf->imgList[0]-1;i++){
		conf->imgList[i*2+1]=conf->imgList[(i+1)*2+1];
		conf->imgList[i*2+2]=conf->imgList[(i+1)*2+2];
	}
	conf->imgList[0]--;
	return conf;
}

/*========
add a period in imgList of conformation conf, at index 
========*/
struct confList*  add_period(struct confList *conf,int index){
	int i=0;
	for(i=conf->imgList[0];i>index;i--){
		conf->imgList[i*2+1]=conf->imgList[(i-1)*2+1];
		conf->imgList[i*2+2]=conf->imgList[(i-1)*2+2];
	}
	conf->imgList[0]++;
	return conf;
}

/*=======================PART (II) : Analysis of multiple trajectories================*/

/*========
Add a conformation to the global conformation list
========*/
struct confList* add_conf(struct MyModel* molecule,struct confList** pred, int numImg){
	struct confList* prev=*pred;
	//Create the current conformation G
	struct confList* courState = create_state(molecule);
	//Compare to the previous one
	int transf=0;
	if(!verif_isom(prev,courState,molecule->level,molecule->nbAtom,&transf,molecule->bondDyn, molecule->nbAtomDyn)){ //Gi and Gi-1 are not isomorphic 
		if(bit_1(transf,0) || bit_1(transf,1) || bit_1(transf,2)) molecule->changType= add_1(molecule->changType,0); // H-bond
		if(bit_1(transf,3) || bit_1(transf,4)) molecule->changType=add_1(molecule->changType,1); //Covalent 
		if(bit_1(transf,5) || bit_1(transf,6)) molecule->changType=add_1(molecule->changType,2); // Intermolecular 
		if(bit_1(transf,7) || bit_1(transf,8)) molecule->changType=add_1(molecule->changType,3); // Organometallic 

		struct confList* state= verif_state(molecule->conformations,courState,molecule->nbAtom,molecule->level,&transf,molecule->bondDyn, molecule->nbAtomDyn); //Verify if the conformation already exists
		if(state == NULL){ //New conformation
			state = add_state(molecule,courState,numImg);
		}
		else
			free_confList(courState,molecule->nbIon,molecule->nbMetal,molecule->nbAtom);
		//Add in the file the new current conformation
	  	char outputFN[270];
	  	sprintf(outputFN,"%s%s",molecule->inputDir,"/traj_conf");
	  	strcat(outputFN,"/");
	  	strcat(outputFN,molecule->inputFNList[molecule->num_F]);
	  	outputFN[strlen(outputFN)-strlen(strchr(molecule->inputFNList[molecule->num_F],'.'))]='\0';
	  	strcat(outputFN,".conf");
	  	FILE *conFile = fopen(outputFN,"a");
	  	fprintf(conFile, "%d \t %d\n", numImg , state->name);
	  	fclose(conFile);			
		//Change the previous conformation
		*pred=state;	
	}
	else
		free_confList(courState,molecule->nbIon,molecule->nbMetal,molecule->nbAtom);
	return molecule->conformations;
}
