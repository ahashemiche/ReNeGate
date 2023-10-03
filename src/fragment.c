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
 * \file fragment.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief Functions to analyse fragments.
 * \details 
 *
 * This file contains all functions used to analyse the fragments of the molecular systems. 
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
#include "memory.h"
#include "comfunct.h"
#include "check.h"
#include "fragment.h"
 
/*========================Functions=========================*/

/*========
Get the list of fragments that have been observed along the trajectory
========*/
void get_frag_dynamics(struct MyModel* molecule){
	//Get the list of fragments 	
	struct confList* head= molecule->conformations;
	molecule->nbFrg =0;
	molecule->fragments=NULL;
	struct frgList *lastfrag=NULL;
	while(head){
		int **CC = allocate_matrix(head->nbAtom,head->nbAtom,0,"CC"); //Adjacency matrix
		int **listfrag = allocate_matrix(head->nbAtom,head->nbAtom+1,-1,"listfrag"); 
		int nbfrag=0;
		init_adj_conf_matrix(CC, head->nbAtom, head ,molecule->level,'f');//Initialize adjacency matrix depending on the level (type of bonds taken into account)		
		// if(head->name==1) display_matrix(CC,head->nbAtom,head->nbAtom,);
		nbfrag = get_cc(CC,molecule->nbAtom,true,listfrag);
		head->nbFrg = nbfrag;
		for(int i=0;i<nbfrag;i++){
			//Construct the fragment 
			struct frgList *currfrg =  create_fragment(head,listfrag,i,molecule->nbConf);
			//Get the isomorphism
			struct frgList *frg = exist_frag(currfrg,molecule->fragments);
			//get the canonical form for each fragment. ?
			if(frg==NULL){//New fragment 
				molecule->nbFrg++;
				currfrg->name=molecule->nbFrg;
				currfrg->conftab[head->name-1]=1;
				if(molecule->fragments==NULL) molecule->fragments=currfrg; //first fragments
				else lastfrag->suiv = currfrg  ;
				lastfrag=currfrg;	
				if(molecule->nbInputF==1)	
					save_frg(molecule,currfrg);
			}
			else{ //Existent fragment 
				//Update the number of occurrence and the conf that is 
				frg->nbocc++;
				frg->conftab[head->name-1]=1;
				free_frgList(currfrg,head->nbAtom);
			}
		}
		free_matrix(CC,head->nbAtom);
		free_matrix(listfrag,head->nbAtom);
		head = head->suiv;	
	}
	if(molecule->nbInputF==1) //Single analysis
		save_frg_list(molecule->fragments,molecule->nbFrg,molecule->nbConf,get_FN(molecule->inputDir,molecule->inputFN,"_frg.txt\0"));
	else //Multiple analysis
		save_frg_list(molecule->fragments,molecule->nbFrg,molecule->nbConf,"occ_frag.txt");

}

/*========
Create a fragment based on the fragment with index "index" in the listfrag
	int name ; //!< Label of the fragment, an unique number is assigned to each fragment 
	struct atom* img ; //!< Pointer to the list of atoms of the conformation (molecular system)
	int nbAtom; //!< Number of atoms in the fragment
	int **B; //!< Adjacency matrix of bonds involved in the fragment
	int nbocc; //!< Number of occurrence of the fragment
	int *conftab; //!< List of conformers where the fragment appears 
	struct frgList *suiv; //!<  Pointer to next conformation 

========*/
struct frgList* create_fragment(struct confList* conf,int **listfrag,int index,int nbConf){ 
	struct frgList *newfrg = allocate_frgList(conf->nbAtom,nbConf,conf->nbAtom,"newfrg"); 
	int i=0;
	int j=0;
	//Save the number of atoms 
	newfrg->nbAtom=conf->nbAtom;
	newfrg->size= listfrag[index][conf->nbAtom];
	newfrg->nbocc=1;
	//Save the Cartesian coordinates	
	for (i = 0; i < conf->nbAtom; i++){
		newfrg->img[i].atomType=conf->img[i].atomType;
		newfrg->img[i].atomName[0]=conf->img[i].atomName[0];		
		newfrg->img[i].atomName[1]=conf->img[i].atomName[1];		
		newfrg->img[i].atomName[2]=conf->img[i].atomName[2];		
		newfrg->img[i].x=conf->img[i].x;		
		newfrg->img[i].y=conf->img[i].y;		
		newfrg->img[i].z=conf->img[i].z;		
		newfrg->img[i].bondMax=0;
		if(listfrag[index][i]==1){
			newfrg->img[i].atomIn=true;
		}
		else{
			newfrg->img[i].atomIn=false;			
		}
	}
	//Put the bonds
	//Covalent bond
	for(i=0;i<conf->nbAtom;i++){
		for(j=i+1;j<conf->nbAtom;j++){
			if(listfrag[index][i]==1 && listfrag[index][j]==1){
				if(conf->CB[i][j]>0){//Covalent bond
					newfrg->B[i][j]=1; 
				}
				else{ 
					if(conf->HB[i][j]==-5){//Hydrogen bond
							int D=-1;
							check_ptr(conf,i,j,conf->nbAtom,&D);							
							newfrg->B[i][j]=2; 	
							newfrg->B[D][D]=6;
					}
					else{//Intermolecular or organometallic interactions
						int ion =get_index_ion(conf->IB,conf->nbIon,conf->nbAtom,i);
						if(ion!=-1){ //Check an intermolecular interaction with i 
							if(conf->IB[ion][j]>0){
								newfrg->B[i][j]=3;
							}
						}
						else{
							ion= get_index_ion(conf->IB,conf->nbIon,conf->nbAtom,j);
							if(ion!=-1){ //Check an intermolecular interaction with j
								if(conf->IB[ion][i]>0){
									newfrg->B[i][j]=3;
								}
							}
							else{
								ion =get_index_ion(conf->MB,conf->nbMetal,conf->nbAtom,i);
								if(ion!=-1){ //Check an organometallic interaction with i 
									if(conf->MB[ion][j]>0){
										newfrg->B[i][j]=4;
									}
								}
								else{
									ion= get_index_ion(conf->MB,conf->nbMetal,conf->nbAtom,j);
									if(ion!=-1){ //Check an organometallic interaction with j
										if(conf->MB[ion][i]>0){
											newfrg->B[i][j]=4;
										}
									}
								}
							}							
						}
					}
				}
			}
		}
	}
	newfrg->suiv=NULL;		
	return newfrg;
}

/*========
Check if the currfrg fragment already exists using the isomorphism
========*/
struct frgList* exist_frag(struct frgList* currfrg, struct frgList* fragments){
	struct frgList* head= fragments;
	while(head){
		if(currfrg->size==head->size){//The fragments should have the same number of atoms
			construct_dreadnaut_file_frg(currfrg, head,"graph.dreadnaut");
			if(verif_dreadnaut_isom( "graph.dreadnaut")){
				return head;
			}
		}
		head = head->suiv;
	}
	return NULL;
}

/*========
construct a dreadnaut file for the two conformation graphs, in order to apply the nauty program for the isomorphism.
========*/
void construct_dreadnaut_file_frg(struct frgList* frgA,struct frgList* frgB, char inputFN[]){
	FILE *outputF =fopen(inputFN,"w");
	int i;
	int j;

	//Two graphs , vertices start with 1
	fprintf(outputF,"$=%d\n",1 );
	//Save first graph for nauty test
	fprintf(outputF, "n=%d g \n",frgA->nbAtom);
	//Save the graph from the adjacency matrix CCA	  
	for(i=0;i<frgA->nbAtom;i++){
		fprintf(outputF, "%d: ",i+1 );
		//Save neighboors
		for(j=0;j<frgA->nbAtom;j++){
			if(frgA->B[i][j]>0 || frgA->B[j][i]>0 )
				fprintf(outputF, "%d ",j+1);
		}
		//Go to the next vertex (atom)
		if(i==frgA->nbAtom-1)
			fprintf(outputF, ".\n");
		else
			fprintf(outputF, ";\n");
	}
	//Save the partitions for the first graph
	save_partition(outputF,frgA->img, frgA->nbAtom);

	//Get and save the canonical graph h (and copy it to h') for the first graph
	fprintf(outputF, ">%s.%s\n",inputFN,"tmp" );
	fprintf(outputF, "%s\n","c x @" );

	//Save the canonical form of the first graph
	fprintf(outputF, ">%s.%s\n",inputFN,"canonical1" );
	fprintf(outputF, "%s\n", "b");

	//Save the second graph graph for nauty test
	fprintf(outputF, "n=%d g \n",frgB->nbAtom);

	//Save the graph from the adjacency matrix CCA	  
	for(i=0;i<frgB->nbAtom;i++){
		fprintf(outputF, "%d: ",i+1 );
		//Save neighboors
		for(j=0;j<frgB->nbAtom;j++){
			if(frgB->B[i][j]>0 || frgB->B[j][i]>0)
				fprintf(outputF, "%d ",j+1);
		}
		//Go to the next vertex (atom)
		if(i==frgB->nbAtom-1)
			fprintf(outputF, ".\n");
		else
			fprintf(outputF, ";\n");
	}

	//Save the partitions for the second graph
	save_partition(outputF,frgB->img, frgB->nbAtom);

	//Get and save the canonical graph h for the second graph
	fprintf(outputF, ">%s.%s\n",inputFN,"tmp");
	fprintf(outputF, "%s\n","x" );

	//Save the canonical form of the second graph
	fprintf(outputF, ">%s.%s\n",inputFN,"canonical2");
	fprintf(outputF, "%s\n", "b");

	//Compare the two graphs 
	fprintf(outputF, ">%s.%s\n",inputFN,"isom" );
	fprintf(outputF, "%s\n", "##");

	//Close file
	fclose(outputF);
}


/*========
Save the mixed graphs of fragments 
========*/
void save_frg(struct MyModel* molecule, struct frgList* frg){
	//create the graph of the fragment using graphViz and save the snapshot
	char resFile[512]="";
    	char graphFN[270]="";
	FILE* outputG=NULL; 
	FILE* outputF=NULL;	
	//Get the graph		
	sprintf(graphFN,"%s/%s",molecule->inputDir,"frgGraph.gv");
  	if((outputG = fopen(get_FN(molecule->inputDir,molecule->inputFN,"_frg.gv\0"),"w")) == NULL){ printf("Can\'not open file to save graphs of fragments \n"); exit(-2);}

	sprintf(resFile,"%s%d.xyz",get_FN(molecule->inputDir,molecule->inputFN,"_frag_xyz/coordinates_frg\0"),frg->name);
  	if((outputF = fopen(resFile,"w")) == NULL){ printf("Can\'not open file to save xyz of fragments \n"); exit(-2);}
	
	fprintf(outputG, "digraph G {\n");
	fprintf(outputG, "label=\"Fragment %d (%d) \";\n",frg->name, frg->nbocc);
	fprintf(outputG, "node [style=filled];\n");
	fprintf(outputG, "graph [bgcolor=transparent];\n");
	fprintf(outputG, "node [shape = circle, fontsize=12];\n");

	int i=0;
	int j=0;
	//Save size of fragment in the xyz file
	fprintf(outputF, "%d\n", frg->size);
	fprintf(outputF, "fragment = %d\n",frg->name ); 

	//Draw atoms except hydrogen atoms
	for(i=0;i<frg->nbAtom;i++){
		if(frg->img[i].atomIn){
			draw_atom(outputG,frg->img,frg->B, frg->nbAtom,i);
			//Save the xyz of the atom 
			fprintf(outputF,"%s \t %lf \t %lf \t %lf \n",frg->img[i].atomName, frg->img[i].x, frg->img[i].y, frg->img[i].z);
		}
	}
	//Make the graph according the bonds
	for(i=0;i<frg->nbAtom;i++){
		for(j=i+1;j<frg->nbAtom;j++){
			// if(frg->B[i][j]>0){
			if(frg->B[i][j]>0){
				//Covalent bond
				if(frg->B[i][j]==1 && strcmp(frg->img[i].atomName,"H")!=0 && strcmp(frg->img[j].atomName,"H")!=0){
					fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black];\n",frg->img[i].atomName, num_occ(frg->img, i, frg->nbAtom),frg->img[j].atomName, num_occ(frg->img, j, frg->nbAtom));
				}
				//Hydorgen bond
				if(frg->B[i][j]==2){
					if(frg->B[i][i]==6){
						fprintf(outputG, "\"%s%d\"->\"%s%d\"[color=red , style=dashed];\n",frg->img[i].atomName, num_occ(frg->img, i, frg->nbAtom),frg->img[j].atomName, num_occ(frg->img, j, frg->nbAtom));							
					}
					else{
						fprintf(outputG, "\"%s%d\"->\"%s%d\"[color=red , style=dashed];\n",frg->img[j].atomName, num_occ(frg->img, j, frg->nbAtom),frg->img[i].atomName, num_occ(frg->img, i, frg->nbAtom));							
					}
				}
				//Intermolecular interaction
				if(frg->B[i][j]==3){
					fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = blue, style=dashed];\n",frg->img[i].atomName, num_occ(frg->img, i, frg->nbAtom),frg->img[j].atomName, num_occ(frg->img, j, frg->nbAtom));
					if(strcmp(frg->img[j].atomName,"H")==0 ||  strcmp(frg->img[i].atomName,"H")==0){
						int index; if(strcmp(frg->img[j].atomName,"H")==0) index = j; else index =i;
						int d=get_cov_index(frg->B,frg->nbAtom,index);
						if(d!=-1){
							fprintf(outputG, "\"%s%d\"[fillcolor=gray87, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",frg->img[index].atomName, num_occ(frg->img, index, frg->nbAtom),frg->img[index].atomName);
							fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",frg->img[d].atomName, num_occ(frg->img, d, frg->nbAtom),frg->img[index].atomName, num_occ(frg->img, index, frg->nbAtom) );
						}
						else{
							printf("Error in cov index with hydrogen %d\n",j); exit(-1);
						}
					}				
				}
				//Organometallic interaction
				if(frg->B[i][j]==4){
					fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = deeppink1,style=dashed];\n",frg->img[i].atomName, num_occ(frg->img, i, frg->nbAtom),frg->img[j].atomName, num_occ(frg->img, j, frg->nbAtom));
					if(strcmp(frg->img[j].atomName,"H")==0 ||  strcmp(frg->img[i].atomName,"H")==0){
						int index; if(strcmp(frg->img[j].atomName,"H")==0) index = j; else index =i;
						int d=get_cov_index(frg->B,frg->nbAtom,index);
						if(d!=-1){
							fprintf(outputG, "\"%s%d\"[fillcolor=gray87, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",frg->img[index].atomName, num_occ(frg->img, index, frg->nbAtom),frg->img[index].atomName);
							fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",frg->img[d].atomName, num_occ(frg->img, d, frg->nbAtom),frg->img[index].atomName, num_occ(frg->img, index, frg->nbAtom) );
						}
						else{
							printf("Error in cov index with hydrogen %d\n",j); exit(-1);
						}
					}
				}
			}
		}
	}
	fprintf(outputG, "}\n");
	fclose(outputG);
	fclose(outputF);

	//Draw the graph of transition using graphViz
	strcpy(resFile,"neato -Tpng  ");
	strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_frg.gv -o \0"));
	strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_frag/fragment\0"));
	char a[3];
	sprintf(a,"%d",frg->name);
	strcat(resFile,a); 
	strcat(resFile,".png\0");	
	system(resFile);
}

/*========
Draw an atom in output file 
========*/
void draw_atom(FILE *output, struct atom* img, int **CB, int nbAtom, int index){
	//Oxygen in red
	if(strcmp(img[index].atomName,"O")==0){
		if(get_nb_CB_type(img,index,CB,nbAtom,"Si")==1)
			fprintf(output, "\"%s%d\"[fillcolor=goldenrod1, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);
		else	
			fprintf(output, "\"%s%d\"[fillcolor=red, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);
	}
	//Neitrogen in blue
	if(strcmp(img[index].atomName,"N")==0)
		fprintf(output, "\"%s%d\"[fillcolor=blue, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);
	//Carbon in gray
	if(strcmp(img[index].atomName,"C")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=gray32, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);
	//Lithium in blueviolet
	if(strcmp(img[index].atomName,"Li")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=blueviolet, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);
	//Chlore in green
	if(strcmp(img[index].atomName,"Cl")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=green, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);
	//Sulphur in cyan
	if(strcmp(img[index].atomName,"S")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=cyan, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);
	//Manganese in deeppink1 
	if(strcmp(img[index].atomName,"Mn")==0){
		fprintf(output, "\"%s%d\"[fillcolor=deeppink1, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	} 
	//Potassium in midnightblue 
	if(strcmp(img[index].atomName,"K")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=midnightblue, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	//Boron in lightpink  
	if(strcmp(img[index].atomName,"B")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=lightpink, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	//Brom in lightgreen  
	if(strcmp(img[index].atomName,"Br")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=seagreen1 , fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	//Phosphorus in dodgerblue1  
	if(strcmp(img[index].atomName,"P")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=dodgerblue1 , fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	//Flourine in darkgoldenrod1  
	if(strcmp(img[index].atomName,"F")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=darkgoldenrod1 , fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	//Sodium in  steelblue 
	if(strcmp(img[index].atomName,"Na")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=steelblue, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	//Ruthinium in deeppink3 
	if(strcmp(img[index].atomName,"Ru")==0)
		fprintf(output, "\"%s%d\"[fillcolor=deeppink3, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	//Iodine in chocolate1  
	if(strcmp(img[index].atomName,"I")==0) 
		fprintf(output, "\"%s%d\"[fillcolor=chocolate1, fontcolor=white, fontname=\"blod, bold\", label=\"%s\"];\n",img[index].atomName, num_occ(img, index, nbAtom),img[index].atomName);					
	 	
}

/*========
Save fragments list, their number occurrence, conformers where they appear, etc. 
========*/
void save_frg_list(struct frgList* fragments, int nbFrg, int nbConf, char *outputFN ){
	FILE* outputF=NULL; 	
  	struct frgList* head=fragments ;
  	if((outputF =fopen(outputFN,"w")) == NULL){ printf("Can\'not open file to save fragments \n"); exit(-2);}

  	fprintf(outputF, "%d\n",nbFrg );
  	fprintf(outputF, "Name \t size \t #occ \t conf.\n");
  	while(head){
  		fprintf(outputF, "%d \t %d \t %d \t ",head->name, head->size, head->nbocc );
  		//Save conformers list where the fragment appears
  		fprintf(outputF, "(");
  		for(int i=0;i<nbConf;i++){
  			if(head->conftab[i]==1){
  				fprintf(outputF, "%d, ",i+1);
  			}
  		}
  		fprintf(outputF, ")");
  		fprintf(outputF, "\n");
  		head = head->suiv;
  	}

  	fclose(outputF);
}
/*========
Identify the motifs that occured in an interface based on the directed paths formed by the HBs between W-W , W-S and S-S 
========*/
void get_motifs(struct MyModel *molecule){
	//Get the list of fragments 	
	struct confList* head= molecule->conformations;
	while(head){
		int **CC = allocate_matrix(head->nbAtom,head->nbAtom,0,"CC"); //Adjacency matrix
		//Initialize adjacency matrix depending on the level, and take into account the direction of the hydrogen bond		
		init_adj_conf_matrix(CC, head->nbAtom, head ,molecule->level,'p');
		
		//Search the different paths 
		int **listMotif = allocate_matrix(head->nbAtom,head->nbAtom+1,-1,"list_motif"); 
		// int **list_motif=NULL;
		int nbmotif= get_dir_path(CC, head->nbAtom, listMotif);
		
		for(int i=0;i<nbmotif;i++){
			printf("motif %d : size %d\n",i,listMotif[i][head->nbAtom]);

			struct frgList *currfrg = create_fragment(head,listMotif,i,molecule->nbConf);
			currfrg->name= (head->name*10)+i;
			save_frg(molecule,currfrg);
			free_frgList(currfrg,head->nbAtom);			
		}

		//Realse memory 
		free_matrix(CC,head->nbAtom);
		free_matrix(listMotif,head->nbAtom);

		head =head->suiv;
	}
}
