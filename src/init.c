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
 * \file init.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief All initialization and some trajectory's file treatments.
 * \details This file contains all functions used to initialize the model such as:
 *
 *		- Initializing variables (lists, adjacency matrices, etc.)
 *		- Reading parameters and constants
 *		- Creating of results files and directories
 *		- Computing the trajectory size (number of snapshots)
 *		- Computing the molecular system size (number of atoms)
 *		- Computing the number of Hydrogen atoms
 *		- Computing the number of ion atoms
 *		- Get the chemical formula of the molecular system ( CxHxNxOx..., where x is a number) 
 *		- Read one snapshot from the input file 
 *    - Treating one snapshot
 *    - Treating one trajectory
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
#include <ctype.h>
#include "constants.h"
#include "struct.h"
#include "comfunct.h"
#include "memory.h"
#include "io.h"
#include "covalent.h"
#include "hbond.h"
#include "organometallic.h"
#include "intermolecular.h" 
#include "conformation.h"
#include "check.h"
#include "init.h"

/*========================Functions=========================*/

/*========
Initialize variables :  lists, adjacency matrices , trajectory size, snapshot size, etc. 
This function is also used to create the results directories and files. 
========*/

void init_model(struct MyModel* molecule){
//Save the directory name
  	char inputF[512]; 
  	sprintf(inputF,"%s/%s",molecule->inputDir,molecule->inputFN);  	
//Initialize the reference snapshot    
    molecule->num_imgRef = 0;
    molecule->num_img = 0;
    molecule->OrbitHSize=0;  	
    molecule->nbAtom =get_nbAtom(inputF);
    molecule->nbHAtom =0;
    molecule->nbIon =0;
    molecule->nbMetal =0;    
    molecule->nbAtomChange=false; //we suppose that we have the same number of atoms
//Get the #molecules (size of the trajectory)
  	molecule->nbMolecule =get_nbMolecule(inputF); 
//Initialize the covalent bonds matrix at '0' & Donor-acceptor at -1 &  Intermolecular bonds matrix at -1
  	allocate_matrix_model(molecule,molecule->nbAtom);
   //Organometallic bonds
    molecule->nbAtomDyn=molecule->nbAtom;
    molecule->bondDyn=allocate_matrix(molecule->nbAtom,molecule->nbAtom,-1,"bondDyn");
    

//Initialize the isthmus list at -1
  	molecule->isthList[0]=0;
  	int i=0; 
  	for (i = 1; i < NB_MaxIsthm; i++) molecule->isthList[i]=-1;	
//If global analysis of multiple files
  if(molecule->nbInputF>1){
	//Create folders for analysis
	if(molecule->num_F==0){ 
		//Create a folder that contains a file for each trajectory describing the sequence of conformations explored
		  create_res_folder(molecule->inputDir,"traj_conf");
		//Create a folder that contains a file for each trajectory describing the path of conformations explored
  		create_res_folder(molecule->inputDir,"traj_conf_graph");		
  	//Create a folder that contains a file for each global graph of conformations explored
  		create_res_folder(molecule->inputDir,"global_conf_graph");
		//Create a folder that contains a file for Cartesian coordinates for each conformation explored 
  		create_res_folder(molecule->inputDir,"conf_xyz");
		//Create a folder that contains a file for each conformation explored describing the bonds
  		create_res_folder(molecule->inputDir,"conformations");
		//Create a folder that contains a file for the corresponding graph of each conformation explored 
  		create_res_folder(molecule->inputDir,"conf_graph");
    //Create a folder to put the fragements structures
      create_res_folder(molecule->inputDir,"_frag"); 
    //Create a folder to put the fragements structures
      create_res_folder(molecule->inputDir,"_frag_xyz"); 
  	}
	char outputFN[512];	
	sprintf(outputFN,"%s/%s/%s",molecule->inputDir,"traj_conf",molecule->inputFNList[molecule->num_F]);
	outputFN[strlen(outputFN)-strlen(strrchr(molecule->inputFNList[molecule->num_F],'.'))]='\0';
	strcat(outputFN,".conf"); 
	FILE *conFile = fopen(outputFN,"w");
	fprintf(conFile, "%s\n",molecule->inputFNList[molecule->num_F]);	
	fprintf(conFile, "%d\n",molecule->nbMolecule);
	fprintf(conFile, "%d\n", molecule->nbAtom);
	fclose(conFile);
  }//Single analysis
  else{
  //Create a folder to put results    	
	char resDir[512];
	strcpy(resDir,inputF);
	resDir[strlen(resDir)-strlen(strrchr(resDir,'.'))]='\0';
	//Check if the directory already exists
    DIR *rep=NULL;
    if((rep =opendir(resDir))){ //Directory already exists, so delete it
      closedir(rep);
      del_dir(resDir);
    }
    else{
      closedir(rep);
    }
    mkdir(resDir,0777); 
  //Create a folder to put the covalent bonds
    mkdir(get_FN(molecule->inputDir,molecule->inputFN,"_cov\0"),0777); 
  //Create a folder to put the conformations structures
    mkdir(get_FN(molecule->inputDir,molecule->inputFN,"_conf\0"),0777);   	
  //Create a folder to put the transitional state structures
    mkdir(get_FN(molecule->inputDir,molecule->inputFN,"_trans\0"),0777);     
  //Create a folder to put the xyz files for each conformation explored
    mkdir(get_FN(molecule->inputDir,molecule->inputFN,"_xyz\0"),0777); 
  //Create a folder to put the xyz files for each conformation explored with extra columns 
    mkdir(get_FN(molecule->inputDir,molecule->inputFN,"_2_xyz\0"),0777); 
  //Create a folder to put the fragements structures
    mkdir(get_FN(molecule->inputDir,molecule->inputFN,"_frag\0"),0777); 
  //Create a folder to put the fragements structures
    mkdir(get_FN(molecule->inputDir,molecule->inputFN,"_frag_xyz\0"),0777); 
 //Create a folder that contains a file for each conformation explored describing the bonds
    mkdir(get_FN(molecule->inputDir,molecule->inputFN,"_bonds\0"),0777);
  //Initialize event's file
    save_event(molecule,1,0,'-' ,'w');
 //Create file to save periods of conformer
    sprintf(molecule->periodFile,"%s",get_FN(molecule->inputDir,molecule->inputFN,"_periods.dat\0"));

  }
}

/*========
Compute the size of the molecular system (#atoms , also represents the size of one snapshot)
========*/
int get_nbAtom(char inputFN[]){ 
	FILE *traj; //The Cartesian coordinates file (.xyz)
	int nbAtom ; //The size of molecule
  	
//Open the coordinate file
  	if((traj = fopen(inputFN,"r")) == NULL){ 
    	printf("Cannot open file  %s to get molecule size\n", inputFN);
    	exit(-1);
  	}
//Set the file position to the beginning of the file
	rewind(traj); 
//Get the #atoms 
	if(fscanf(traj,"%d", &nbAtom)!=1){
      	printf("Cannot read number of atoms in trajectory .xyz file  (nbAtom)\n");
      	exit(-1);
   	}
//Close the coordinate file
	fclose(traj);
   	return nbAtom;
} 

/*========
Compute the size of trajectory (#snapshots)
========*/
int get_nbMolecule(char inputFN[]){ 
	FILE *traj; //The Cartesian coordinates of atoms (.xyz)
	int nbMolecule=0 ; // size of trajectory
//Open the coordinate file
  	if((traj = fopen(inputFN,"r")) == NULL){ 
    	printf("Cannot open file %s to get trajectory size\n", inputFN);
    	exit(-1);
  	}
  	int nbAtom;
  	while (fscanf(traj,"%d",&nbAtom) == 1) {
  		nbMolecule++;
  		skip_return(traj,nbAtom+2);
  	}
//Close the coordinate file
	fclose(traj);
  return nbMolecule;
} 

/*========
Get #ion atoms on the molecular system
========*/
int get_ion_atoms(struct MyModel *molecule){
	int i =0; int j=0 ;
	while(i<molecule->nbAtom) {if(strcmp(molecule->img[i].atomName,"Li")==0 ||  strcmp(molecule->img[i].atomName,"Na")==0 || strcmp(molecule->img[i].atomName,"Ar")==0 || strcmp(molecule->img[i].atomName,"K")==0 || strcmp(molecule->img[i].atomName,"I")==0 || strcmp(molecule->img[i].atomName,"Cl")==0 || strcmp(molecule->img[i].atomName,"Br")==0 || strcmp(molecule->img[i].atomName,"F")==0 ) j++; i++; }
	return j ;
}
/*========
Get #metal atoms on the molecular system
========*/
int get_metal_atoms(struct MyModel *molecule){
  int i =0; int j=0 ;
  while(i<molecule->nbAtom) {if(strcmp(molecule->img[i].atomName,"Mn")==0 || strcmp(molecule->img[i].atomName,"Ru")==0 || strcmp(molecule->img[i].atomName,"Au")==0) j++; i++; }
  return j ;
}
/*======
Transform chain into  lower letter 
======*/
void lowchain(char *chain){
    int i = 0;
    for (i = 0 ; chain[i] != '\0' ; i++){
        chain[i] = tolower(chain[i]);
        // chain[i] = toupper(chain[i]);
    }
}
/*========
Read one snapshot (the current one) and save type of atoms with positions.
* COMPLEXITY: O(nbAtom)                            
========*/
void get_molecule(FILE *traj,struct MyModel *molecule){
  int i ;
	double distMax;
	molecule->max1=0;
	molecule->max2=0;
  molecule->nbAtomChange = false;
	int nbA;
	int nbI=0;
  int nbM=0;

	double CMn[3];//Center of Mass (x,y,z) of current image
	char atomType[10]="";
//Initialize the coordinates  
	CMn[0]=0.00;
	CMn[1]=0.00;
	CMn[2]=0.00;
//Check the number of atoms
	fscanf(traj,"%d",&nbA);
	if (nbA != molecule->nbAtom || molecule->sysType[0]=='w') { //if we change number of atoms or it's interface
    //Release memory allocated for the previous snapshot
    free_matrix_model(molecule,molecule->nbAtom);   
    //Allocate a new memory for atom's coordinates with the new number of atoms
    allocate_matrix_model(molecule,nbA);
    molecule->nbAtom = nbA;
    molecule->nbAtomChange = true;
  }
//Skip the two first lines
	skip_return(traj,2);
//Copy the content of the coordinate-file
  // printf("nb atom %d\n",molecule->nbAtom);
  for(i=0; i < molecule->nbAtom; i++){
    // printf("atom %d\n",i);
    //Read type of atom:
    fscanf(traj,"%s",atomType);
    if(strlen(atomType)>=1 && strlen(atomType)<=3){//correct        
    	strcpy(molecule->img[i].atomName,atomType);
      //Check if the first letter of the atom is in capital letter or not. 
      lowchain(molecule->img[i].atomName);
      if (atomType[0]  >= 97 &&  atomType[0] <= 122) atomType[0] = atomType[0] - 32;
      molecule->img[i].atomType=atomType[0];  
      molecule->img[i].atomName[0] = atomType[0];
    }  
    else{
    	printf("atom type incorrect at image %d\n",molecule->num_img );
      exit(-1);        	
    }

    //Read the x,y,z coordinates
  	if(fscanf(traj,"%lf", &molecule->img[i].x)!=1){
  		printf("cannot read the x coordinate at image %d\n", molecule->num_img);
      exit(-1);        	       		
  	}
  	if(fscanf(traj,"%lf", &molecule->img[i].y)!=1){
  		printf("cannot read the y coordinate at image %d\n", molecule->num_img);
  		exit(-1);
  	}
  	if(fscanf(traj,"%lf", &molecule->img[i].z)!=1){
  		printf("cannot read the z coordinate at image %d\n", molecule->num_img);
  		exit(-1);
    }
  	skip_return(traj,1); //go to the next line 
	  //Set the valence of current atom (maximum number of covalent bonds that an atom can form) 
	    molecule->img[i].bondMax = 0 ; 
  	//switch(molecule->img[i].atomType){
			if(strcmp(molecule->img[i].atomName,"C")==0) 
        molecule->img[i].bondMax = 4 ; 

			if(strcmp(molecule->img[i].atomName,"N")==0) 
        molecule->img[i].bondMax = 3 ; 
      
			if(strcmp(molecule->img[i].atomName,"O")==0)  
        molecule->img[i].bondMax = 2 ; 
      
			if(strcmp(molecule->img[i].atomName,"H")==0) 
        molecule->img[i].bondMax = 1 ; 
      
      if(strcmp(molecule->img[i].atomName,"S")==0) 
        molecule->img[i].bondMax = 6 ; 

      if(strcmp(molecule->img[i].atomName,"Si")==0) 
        molecule->img[i].bondMax = 4 ; 
      
      if(strcmp(molecule->img[i].atomName,"Mn")==0) 
        molecule->img[i].bondMax = 7 ; 
      
      if(strcmp(molecule->img[i].atomName,"B")==0) 
        molecule->img[i].bondMax = 3 ; 
      
      if(strcmp(molecule->img[i].atomName,"F")==0) 
        molecule->img[i].bondMax = 1 ; 
      
      if(strcmp(molecule->img[i].atomName,"P")==0) 
        molecule->img[i].bondMax = 5 ; 
      
      if(strcmp(molecule->img[i].atomName,"Ru")==0) 
        molecule->img[i].bondMax = 5 ; 
      
      if(strcmp(molecule->img[i].atomName,"Zn")==0) 
        molecule->img[i].bondMax = 4 ; 

      if(strcmp(molecule->img[i].atomName,"Au")==0) 
        molecule->img[i].bondMax = 0 ; 
			// default  : molecule->img[i].bondMax = 0;

      if(strcmp(molecule->img[i].atomName,"Cl")==0 ){
       			if(molecule->sysType[0]=='w') 
              molecule->img[i].bondMax = 4;
            else
             molecule->img[i].bondMax = 0;

      }
		//}
		if(strcmp(molecule->img[i].atomName,"Br")==0 ||  strcmp(molecule->img[i].atomName,"Na")==0 ||  strcmp(molecule->img[i].atomName,"Li")==0 ||  strcmp(molecule->img[i].atomName,"Ar")==0 ||   strcmp(molecule->img[i].atomName,"K")==0 ||  strcmp(molecule->img[i].atomName,"I")==0 ){
			molecule->img[i].bondMax = 0;
		}

	   //Check if there is an ion
		if( /*(molecule->sysType[0]=='c' || molecule->sysType[0]=='w') && */(molecule->num_img==0 || molecule->nbAtomChange)){  
			//Cation or anion
      if(strcmp(molecule->img[i].atomName,"Li")==0 || strcmp(molecule->img[i].atomName,"Na")==0 || strcmp(molecule->img[i].atomName,"Ar")==0 || strcmp(molecule->img[i].atomName,"K")==0 || strcmp(molecule->img[i].atomName,"I")==0 || strcmp(molecule->img[i].atomName,"F")==0 || strcmp(molecule->img[i].atomName,"Br")==0 || strcmp(molecule->img[i].atomName,"Cl")==0 ){
        molecule->ionBond[nbI][i] = -2 ; // Just to indicate that is not a simple atom
				molecule->ionBond[nbI][molecule->nbAtom] = i ;
				molecule->ionBond[nbI][molecule->nbAtom+1] = 0 ; //Coordinate number (CN) , no bonds for the moment 
				nbI++;
			}
      //Metal
      if(strcmp(molecule->img[i].atomName,"Mn")==0 || strcmp(molecule->img[i].atomName,"Ru")==0 ||  strcmp(molecule->img[i].atomName,"Au")==0 ){
        molecule->metalBond[nbM][i] = -2 ; // Just to indicate that is not a simple atom
        molecule->metalBond[nbM][molecule->nbAtom] = i ;
        molecule->metalBond[nbM][molecule->nbAtom+1] = 0 ; //Coordinate number (CN) , no bonds for the moment         
        nbM++;
      }            
		}		

	  //Calculate the center of mass
		CMn[0]+=molecule->img[i].x;
		CMn[1]+=molecule->img[i].y;
		CMn[2]+=molecule->img[i].z;
	  //Change the (x,y,z) axis
		molecule->img[i].x-=molecule->cm[0];
		molecule->img[i].y-=molecule->cm[1];
		molecule->img[i].z-=molecule->cm[2];			
	  //Calculate the max distance (max1,max2)
		if (molecule->num_img==0 || molecule->nbAtomChange)
			distMax=distance(0,0,0,molecule->img[i].x,molecule->img[i].y,molecule->img[i].z);
		else
			distMax=distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->imgRef[i].x,molecule->imgRef[i].y,molecule->imgRef[i].z);
		if(distMax>molecule->max1){
			molecule->max2=molecule->max1;
			molecule->max1=distMax;
		}
		else{
			if(distMax>molecule->max2) molecule->max2=distMax;		
    }
 	}
  if(molecule->num_img==0 || molecule->nbAtomChange){//To change after
    molecule->nbIon=nbI;     
    molecule->nbMetal=nbM; 
  }   
 	molecule->cm[0]=CMn[0]/molecule->nbAtom;
	molecule->cm[1]=CMn[1]/molecule->nbAtom;
	molecule->cm[2]=CMn[2]/molecule->nbAtom;

	#ifdef DISP_IMG
		display_img(molecule->img, molecule->nbAtom, "img");
	#endif
}

/*========
Change the reference snapshot
========*/
struct atom* get_new_ref(struct atom* img, int size){
	struct atom* newImg = allocate_img(size); // atom : (Type,(x,y,z),bondMax)
	int i=0;
	for(i=0;i<size;i++){
		newImg[i].atomType=img[i].atomType;
		newImg[i].atomName[0]=img[i].atomName[0];		
		newImg[i].atomName[1]=img[i].atomName[1];		
		newImg[i].atomName[2]=img[i].atomName[2];		
		newImg[i].x= img[i].x;
		newImg[i].y=img[i].y;
		newImg[i].z=img[i].z;
		newImg[i].bondMax=img[i].bondMax;
	}
	return newImg;
}

/*========
Treat a snapshot : get bonds dynamics
========*/
bool treat_snapshot(FILE *traj,struct MyModel* molecule){
    bool covChange=false; //Change in covalent bonds
  	bool hbondChange=false; //Change in hydrogen bonds
    bool ionChange=false; //Change in intermolecular bonds
    bool metalChange=false; //Change in organometallic bonds
  //Read Snapshot from .xyz file; construct V and V_H .
    get_molecule(traj,molecule); //  List of Cartesian coordinates
  //Get the number of H atoms on the molecule
    if(molecule->num_img==0 || molecule->nbAtomChange){
      molecule->nbHAtom = get_nbTypeatoms(molecule->img,"H",molecule->nbAtom);  
    }
  //Check the condition of recalculation
    if(molecule->max1+molecule->max2 >=(ALPHA-1)*3.24 || molecule->num_img==0 || molecule->nbAtomChange){//
    //Change the reference snapshot
      if( molecule->num_img != 0 ) free_img(molecule->imgRef);        
      molecule->num_imgRef = molecule->num_img;
      molecule->imgRef = get_new_ref(molecule->img, molecule->nbAtom);
      // save_event(molecule,1,0,'r','a');
    //Recompute the strongest bonds, depends on the system studied. 
      covChange=false;
      if(molecule->sysType[0]=='f' || bit_1(molecule->level,1) || molecule->num_img==0 || molecule->nbAtomChange){
      //Construct E_C by computing covalent bonds
        covChange=get_cov(molecule);    
      //Get the name of molecule if it is the first snapshot
        if(molecule->num_img==0){
          get_molecule_name(molecule);                            
          save_cov_bonds(molecule);
        }
      //Check change in covalent bonds
        // if(covChange){
        //Get the internal bonds : this is used to identify the rotational axes
          // get_inter_bonds(molecule);  
        //Save covalent bonds if there is changes
          // save_cov_bonds(molecule);
        //Get the #fragments : get the number of connected components
          // molecule->nbFrg = get_nb_frg(molecule->covBond,molecule->nbAtom);
        // }      
      }
    //Calculate the orbits of hydrogen (recompute the weaker bonds)   
      get_orbits(molecule,ALPHA);     
    }
  //Recompute the strongest bonds, depends on the system studied :intermolecular interactions
      ionChange=false;
      metalChange=false;
      if(molecule->sysType[0]=='c' || molecule->sysType[0]=='w' || bit_1(molecule->level,2) || bit_1(molecule->level,3) ){
      //Construct A_I / M_I by computing intermolecular interactions
        if(molecule->nbIon>0 && bit_1(molecule->level,2)){
          ionChange=get_ion(molecule); 
        } 
        if(molecule->nbMetal>0 && bit_1(molecule->level,3)){
          metalChange = get_metal_dynamics(molecule);
        }      
      } 
  //Dynamic analysis of hydrogen bonds (proton transfer included); construct E_H 
    if(bit_1(molecule->level,0)){
      hbondChange=get_Hbonds_dynamics(molecule);
    }
    //Return true if there has been a change in bonds.
  	if(covChange || ionChange || hbondChange || metalChange || molecule->num_img==0)
  		return true;
  	else
      return false;
}

/*========
Treat a trajectory : for all snapshots get bonds dynamics.
Browse all snapshots and analyse the conformational dynamics. 
Construct the mixed graphs and identify the transitions between the identified conformations.
========*/
void treat_trajectory(struct MyModel* molecule){
  //Get the current file
  	sprintf(molecule->inputFN,"%s",molecule->inputFNList[molecule->num_F]);
    printf("%d) %s\n", molecule->num_F+1 , molecule->inputFN);
//Create a file to save the periods of appearance of each conf 
    char outputPerdFN[strlen(molecule->inputDir)+strlen("conf_periods.txt")+2];
    sprintf(outputPerdFN,"%s/%s",molecule->inputDir,"conf_periods.txt");    
    FILE *outputPerdF;
    if((outputPerdF = fopen(outputPerdFN,"a")) == NULL){ printf("Can\'not create file %s\n", outputPerdFN); exit(3); }

  //Initialization of the model for the current trajectory
  	init_model(molecule);
  //Conformational analysis of current trajectory
  	//open the coordinate file
    FILE *traj;
    char inputF[512]; 
    sprintf(inputF,"%s/%s",molecule->inputDir,molecule->inputFN);
    if((traj = fopen(inputF,"r")) == NULL){ printf("Can\'not open file %s\n", inputF); exit(3); }
    struct confList* pred=NULL;
    //Save the evolution of conformations: 
    // int currConf ;     
    // int numImg =0;
    while(molecule->num_img< molecule->nbMolecule){
    //Get the conformational dynamics (Isomorphism)
      if(treat_snapshot(traj,molecule)) {
        // if(pred)  currConf= pred->name;
        molecule->conformations = add_conf(molecule,&pred,molecule->num_img); 
        // if(currConf != pred->name){
        //   fprintf(outputPerdF, "%3d \t %5d \t %5d \t %3d \t %s\n", molecule->num_img - numImg, numImg,molecule->num_img , currConf, molecule->inputFNList[molecule->num_F]);
        //   numImg =  molecule->num_img;
        // }
      }
    //Go to next snapshot

      molecule->num_img++;
    //Tracking 
      if(molecule->num_img%1000==0){printf("\r %5d / %d ", molecule->num_img,molecule->nbMolecule-1); fflush(stdout);} 
    }      
  //Go to the next file
    fclose(outputPerdF);
    fclose(traj);	
}

/*========
Get the Empirical formula according the atoms of the molecular system
========*/
void get_molecule_name(struct MyModel* molecule){
  char name[256]=""; //The name of the molecule
	strcat(name,add_atom("C",get_nbTypeatoms(molecule->img,"C",molecule->nbAtom)));
	strcat(name,add_atom("H",get_nbTypeatoms(molecule->img,"H",molecule->nbAtom)));
  strcat(name,add_atom("N",get_nbTypeatoms(molecule->img,"N",molecule->nbAtom)));
	strcat(name,add_atom("O",get_nbTypeatoms(molecule->img,"O",molecule->nbAtom)));
	strcat(name,add_atom("Cl",get_nbTypeatoms(molecule->img,"Cl",molecule->nbAtom)));	
  strcat(name,add_atom("Br",get_nbTypeatoms(molecule->img,"Br",molecule->nbAtom)));
	strcat(name,add_atom("Li",get_nbTypeatoms(molecule->img,"Li",molecule->nbAtom)));
	strcat(name,add_atom("Ar",get_nbTypeatoms(molecule->img,"Ar",molecule->nbAtom)));
  strcat(name,add_atom("S",get_nbTypeatoms(molecule->img,"S",molecule->nbAtom)));
  strcat(name,add_atom("Mn",get_nbTypeatoms(molecule->img,"Mn",molecule->nbAtom)));
  strcat(name,add_atom("B",get_nbTypeatoms(molecule->img,"B",molecule->nbAtom)));
  strcat(name,add_atom("K",get_nbTypeatoms(molecule->img,"K",molecule->nbAtom)));
  strcat(name,add_atom("F",get_nbTypeatoms(molecule->img,"F",molecule->nbAtom)));
  strcat(name,add_atom("P",get_nbTypeatoms(molecule->img,"P",molecule->nbAtom)));
  strcat(name,add_atom("Ru",get_nbTypeatoms(molecule->img,"Ru",molecule->nbAtom)));
  strcat(name,add_atom("Na",get_nbTypeatoms(molecule->img,"Na",molecule->nbAtom)));
  strcat(name,add_atom("Si",get_nbTypeatoms(molecule->img,"Si",molecule->nbAtom)));
  strcat(name,add_atom("I",get_nbTypeatoms(molecule->img,"I",molecule->nbAtom)));
  strcat(name,add_atom("Au",get_nbTypeatoms(molecule->img,"Au",molecule->nbAtom)));
  strcat(name,add_atom("Zn",get_nbTypeatoms(molecule->img,"Zn",molecule->nbAtom)));
  sprintf(molecule->name,"%s",name);
}

/*========
Create chain with name of atom and its number of occurrences.
========*/
char* add_atom(char *atomName, int nbAtom){
	static char chaine[10]="";
	if(nbAtom==0) return "";
	if(nbAtom==1)
	 	sprintf(chaine,"%s",atomName);
	else
	 	sprintf(chaine,"%s%d",atomName,nbAtom);		
	return chaine;
}

