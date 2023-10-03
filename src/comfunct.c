/*****************************************************************************
    ReNeGate

    Copyright (c) 2014-2022 Sana Bougueroua
                  2020-2022 Ali Hashemi
    Please cite:  J. Chem. Phys. 2018, 149 (18), 184102.         (DOI 10.1063/1.5045818 )
		   	
    

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
 * \file comfunct.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief Common functions. 
 * \details 
 *
 * This file contains all common functions used by many parts of the algorithm.
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
 
/*========================Functions=========================*/

/*========
Skip n '\n' in file F
========*/
void skip_return (FILE *F, int n) {
	while (n) {
		while (fgetc(F) != '\n')  ;
		n--;
	}
}

/*========
Skip nb_snapshot if possible  
========*/
void skip_snapshot(FILE* traj, int num_snap, int nb_snaps){
	int cpt_snap=0;
	int size_snap;
	if(num_snap+JUMP_SNAPS<nb_snaps){
		while(cpt_snap<JUMP_SNAPS ){
			fscanf(traj,"%d",&size_snap);
			skip_return(traj,size_snap+2);
			cpt_snap++;
		}
	}
	else{
		printf("Arrived to the end of the trajectory, no suffisient snapshots to jump\n");
	}
}

/*========
Create folder_name directory in path
========*/
void create_res_folder(char path[500], char folder_name[256]){
	char outputFN[756];	
	sprintf(outputFN,"%s/%s",path,folder_name);
	mkdir(outputFN,0777);
}

/*========
Get File name with path
========*/
char* get_FN(char pathFN[], char inputFN[],char* ext){
	static char outputFN[256]="";
	char resFile[256];
	sprintf(outputFN,"%s/%s",pathFN,inputFN);
	outputFN[strlen(outputFN)-strlen(strrchr(inputFN,'.'))]='\0';
	strcpy(resFile,inputFN);
	resFile[strlen(inputFN)-strlen(strrchr(inputFN,'.'))]='\0';
	strcat(resFile,ext);
	strcat(outputFN,"/\0");
	strcat(outputFN,resFile);
	//Return the new file name with extension ext
	return outputFN;
}

/*========
Get the number of files in directory filePath
========*/
int get_nbFiles(char filPath[500]){
	int nbFiles=0;
  	DIR* rep =opendir(filPath); 
  	if(rep == NULL){ printf("Can\'not open directory %s\n", filPath); exit(5); }; 	
  	struct dirent* curFile = NULL;   
  	while ((curFile = readdir(rep)) != NULL){
		//All content of rep expect '.' and '..' is considered 
  		if (strcmp(curFile->d_name, ".") != 0 && strcmp(curFile->d_name, "..") != 0){
			nbFiles++;
		}
	}
	//Close the directory
	closedir(rep);
	//Return the number of files found in filPath		
	return nbFiles;
}

/*========
Get number of trajectories to be analysed
========*/
int get_nbCoordFiles(char filPath[500]){
//Read each trajectory 
  DIR* rep = NULL;
  struct dirent* cordFile = NULL; 
  if((rep =opendir(filPath)) == NULL){ printf("Can\'not open directory %s\n", filPath); exit(5); };
  int num_traj =0;
	//Get the number of files in rep
  	while ((cordFile = readdir(rep)) != NULL){
		if(strrchr(cordFile->d_name, '.') != NULL){ //Check if it is a directory
			//Only ".xyz" and ".xmol" are taken into account, one can add more extension or remove this verification
			if(strcmp(strrchr(cordFile->d_name,'.'),".xyz")==0 || 
		   	   strcmp(strrchr(cordFile->d_name,'.'),".XYZ")==0 || 
		       strcmp(strrchr(cordFile->d_name,'.'),".xmol")==0 || 
		       strcmp(strrchr(cordFile->d_name,'.'),".XMOL")==0){
				num_traj++;
			}
		}
  	}
  //Close the directory
  closedir(rep);
  //Return the number of trajectories found in filPath
  return num_traj;
}

/*========
Get trajectories files names to be analysed
========*/
void get_FileList(struct MyModel *molecule){

  	DIR* rep = NULL;
  	struct dirent* cordFile = NULL; 
    //Allocate memory to save the trajectories files to be analysed
  	molecule->inputFNList= allocate_char_matrix(get_nbCoordFiles(molecule->inputDir),256,"coordFiles");
  	//Associate to each trajectory file a number molecule->num_F, this is used in another function
  	molecule->nbInputF=0;
  	molecule->num_F=0;
  	if((rep =opendir(molecule->inputDir)) == NULL){ printf("Can\'not open directory %s\n", molecule->inputDir); exit(5); };
	//Get list of files
  	while ((cordFile = readdir(rep)) != NULL){
		//Only ".xyz" and ".xmol" are taken into account, one can add more extension or remove this verification	
		if(strrchr(cordFile->d_name, '.') != NULL){ //Check if it is a directory
			if(strcmp(strrchr(cordFile->d_name,'.'),".xyz")==0 || 
		   	   strcmp(strrchr(cordFile->d_name,'.'),".XYZ")==0 || 
		       strcmp(strrchr(cordFile->d_name,'.'),".xmol")==0 || 
		       strcmp(strrchr(cordFile->d_name,'.'),".XMOL")==0){
				sprintf(molecule->inputFNList[molecule->nbInputF],"%s", cordFile->d_name);
				molecule->nbInputF++;
			}
		}
  	}
	//Close the directory
  	closedir(rep);
  	//Check the files list
  	#ifdef DISP_FILE_LIST
  		display_results(molecule);
  	#endif

}

/*========
Delete the file fileName if it exists
========*/
void del_file(char *fileName){
	FILE *f = fopen(fileName, "r") ;
	if(f){
		fclose(f) ;
		chmod(fileName, S_IWRITE) ;
		remove(fileName) ;		
	}
}    

/*========
Delete the directory dirName if it exists, with all its content
========*/
void del_dir(char *dirName){
	DIR *dir;
	struct dirent *ent;
	dir = opendir(dirName) ;
	if(dir){
		while((ent = readdir(dir)) != NULL){
			if ( strcmp(ent->d_name, ".") != 0 &&  strcmp(ent->d_name, "..") != 0){
				char fileName[512];
				sprintf(fileName,"%s/%s",dirName,ent->d_name);
				del_dir(fileName) ;
			}
		}		
		closedir(dir);
	}
	remove(dirName);
}

/*========
Size of the first snapshot
========*/
int snapshot_size(FILE* traj){
	int size;
	if(fscanf(traj,"%d", &size)!=1){
      	printf("Cannot read number of atoms in trajectory .xyz file\n");
      	exit(-1);
   	}
	//Set the file position to the beginning of the file
	rewind(traj); 
	//Return the size of one snapshot = number of atoms + 2 lines of comments
	return size+2;
}

/*========
Skip num_snap snapshots in trajectory  
========*/
// void skip_snapshot(FILE* traj, int num_snap, int size){
// 	int cpt_snap=0;
// 	while(cpt_snap<num_snap){
//   		skip_return(traj,size);
//   		cpt_snap++;
// 	}
// }

/*========
Viewing format (AtomNumber)
========*/
char* format_str(char c, int index){
	static char chaine[4]="";
	if(index<10) sprintf(chaine,"%c%d ",c,index); else  sprintf(chaine,"%c%d",c,index);
	return chaine;
}

/*========
Get a substring from string
========*/
char *str_sub (const char *s, unsigned int start, unsigned int end){
   char *new_s = NULL;
   	if (s != NULL && start < end){
      new_s = malloc (sizeof (*new_s) * (end - start + 2));
      if (new_s != NULL){
		unsigned int i;
		for (i = start; i <= end; i++)
		new_s[i-start] = s[i];
		new_s[i-start] = '\0';
      }
      else{
        fprintf (stderr, "Memory overflow\n");
        exit (-1);
      }
   	}
   return new_s;
}

/*========
Put char c to string word at position i
========*/ 
char* put_char(char *word, int index, char c){	
	static char tmp[256]="";
	int size=strlen(word);
	int i; for(i=0;i<size;i++) tmp[i]=word[i]; 
	tmp[index]=c;
	tmp[size]='\0';
	return tmp;
}

/*========
Make a copy of the  matrix1 into matrix2
========*/
void copy_matrix(int **matrix1,int **matrix2,int size){
	int i; int j;
	for (i = 0; i < size; i++) for (j = 0; j < size; j++) matrix2[i][j]=matrix1[i][j];
}      

/*========
Make a copy of the  vector tab1 into the vector tab2
========*/
void copy_tab(int *tab1,int *tab2,int size){
	int i; 
	for (i = 0; i < size; i++)  tab2[i]=tab1[i];
}      

/*========
Initialise the vector tab of integers with the value val
========*/
void init_tab(int *tab,int size, int val){
	int i; 
	for (i = 0; i < size; i++) tab[i]=val;	
}

/*========
Initialise the vector tab of doubles with the value val
========*/
void init_tab_double(double *tab,int size, double val){
	int i; 
	for (i = 0; i < size; i++) tab[i]=val;	
}

/*========
Initialise the matrix matrix of integers with the value val
========*/
void init_matrix(int **matrix,int sizeL, int sizeC, int val){
	int i; int j;
	for (i = 0; i < sizeL; i++) for (j = 0; j < sizeC; j++) matrix[i][j] = val;
}

/*========
Initialise the matrix matrix of doubles with the value val
========*/
void init_matrix_double(double **matrix,int sizeL, int sizeC, double val){
	int i; int j;
	for (i = 0; i < sizeL; i++) for (j = 0; j < sizeC; j++) matrix[i][j] = val;
}

/*========
Calculate the euclidean distance between two points (atoms)
========*/
double distance (double xi,double yi,double zi,double xj,double yj,double zj){
	return sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj));
}

/*========
Calculate the euclidean distance between two atoms considering the Periodic Boundary Conditions (PBC)
========*/
double distance_pbc (struct atom atomA, struct atom atomB, double boxX,double boxY,double boxZ){
	//Check X difference 
	double xdiff = atomA.x - atomB.x;
	if(xdiff - boxX/2 > EPSILON){ // we take the nearest neighbour 
		xdiff = boxX - xdiff ; 
	}
	else{
		if(xdiff + boxX/2 < EPSILON){
			xdiff = boxX + xdiff ;
		}
	}
	//Check Y difference 
	double ydiff = atomA.y - atomB.y;
	if(ydiff -boxY/2 > EPSILON){ // we take the nearest neighbour 
		ydiff = boxY - ydiff ; 
	}
	else{
		if(ydiff + boxY/2 < EPSILON){
			ydiff = boxY + ydiff ;
		}
	}
	//Check Z difference 
	double zdiff = atomA.z - atomB.z;
	if(zdiff-boxZ/2>EPSILON){ // we take the nearest neighbour 
		zdiff = boxZ - zdiff ; 
	}
	else{
		if(zdiff+ boxZ/2 < EPSILON){
			zdiff = boxZ + zdiff ;
		}
	}
	
	//Compute distance 	
	return sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff);
}

/*========
Calculate the angle between three atoms
========*/
double angle (double xi,double yi,double zi,double xj,double yj,double zj,double xk,double yk,double zk){
	// consider i as the donor atom (heavy atom) ; j the Hydrogen atom  and k the acceptor
	return acos(((xi-xj)*(xk-xj)+(yi-yj)*(yk-yj)+(zi-zj)*(zk-zj))/(distance (xi,yi,zi,xj,yj,zj)*distance (xk,yk,zk,xj,yj,zj)));
}

/*========
Calculate the angle between two atoms considering a periodic box
========*/
double angle_pbc(struct atom atomD, struct atom atomH,struct atom atomA, double boxX,double boxY,double boxZ){
	// consider atomD as the donor atom ; atomH the H and atomA the acceptor
	//AH
	//Check X difference 
	double xdiffAH = atomA.x - atomH.x;
	if(xdiffAH - boxX/2 > EPSILON){ // we take the nearest neighbour 
		xdiffAH = boxX - xdiffAH ; 
	}
	else{
		if(xdiffAH  + boxX/2 < EPSILON){
			xdiffAH = boxX + xdiffAH ;
		}
	}
	//Check Y difference 
	double ydiffAH = atomA.y - atomH.y;
	if(ydiffAH-boxY/2>EPSILON){ // we take the nearest neighbour 
		ydiffAH = boxY - ydiffAH ; 
	}
	else{
		if(ydiffAH + boxY/2 < EPSILON){
			ydiffAH = boxY + ydiffAH ;
		}
	}
	//Check Z difference 		
	double zdiffAH = atomA.z - atomH.z;
	if(zdiffAH-boxZ/2>EPSILON){ // we take the nearest neighbour 
		zdiffAH = boxZ - zdiffAH ; 
	}
	else{
		if(zdiffAH + boxZ/2 <EPSILON){
			zdiffAH = boxZ + zdiffAH ;
		}
	}

	//DH
	//Check X difference 
	double xdiffDH = atomD.x - atomH.x;
	if(xdiffDH - boxX/2>EPSILON){ // we take the nearest neighbour 
		xdiffDH = boxX - xdiffDH ; 
	}
	else{
		if(xdiffDH + boxX/2 < EPSILON){
			xdiffDH = boxX + xdiffDH ;
		}
	}
	//Check Y difference 
	double ydiffDH = atomD.y - atomH.y;
	if(ydiffDH-boxY/2>EPSILON){ // we take the nearest neighbour 
		ydiffDH = boxY - ydiffDH ; 
	}
	else{
		if(ydiffDH + boxY/2 <EPSILON){
			ydiffDH = boxY + ydiffDH ;
		}
	}		
	//Check Z difference 	
	double zdiffDH = atomD.z - atomH.z;
	if(zdiffDH-boxZ/2>EPSILON){ // we take the nearest neighbour 
		zdiffDH = boxZ - zdiffDH ; 
	}
	else{
		if(zdiffDH + boxZ/2 < EPSILON){
			zdiffDH = boxZ + zdiffDH ;
		}
	}
	double ang = acos((xdiffAH*xdiffDH+ydiffAH*ydiffDH+zdiffAH*zdiffDH)/(sqrt(xdiffAH*xdiffAH+ydiffAH*ydiffAH+zdiffAH*zdiffAH)*sqrt(xdiffDH*xdiffDH+ydiffDH*ydiffDH+zdiffDH*zdiffDH)));

	//((xdiffAH*xdiffDH+ydiffAH*ydiffDH+zdiffAH*zdiffDH)/(sqrt(xdiffAH*xdiffAH+ydiffAH*ydiffAH+zdiffAH*zdiffAH)*sqrt(xdiffDH*xdiffDH+ydiffDH*ydiffDH+zdiffDH*zdiffDH)));
	return ang;	
}

/*========
Get the number of connected components using the Breadth-First Search algorithm 
========*/
int get_cc(int** CC, int size, bool save, int **frg){
	int cp = 0;
	int i=0;
	struct path* F = (struct path*)malloc(sizeof(struct path));   //queue : FIFO
	//Check that the request space if allocated without problem
	#ifdef VERIF_ACCES
		if(!F){fprintf(stderr,"Cannot malloc filed for FIFO \n"); exit(1);}
	#endif
	
	F->head = NULL;
	F->queue = NULL;
	// int p[size]; // p[i] = j such as the j lead to put i in the queue
	int color[size]; // color[i]= -1 is vertex i not visited , 0 if visited , 1 there is path from S0 to i
	//Initilization
	for(i=0;i<size;i++){color[i]=-1; /*p[i]=-1 ;*/ /*d[i]=0;*/} 

	//Get the #connected components using the Breadth-first search algorithm 
	for(i=0;i<size;i++){
		//if(CC[i][i]!=0) color[i]=0; //Ignore argon / cation / anion
		if(color[i]==-1){ //Still we have a vertex not marked 
			//Start the BFS
			int d=0; //size of the fragement 
			int x = i; //Start point : i
			int y = -1;
			//d[x]=0 ;  //This is used to  get the distance of a path
			color[x]=0;
			F = add_elt(F,x);
			//Start the BFS
			while(F->head!=NULL){
				x  = get_elt_FIFO(F);
				if(save){frg[cp][x]=1;} //save the element of the fragment
				d++;
				bool end= false;
				while(!(end)){
					y = get_elt_notMarked(CC, size, x , color);
					if(y != -1){
						if(color[y]==-1){
							F =add_elt(F,y);
							color[y]=0;
							// p[y]=x;
							//d[y]=d[y]+1 ;
						}	
						// if(CC[x][y]==2 || CC[y][x]==2){ //It's an hydrogen bond
						// 	*impDA = 1; 		
						// }
					}
					else
						end = true;
				}
			 	color[x]=1;
			}
			if(save) frg[cp][size]=d; // save the size of the fragment
			if(d>1) {
				cp++;
				if(save){
					frg[cp][size]=d;	
				}
			}
			else{
				if(save){
					frg[cp][x]=-1;
				}	
			}
		}
	}
	free(F);
	return cp;	
}

/*========
Search the diffrent directed paths using the Breadth-First Search algorithm 
========*/
int get_dir_path(int** CC, int size, int **motifs){
	int cp = 0;
	int i=0;
	struct path* F = (struct path*)malloc(sizeof(struct path));   //queue : FIFO
	//Check that the request space if allocated without problem
	#ifdef VERIF_ACCES
		if(!F){fprintf(stderr,"Cannot malloc filed for FIFO \n"); exit(1);}
	#endif
	
	F->head = NULL;
	F->queue = NULL;
	// int p[size]; // p[i] = j such as the j lead to put i in the queue
	int color[size]; // color[i]= -1 is vertex i not visited , 0 if visited , 1 there is path from S0 to i
	//Initilization
	for(i=0;i<size;i++){color[i]=-1; /*p[i]=-1 ;*/ /*d[i]=0;*/} 

	//Get the #connected components using the Breadth-first search algorithm 
	for(i=0;i<size;i++){
		//if(CC[i][i]!=0) color[i]=0; //Ignore argon / cation / anion
		if(color[i]==-1){ //Still we have a vertex not marked 
			//Start the BFS
			int d=0; //size of the fragement 
			int x = i; //Start point : i
			int y = -1;
			//d[x]=0 ;  //This is used to  get the distance of a path
			color[x]=0;
			F = add_elt(F,x);
			//Start the BFS
			while(F->head!=NULL){
				x  = get_elt_FIFO(F);
				motifs[cp][x]=1;//save the element of the fragment
				d++;
				bool end= false;
				while(!(end)){
					y = get_elt_notMarked2(CC, size, x , color);
					if(y != -1){
						if(color[y]==-1){
							F =add_elt(F,y);
							color[y]=0;
							// p[y]=x;
							//d[y]=d[y]+1 ;
						}	
						// if(CC[x][y]==2 || CC[y][x]==2){ //It's an hydrogen bond
						// 	*impDA = 1; 		
						// }
					}
					else
						end = true;
				}
			 	color[x]=1;
			}
			motifs[cp][size]=d; // save the size of the fragment
			if(d>1) {
				cp++;
				motifs[cp][size]=d;	
			}
			else{
				motifs[cp][x]=-1;
			}
		}
	}
	free(F);
	printf("cp=%d\n",cp );
	return cp;	
}

/*========
Initialize the adjacency matrix with the different bonds
========*/
void init_adj_conf_matrix(int** CC, int nbAtom,struct confList* conf, int level, char type){
	int i;
	int j;
	//if((bit_1(level,0) && type=='i')|| type=='f'){ //check 1/0 with -1 , i for isomorphism and f for fragment 
		for (i = 0; i < nbAtom; i++){
			for (j = i+1; j < nbAtom; j++){
				if(conf->img[i].atomIn && conf->img[j].atomIn){
					if((conf->HB[i][j]==-5 || (conf->CB[i][j]>0 && conf->img[i].atomType!='H' && conf->img[i].atomType!='H'))){
						if(conf->HB[i][j]==-5){
							if(type=='p'){//search directed paths
								int D=-1;
								check_ptr(conf,i,j,conf->nbAtom,&D);							
								if(D==i){
									CC[i][j]=1;
								}
								else{
									CC[j][i]=1;
								}
							}	
							else{
								CC[i][j]=1;
								CC[j][i]=1;							
							}
						}
						else{
							CC[i][j]=1;
							CC[j][i]=1;
						}
					}
				}
			}
		}
	//}
	if((bit_1(level,2) /*&& type=='i'*/)/*|| type=='f')*/){ //check 1/0 with -1
		for(i=0;i<conf->nbIon;i++){
			if(conf->IB[i][nbAtom]!=-1 && conf->IB[i][nbAtom+1]>0 ){
				for(j=0;j<nbAtom;j++){
					if(conf->IB[i][j]>0 && !((strcmp(conf->img[j].atomName,"Ar")==0  || strcmp(conf->img[j].atomName,"Cl")==0  || strcmp(conf->img[j].atomName,"Na")==0 || strcmp(conf->img[j].atomName,"Li")==0 || strcmp(conf->img[j].atomName,"K")==0 ||strcmp(conf->img[j].atomName,"I")==0 || strcmp(conf->img[j].atomName,"F")==0 || strcmp(conf->img[j].atomName,"Br")==0) && get_index_ion(conf->IB,conf->nbIon,conf->nbAtom,j) < i) && conf->img[conf->IB[i][conf->nbAtom]].atomIn && conf->img[j].atomIn ){ //Intermolecular bond
						CC[conf->IB[i][conf->nbAtom]][j]=1;
						CC[j][conf->IB[i][conf->nbAtom]]=1;
					}		
				}
			}
		}
	}
	if((bit_1(level,3)/* && type=='i'*/)/*|| type=='f'*/){ //check 1/0 with -1
		for(i=0;i<conf->nbMetal;i++){
			if(conf->MB[i][nbAtom]!=-1 && conf->MB[i][nbAtom+1]>0 ){
				for(j=0;j<nbAtom;j++){
					if(conf->MB[i][j]>0 && !((strcmp(conf->img[j].atomName,"Mn")==0 ||strcmp(conf->img[j].atomName,"Ru")==0 || strcmp(conf->img[j].atomName,"Au")==0 ) && get_index_ion(conf->MB,conf->nbMetal,conf->nbAtom,j) < i) && conf->img[conf->MB[i][conf->nbAtom]].atomIn && conf->img[j].atomIn){ //Intermolecular bond
						CC[conf->MB[i][conf->nbAtom]][j]=1;
						CC[j][conf->MB[i][conf->nbAtom]]=1;
					}		
				}
			}
		}
	}
}



/*========
Get the number of the connected components (cc) and save the distribution. This function is used in the interface analysis module.
Some treatments are performed here about the structure of each connected component, mainly, size of the connected component. 
========*/
int get_cc1(int** CC, int size, char outputFN[], int *atomIgnored, struct CCModel* ccDistrib/*,int nbHAtom*/){
	int cp = 0;
	int i=0;
	init_tab(ccDistrib->bigCCtab,size,-1);
	ccDistrib->smallCC=size;
	struct path* F = (struct path*)malloc(sizeof(struct path));   //queue : FIFO
	#ifdef VERIF_ACCES
		if(!F){fprintf(stderr,"Cannot malloc filed for FIFO \n"); exit(1);}
	#endif
	
	F->head = NULL;
	F->queue = NULL;
	//File to save CCs sizes and draw them with gnuplot
	char inputFN[500]="cc.txt";
    	FILE* inputF=fopen(inputFN,"w");	
	if(inputF == NULL){ printf("Can\'not open file  cc.txt\n"); exit(4);}  

	//The output file name 
	char outputG[500];
	sprintf(outputG,"%s_cc.png",outputFN) ;

	//int p[size]; // p[i] = j such as the j lead to put i in the queue
	int color[size]; // color[i]= -1 is vertex i not visited , 0 if visited , 1 there is path from S0 to i
	//Initilization
	for(i=0;i<size;i++){color[i]=-1;} //p[i]=-1 ; //d[i]=0;

	//Get the #connected components using the Breadth-first search algorithm 
	int bigCC=0;
	for(i=0;i<size;i++){
		//if(CC[i][i]!=0) color[i]=0; //Ignore argon / cation / anion
		if(color[i]==-1){ //Still we have a vertex not marked 
			int ccsize=0;
			int bigCCtab[size];
			init_tab(bigCCtab,size,-1);
			//Start the BFS
			int x = i; //Start point : i
			int y = -1;
			//d[x]=0 ; 
			color[x]=0;
			F = add_elt(F,x);
			//Start the BFS
			while(F->head!=NULL){
				x  = get_elt_FIFO(F);
				bigCCtab[ccsize]=x;
				ccsize++;
				bool end= false;
				while(!(end)){
					y = get_elt_notMarked(CC, size, x , color);
					if(y != -1){
						F =add_elt(F,y);
						color[y]=0;
					}
					else
						end = true;
				}
			 	color[x]=1;
			}	
			//save only CC with size > CC_SIZ_MIN 
			if(ccsize > CC_SIZ_MIN){
				// printf("new cc %d , size %d , from i=%d\n", cp+1,ccsize,i );
				cp++;
				//Save the name of the connected component (number) and its size
				fprintf(inputF, "%d \t %d \t \"%d\" \n",cp, ccsize ,  cp);
				//Increment the number of the connected component with size "ccsize/3"
				ccDistrib->tabCC[ccsize]++; 
				//Get the size of the biggest connected component
				ccDistrib->bigCC= Max(ccDistrib->bigCC,ccsize); 
				ccDistrib->smallCC = Min(ccDistrib->smallCC,ccsize);	
				if(bigCC<ccsize){
					bigCC= Max(bigCC,ccsize);					
					//Update the copy bigCCtab by setting the new biggest connected component
					copy_tab(bigCCtab, ccDistrib->bigCCtab,ccsize);
				}	
			}
			else{ //compute water molecules ignored  

				*atomIgnored +=ccsize;
			}
			ccsize=0;
		}
	}
	free(F);
	fclose(inputF);

	//Save the local biggest connected component 
	ccDistrib->averageSize+=bigCC;
	ccDistrib->tabCC[0]=bigCC;
	// printf("tabCC[0]=%d , bigcc=%d\n",ccDistrib->tabCC[0],bigCC);
	//Draw the number of the connected components according to their size using gnuplot
	draw_interface_distribution(inputFN, outputG, (int)(size)+1, "Distribtion of water molecule according CC.","C.C. label","C.C. size (%)",1) ;
	//Return the number of connected components in the graph
	return cp;	
}

/*========
Get the list of atoms involved in the conformational search
========*/
void get_atom_in_conf(struct MyModel *molecule,struct atom* img){
	int i=0;
	int j=0;
	int** CC = allocate_matrix(molecule->nbAtom,molecule->nbAtom,0,"CC"); //Adjacency matrix
	
	//Covalent & hydrogen bonds
	for (i = 0; i < molecule->nbAtom; i++){
		CC[i][i]=0; //All nodes are not marked
		for (j = i+1; j < molecule->nbAtom; j++){
			if(molecule->covBond[i][j]>0 ){ 
				CC[i][j]=1;
				CC[j][i]=1;
			}
			if(molecule->donHacc[i][j]==-5 && img[i].atomType!='H' && img[j].atomType!='H'){
					CC[i][j]=1;
					CC[j][i]=1;
			}

		}
	} 

	//Intermolecular interactions
	for(i=0;i<molecule->nbIon;i++){
		if(molecule->ionBond[i][molecule->nbAtom]!=-1 && molecule->ionBond[i][molecule->nbAtom+1]>0 ){
			for(j=0;j<molecule->nbAtom;j++){
				if(molecule->ionBond[i][j]>0 && !((strcmp(molecule->img[j].atomName,"Ar")==0 || strcmp(molecule->img[j].atomName,"Li")==0  ||strcmp(molecule->img[i].atomName,"Na")==0 || strcmp(molecule->img[j].atomName,"K")==0 || strcmp(molecule->img[j].atomName,"F")==0) && get_index_ion(molecule->ionBond,molecule->nbIon,molecule->nbAtom,j) < i)){ //Intermolecular bond
					CC[molecule->ionBond[i][molecule->nbAtom]][j]=1;
					CC[j][molecule->ionBond[i][molecule->nbAtom]]=1;
				}		
			}
		}
	}

	//Organometallic interactions
	for(i=0;i<molecule->nbMetal;i++){
		if(molecule->metalBond[i][molecule->nbAtom]!=-1 && molecule->metalBond[i][molecule->nbAtom+1]>0 ){
			for(j=0;j<molecule->nbAtom;j++){
				if(molecule->metalBond[i][j]>0 && !((strcmp(molecule->img[j].atomName,"Mn")==0 || strcmp(molecule->img[j].atomName,"Ru")==0 || strcmp(molecule->img[j].atomName,"Au")==0) && get_index_ion(molecule->metalBond,molecule->nbMetal,molecule->nbAtom,j) < i)){ //Intermolecular bond
					CC[molecule->metalBond[i][molecule->nbAtom]][j]=1;
					CC[j][molecule->metalBond[i][molecule->nbAtom]]=1;
				}		
			}
		}
	}
	//Get the connected components and update atomIn (img)
	update_atom_in(CC,img,molecule->nbAtom);

	free_matrix(CC,molecule->nbAtom);
}

/*========
Get the atoms in the interface layer
========*/
void get_atom_in_interface(struct confList *conf){
	int i=0;
	int j=0;
	// struct confList *conf = molecule->conformations;
	double cmz = center_mass_axis(conf->img,conf->nbAtom,'z');
	//check water molecules
	for(i=0;i<conf->nbAtom;i++){
		conf->img[i].atomIn=false;
		if(strcmp(conf->img[i].atomName,"O")==0 && get_nb_CB_type(conf->img,i,conf->CB,conf->nbAtom,"H")==2){
			//Check that Z>0 (up layer) 
			if(conf->img[i].z >=cmz){
				conf->img[i].atomIn=true;	
			}
		}
	}	
	// //Check all the oxygen that are involved in hydrogen bond
	// for(i=0;i<conf->nbAtom;i++){
	// 	if(strcmp(conf->img[i].atomName,"O")==0 && check_inv_HB(conf->HB,conf->nbAtom,i)){
	// 		if(conf->img[i].z >=cmz){
	// 			conf->img[i].atomIn=true;	
	// 		}
	// 	}
	// }

	// Check the Oxygen attached to only 1 Si => Oint and are involved in hydrogen bond  
	for(i=0;i<conf->nbAtom;i++){
		if(strcmp(conf->img[i].atomName,"O")==0 && (get_nb_CB_type(conf->img,i,conf->CB,conf->nbAtom,"Si")==1 /*|| get_nb_CB_type(conf->img,i,conf->CB,conf->nbAtom,"Si")==2*/) && check_inv_HB(conf->HB,conf->nbAtom,i)){
			if(conf->img[i].z >=cmz){
				conf->img[i].atomIn=true;	
			}
		}
	}
	// Check the Si and H atoms that are attached to Oxygen Oint  
	for(i=0;i<conf->nbAtom;i++){
		if(strcmp(conf->img[i].atomName,"Si")==0 || strcmp(conf->img[i].atomName,"H")==0){
			// printf("%s\n",conf->img[i].atomName,"Si");
			for(j=0;j<conf->nbAtom;j++){
				if(strcmp(conf->img[j].atomName,"O")==0 && conf->img[j].atomIn){
					if(conf->CB[i][j]==1 || conf->CB[j][i]==1){
						conf->img[i].atomIn = true;
						// printf("IN\n");
					}
				}
			}

		}
		if( (strcmp(conf->img[i].atomName,"Cl")==0 ||  strcmp(conf->img[i].atomName,"K")==0) &&  conf->img[i].z >=cmz){
			conf->img[i].atomIn = true;
		}

	}	

}
/*========
Get The statistics on the conformers in order to create a graph of transitions
========*/
void get_stat(struct MyModel *molecule){
	int i=0;
	int j=0;
	struct confList *conf = molecule->conformations;
	int cptHBWW=0;// #HB water-water
	int cptHBWS=0; // #HB water-solid
	int cptHBSS=0; // #HB solid-solid

	//Get the number of HB per type
	FILE *outputF=fopen(get_FN(molecule->inputDir,molecule->inputFN,"_stat.txt"),"w");
	while(conf){
		for (i = 0; i < conf->nbAtom; i++){
			for(j=i+1;j<conf->nbAtom;j++){
				if(conf->img[i].atomIn && conf->img[j].atomIn){
					if(strcmp(conf->img[i].atomName,"O")==0 && strcmp(conf->img[j].atomName,"O")==0 && (conf->HB[i][j]==-5 || conf->HB[j][i]==-5)){
						if(get_nb_CB_type(conf->img,i,conf->CB,conf->nbAtom,"Si")==1){ //solid 
							if(get_nb_CB_type(conf->img,j,conf->CB,conf->nbAtom,"Si")==1){ //solid-solid HB
								cptHBSS++;
							}
							else{//water-solid HB
								cptHBWS++;
							}
						}
						else{ //water-solid
							if(get_nb_CB_type(conf->img,j,conf->CB,conf->nbAtom,"Si")==1){ //solid-water HB
								cptHBWS++;
							}
							else{//water-water HB
								cptHBWW++;
							}
						}
					}	
				}
			}
		}
		//Get the number of connected component
		int **CC=allocate_matrix(conf->nbAtom,conf->nbAtom,0,"CC-interf");
		init_adj_conf_matrix(CC, conf->nbAtom, conf, molecule->level,'f');
		int cptCC = get_cc(CC,conf->nbAtom,false,NULL);
		//Get the number of cycles 
		//Get #atomsIn
		int cptAtomIN = get_nb_atomIn(conf->img,conf->nbAtom); 
		//Get #bonds 
		int cptBonds = get_nb_bonds(CC,conf->img,conf->nbAtom);
		//Get the number of rings
		int cptRing = cptBonds - cptAtomIN + cptCC ; 
		fprintf(outputF,"conf%d  \t%d \t%d \t%d \t%d \t%d \t%d \t%d \n",conf->name, cptHBWW,cptHBWS,cptHBSS,cptCC, cptBonds ,cptAtomIN, cptRing); 
		free_matrix(CC,conf->nbAtom);
		conf=conf->suiv;
	}	
	fclose(outputF);
}
/*========
Get the number of atomIn that are activated
========*/
int get_nb_atomIn(struct atom *img, int size){
	int cpt=0;
	for (int i = 0; i < size; i++){
		if(img[i].atomIn && (strcmp(img[i].atomName,"O")==0 || strcmp(img[i].atomName,"Si")==0 || strcmp(img[i].atomName,"K")==0 || strcmp(img[i].atomName,"Cl")==0 ) ) cpt++;
	}
	return cpt;
}

/*========
Get the number of bonds (only atomIn)
========*/
int get_nb_bonds(int **CC, struct atom *img, int size){
	int cpt=0;
	for (int i = 0; i < size; i++){
		for (int j = i+1; j < size; j++){
			if(img[i].atomIn && img[j].atomIn && CC[i][j]>0){
				if((strcmp(img[i].atomName,"O")==0 || strcmp(img[i].atomName,"Si")==0 || strcmp(img[i].atomName,"K")==0 || strcmp(img[i].atomName,"Cl")==0 ) && (strcmp(img[j].atomName,"O")==0 || strcmp(img[j].atomName,"Si")==0 || strcmp(img[j].atomName,"K")==0 || strcmp(img[j].atomName,"Cl")==0 ))
					cpt++;
			}
		}
	}
	return cpt;
}

/*========
Get the number of covalent bonds of atom i with atom of type "type"
========*/
int get_nb_CB_type(struct atom *img,int index ,int **CB,int size,char *type){
	int j=0;
	int cpt = 0; 
	//Research by lines
	while (j<index){
		if(CB[j][index]>0 && strcmp(img[j].atomName,type)==0){
			cpt++;
			j++;
		}
		else
			j++;
	}
	if(j==index) j++;
	//Research by columns
	while (j<size){
		if(CB[index][j]>0 && strcmp(img[j].atomName,type)==0){
			cpt++;
			j++;
		}
		else
			j++;
	} 
	return cpt;

}

/*========
Check if atom i is involved in a hydrogen bond 
========*/
bool check_inv_HB(int **HB, int size, int index){
	int j=0;
	//Research by lines
	while (j<index){
		if(HB[j][index]==-5)
			return true;
		else
			j++;
	}
	if(j==index) j++;
	//Research by columns
	while (j<size){
		if(HB[index][j]==-5)
			return true;
		else
			j++;
	} 
	return false;
}

/*========
Center the atoms in the box 
========*/
void center_atoms(char *outputFN,struct atom* img, int size, double boxX, double boxY, double boxZ ){
	int i;
	FILE* outputF=fopen(outputFN,"w");	
	if(outputF == NULL){ printf("Can\'not open file  %s\n",outputFN); exit(4);} 
	int cptatom=0;
	for(i=0;i<size;i++) if(img[i].atomIn) cptatom++; 

	fprintf(outputF,"%d\n",cptatom);
	fprintf(outputF,"center box\n");
	for(i=0;i<size;i++){
		if(img[i].atomIn){
			fprintf(outputF,"%s\t",img[i].atomName);
			//Check the X value 
			if(img[i].x < -boxX/2.00){
				fprintf(outputF,"%lf \t",img[i].x+boxX);
			}
			else{
				if(img[i].x > boxX/2.00){
					fprintf(outputF,"%lf \t",img[i].x-boxX);
				}
				else{
					fprintf(outputF,"%lf \t",img[i].x);
				}
			}

			//Check the Y value 
			if(img[i].y < -boxY/2.00){
				fprintf(outputF,"%lf \t",img[i].y+boxY);
			}
			else{
				if(img[i].y > boxY/2.00){
					fprintf(outputF,"%lf \t",img[i].y-boxY);
				}
				else{
					fprintf(outputF,"%lf \t",img[i].y);
				}
			}

			//Check the Z value 
			if(img[i].z < -boxZ/2.00){
				fprintf(outputF,"%lf \n",img[i].z+boxZ);
			}
			else{
				if(img[i].z > boxZ/2.00){
					fprintf(outputF,"%lf \n",img[i].z-boxZ);
				}
				else{
					fprintf(outputF,"%lf \n",img[i].z);
				}
			}
		}
	}
	fclose(outputF); 

}
/*========
Center the atoms in the box 
========*/
double center_mass_axis(struct atom* img, int size, char axis){
	int i;
	double cm=0.00; 
	for(i=0;i<size;i++){
		switch (axis){
			case 'x':
				cm +=img[i].x;	
			break;
			case 'y':
				cm +=img[i].y;	
			break;
			case 'z':
				cm +=img[i].z;	
			break;

			default:
				printf("ERROR in axis \n");
			break;
		}
	}
	return cm/(double)(size);
}

/*========
Browse the fragments based on CC adjacency matrix and update the atomIn field in img structure 
========*/
void update_atom_in(int **CC, struct atom* img, int size){
	int i=0;
	struct path* F = (struct path*)malloc(sizeof(struct path));   //queue : FIFO
	//Check that the request space if allocated without problem
	#ifdef VERIF_ACCES
		if(!F){fprintf(stderr,"Cannot malloc filed for FIFO \n"); exit(1);}
	#endif
	
	F->head = NULL;
	F->queue = NULL;
	//int p[size]; // p[i] = j such as the j lead to put i in the queue
	int color[size]; // color[i]= -1 is vertex i not visited , 0 if visited , 1 there is path from S0 to i
	//Initilization
	for(i=0;i<size;i++){color[i]=-1;} //p[i]=-1 ; //d[i]=0;

	//Get the #connected components using the Breadth-first search algorithm 
	for(i=0;i<size;i++){
		//We choose S0 as cation, anion or metal atom // Try to automatize this part 
		// if(	strcmp(img[i].atomName,"N")==0){ //|| 
		if(	strcmp(img[i].atomName,"Mn")==0 || strcmp(img[i].atomName,"F")==0 || strcmp(img[i].atomName,"Ru")==0 ){
			// strcmp(img[i].atomName,"Li")==0	||
			// strcmp(img[i].atomName,"Ar")==0	||
			// strcmp(img[i].atomName,"K")==0	){
			if(color[i]==-1){ //Still we have a vertex not marked 
				//Start the BFS
				int x = i; //Start point : i
				int y = -1;
				color[x]=0;
				F = add_elt(F,x);
				//Start the BFS
				while(F->head!=NULL){
					x  = get_elt_FIFO(F);
					bool end= false;
					while(!(end)){
						y = get_elt_notMarked(CC, size, x , color);
						if(y != -1){
							F =add_elt(F,y);
							color[y]=0;
						}
						else
							end = true;
					}
					img[x].atomIn = true; // The atom belongs to a fragment containing an anion/cation or metal atom 
				 	color[x]=1;
				}			
			}
		}	
	}
	free(F);
}

/*========
Add element to de queue F
========*/
struct path* add_elt(struct path* F , int x){
	struct path_elt* newElt = (struct path_elt*)malloc(sizeof(struct path_elt));
	#ifdef VERIF_ACCES
		if(!newE){fprintf(stderr,"Cannot malloc filed for new element on FIFO %d\n",head->name); exit(1);}
	#endif
	newElt->index=x; 
	newElt->suiv = NULL;
	if(F->queue !=NULL)
		F->queue->suiv=newElt;
	else
		F->head = newElt;
	F->queue = newElt;
	return F;
}

/*========
Get the first element of the queue
========*/
// int get_elt(struct path* F){
// 	struct path_elt* elt = F->head;
// 	int x = elt->index ;
// 	F->head = F->head->suiv ; 
// 	if(F->head ==NULL) F->queue = NULL;
// 	free(elt);
// 	return x;
// }

/*========
Get a covalent bond not marked
========*/
int get_elt_notMarked(int** CC, int size, int index , int color[]){
	int j = 0; 
	//Research by lines
	while (j<index){
		if(CC[j][index]>0 && color[j]==-1)//bond not marked
			return j;
		else
			j++;
	}
	if(j==index) j++;
	//Research by columns
	while (j<size){
		if(CC[index][j]>0 && color[j]==-1)//bond not marked
			return j;
		else
			j++;
	} 
	return -1;
}

/*========
Get a bond not marked but in directed way
========*/
int get_elt_notMarked2(int** CC, int size, int index , int color[]){
	int j = 0; 
	while(j<size){
		if(CC[index][j]>0 && color[j]==-1 && j!=index)//bond not marked
			return j;
		else
			j++;

	}
	return -1;
}

/*========
Calculate the order of atom in the atomList, according to the type.
========*/
int num_occ(struct atom* img, int atomA, int nbAtom){
	int i= 0 ; int k=1;
	while(i!=atomA && i<nbAtom){
		if (strcmp(img[i].atomName,img[atomA].atomName)==0) k++ ; 
			i++;
		}
	if(i<nbAtom) return k ; else return -1;
}

/*========
Get the number of fragments (connected components) in terms of covalent
========*/
int get_nb_frg(int **covBond,int size){
	int i=0;
	int j=0;
	//int CC[size][size];
	int** CC = allocate_matrix(size,size,0,"CC"); //Adjacency matrix
	for (i = 0; i < size; i++){
		CC[i][i]=covBond[i][i]; //All nodes are not marked
		for (j = i+1; j < size; j++){
			if(covBond[i][j]>0){
				CC[i][j]=1;
				CC[j][i]=1;
			}
		}
	} 
	int nbFrg = get_cc(CC,size,false,NULL);
	free_matrix(CC,size);
	return nbFrg;
}

/*========
Get the number of fragments in conformation conf and save these fragments (connected components) in ccDistrib
========*/
int get_interface_frg(struct confList* conf,int size, char outputFN[],int *atomIgnored,struct CCModel* ccDistrib){
	int i=0;
	int j=0;
	//int CC[size][size];
	int** CC = allocate_matrix(size,size,0,"CC"); //Adjacency matrix
	for (i = 0; i < size; i++){
		CC[i][i]=0; //All nodes are not marked, hydrogen atom are ignored 
		for (j = i+1; j < size; j++){
			if(conf->CB[i][j]>0 && conf->img[i].atomType!='H' && conf->img[j].atomType!='H'){ 
				CC[i][j]=1;
				CC[j][i]=1;
			}
		}
	} 
	for(i=0;i<conf->nbAtom;i++){
		for(j=i+1;j<conf->nbAtom;j++){ //Update
			if(conf->HB[i][j]==-5 && conf->img[i].atomType!='H' && conf->img[j].atomType!='H'){
					CC[i][j]=1;
					CC[j][i]=1;
			}
		}					
	}
	for(i=0;i<conf->nbIon;i++){
		if(conf->IB[i][conf->nbAtom]!=-1 && conf->IB[i][conf->nbAtom+1]>0 ){
			for(j=0;j<conf->nbAtom;j++){
				if(conf->IB[i][j]>0 && !((strcmp(conf->img[j].atomName,"Ar")==0  || strcmp(conf->img[j].atomName,"Cl")==0  || strcmp(conf->img[j].atomName,"Na")==0 || strcmp(conf->img[j].atomName,"Li")==0 || strcmp(conf->img[j].atomName,"K")==0 ||strcmp(conf->img[j].atomName,"I")==0 || strcmp(conf->img[j].atomName,"F")==0 || strcmp(conf->img[j].atomName,"Br")==0) && get_index_ion(conf->IB,conf->nbIon,conf->nbAtom,j) < i) && conf->img[conf->IB[i][conf->nbAtom]].atomIn && conf->img[j].atomIn ){ //Intermolecular bond
					CC[conf->IB[i][conf->nbAtom]][j]=1;
					CC[j][conf->IB[i][conf->nbAtom]]=1;
				}		
			}
		}
	}

	int nbFrgHB = get_cc1(CC,size,outputFN,atomIgnored,ccDistrib/*,get_nbTypeatoms(conf->img,"H",conf->nbAtom)*/);
	free_matrix(CC,size);
	return nbFrgHB;
}

/*========
Get the number of atoms of type atomName within the molecular system
========*/
int get_nbTypeatoms(struct atom* img, char *atomName,int size){
	int i =0; int j=0 ;
	while(i<size) {if(strcmp(img[i].atomName,atomName)==0) j++; i++; }
	return j ;
}

/*========
Get index of occurrence occ of an atom with type atomName
========*/
int get_index(struct MyModel *molecule,char *atomName, int occ){
	int i = 0 ;
	while(i<molecule->nbAtom)
		if (strcmp(atomName,molecule->imgRef[i].atomName)==0 && num_occ(molecule->imgRef, i, molecule->nbAtom)==occ) return i; else i++;
	return  -1 ;
}

/*========
Get the index of an ion in the IB matrix
========*/
int get_index_ion(int **IB,int size1, int size2,int index){
	int i=0;
	while(i < size1){
		if(IB[i][size2]==index) return i; else i++;
	}
	return -1;
}

/*========
Check if atoms with indexes index1 and index2 can forms a hydrogen bond
========*/
bool pair_valid(struct MyModel* molecule, int index1, int index2){
	if(index1==index2)
		return false;
	if(strcmp(molecule->img[index1].atomName,"H")==0 || strcmp(molecule->img[index2].atomName,"H")==0)//one of atoms is an Hydrogen
		return false;
	//A carbon atom doesn't form an H-bond 
	if(strcmp(molecule->img[index1].atomName,"C")==0 || strcmp(molecule->img[index2].atomName,"C")==0 )//acceptor atoms is a carbon  //strcmp(molecule->img[index1].atomName,"C")==0 || 
		return false;
	//Cl
	if(strcmp(molecule->img[index1].atomName,"Cl")==0 || strcmp(molecule->img[index2].atomName,"Cl")==0)//one of the atoms is a carbon
		return false;
	//Li   
	if(strcmp(molecule->img[index1].atomName,"Li")==0 || strcmp(molecule->img[index2].atomName,"Li")==0)
		return false;
	//Ar
	if(strcmp(molecule->img[index1].atomName,"Ar")==0 || strcmp(molecule->img[index2].atomName,"Ar")==0)
		return false;
	//S
	if(strcmp(molecule->img[index1].atomName,"S")==0 || strcmp(molecule->img[index2].atomName,"S")==0)
		return false;
	//Mn
	if(strcmp(molecule->img[index1].atomName,"Mn")==0 || strcmp(molecule->img[index2].atomName,"Mn")==0)
		return false;
	//Ru
	if(strcmp(molecule->img[index1].atomName,"Ru")==0 || strcmp(molecule->img[index2].atomName,"Ru")==0)
		return false;
	//B
	if(strcmp(molecule->img[index1].atomName,"B")==0 || strcmp(molecule->img[index2].atomName,"B")==0)
		return false;
	//K
	if(strcmp(molecule->img[index1].atomName,"K")==0 || strcmp(molecule->img[index2].atomName,"K")==0)
		return false;
	//F
	if(strcmp(molecule->img[index1].atomName,"F")==0 || strcmp(molecule->img[index2].atomName,"F")==0)
		return false;
	//P
	if(strcmp(molecule->img[index1].atomName,"P")==0 || strcmp(molecule->img[index2].atomName,"P")==0)
		return false;
	//Na
	if(strcmp(molecule->img[index1].atomName,"Na")==0 || strcmp(molecule->img[index2].atomName,"Na")==0)
		return false;
	//I
	if(strcmp(molecule->img[index1].atomName,"I")==0 || strcmp(molecule->img[index2].atomName,"I")==0)
		return false;
	//Br
	if(strcmp(molecule->img[index1].atomName,"Br")==0 || strcmp(molecule->img[index2].atomName,"Br")==0)
		return false;
	//Au
	if(strcmp(molecule->img[index1].atomName,"Au")==0 || strcmp(molecule->img[index2].atomName,"Au")==0)
		return false;

	//No hydrogen atoms between the 2 atoms
	if(!exist_atom(molecule,index1,'H') /*&& !exist_H(molecule,index2)*/) 
		return false;
	//Acceptor saturated doesn't accept HBond
	if(molecule->img[index2].atomType=='N' && (get_H_num(molecule,index2)==1 || get_H_num(molecule,index2)==3) && exist_atom(molecule,index1,'C')) 
		return false;
	// if(exist_H(molecule,index1) && exist_H(molecule,index2) ) //Both atoms has an hydrogen atom
	// 	return false;
	//Atoms are at a distance less than 3
	if(atoms_near(molecule,index1,index2)) 
		return false;
	//if(!(feuille(molecule,index1)) && !(feuille(molecule,index2)) )
		//return false;
	return true;
}

/*========
Check if the atoms are at a distance less than 3 (in terms of separated bonds, i.e. edges in the graph)
========*/
bool atoms_near(struct MyModel* molecule, int index1, int index2){
	if(molecule->covBond[index1][index2]>0) //Atoms are covalently bonded 
		return true;
	else{ //check if the distance between index1 and index2 is greater than 2
		int j=0;
		while (j<index1)
			if(molecule->covBond[j][index1]>0 && molecule->img[j].atomType!='H')
				if(molecule->covBond[j][index2]>0 || molecule->covBond[index2][j]>0 )
					return true;
				else
					j++;
			else
				j++;
		//research by clowns
		while (j<molecule->nbAtom)
			if(molecule->covBond[index1][j]>0 && molecule->img[j].atomType!='H')
				if(molecule->covBond[j][index2]>0 || molecule->covBond[index2][j]>0 )
					return true;
				else
					j++;
			else
				j++;		
	}
	return false;
}

/*========
Check if there is a covalent bond between index and atom with type atomType
========*/
bool exist_atom(struct MyModel* molecule, int index, char atomType){
	int j=0;
	//research by lines
	while (j<index)
		if(molecule->covBond[j][index]==1 && molecule->img[j].atomType==atomType)
			return true;
		else
			j++;
	//research by clowns
	while (j<molecule->nbAtom)
		if(molecule->covBond[index][j]==1 && molecule->img[j].atomType==atomType)
			return true;
		else
			j++;
	return false;
}

/*========
Get the number of hydrogen atoms related to the atom index
========*/
int get_H_num(struct MyModel* molecule, int index){
	int cpt=0;
	int j=0;
	//research by lines
	while (j<index){
		if(molecule->covBond[j][index]!=0 && molecule->img[j].atomType=='H') cpt++;
		j++;
	}
	//research by columns
	while (j<molecule->nbAtom){
		if(molecule->covBond[index][j]!=0 && molecule->img[j].atomType=='H') cpt++;
		j++;
	} 
	return cpt;
}

/*========
Get the index of the atom covalently bonded to hydrogen with index "index"
========*/
int get_cov_index (int **covbond,int size, int index){
	int i=0;
	//research by lines
	while (i<index){
		if(covbond[i][index]==1) 
			return i;
		else
			i++;
	}
	i++;
	//research by columns
	while (i<size){
		if(covbond[index][i]==1 )
			return i;
		else
			i++;
	} 
	return -1;
}


/*========
Check if there is a covalent bond between index and an hydrogen atom
Note: this function is not used any more, it has been generalised to the "exist_atom" function.
========*/
bool exist_H(struct MyModel* molecule, int index){
	int j=0;
	//research by lines
	while (j<index)
		if(molecule->covBond[j][index]==1 && molecule->img[j].atomType=='H')
			return true;
		else
			j++;
	//research by clowns
	while (j<molecule->nbAtom)
		if(molecule->covBond[index][j]==1 && molecule->img[j].atomType=='H')
			return true;
		else
			j++;
	return false;
}

/*========
Find the orientation of (donor,acceptor). Simple H-bond or with proton transfer
Note : last update: March 2017
========*/
int get_Hbond_sens(int **HB,int size, int A, int D){
	int i=0;
	int k=-1;
	while(i<size){
		if(HB[A][i]==-3 && HB[i][i]==D) return 3; //Proton transfer : donor= CB[i][i]
		//HB[A][i]==1 (before)
		if((HB[A][i]==1 || HB[A][i]==0) && (HB[D][i]==-4)) k=2; //Simple H-bond
		i++;
	}
	return k;
}

/*========
Check if the H-bond is formed for a significant period of time along trajectory.
Note : this function is not used yet. It will be used in order to ignore the vibrational motion around a HB.
========*/
bool check_formed_Hbond(struct MyModel *molecule, int index){
	int j=0;
	bool verif=false;
	j=0;
	int resTime=0; //Calculate number of snapshots where H-bond is present
	while(j<molecule->bondHdyn[index].nChange){
		if(j==molecule->bondHdyn[index].nChange-1 && molecule->bondHdyn[index].imgList[j*3+2]==-1){
			molecule->bondHdyn[index].imgList[j*3+2]=molecule->nbMolecule-1;
		}
		if((molecule->bondHdyn[index].imgList[j*3]==1 || molecule->bondHdyn[index].imgList[j*3]==3) && molecule->bondHdyn[index].imgList[j*3+2]-molecule->bondHdyn[index].imgList[j*3+1] >(molecule->tshVal.pourtMinResT*molecule->nbMolecule/100)){
			verif= true; 	

		} 
		resTime =resTime+ (molecule->bondHdyn[index].imgList[j*3+2]-molecule->bondHdyn[index].imgList[j*3+1]);
		// printf("index=%d, j=%d, %d , %d , nbchange=%d , state=%d\n",index, j, molecule->bondHdyn[index].imgList[j*3+2],molecule->bondHdyn[index].imgList[j*3+1],molecule->bondHdyn[index].nChange,molecule->bondHdyn[index].imgList[j*3] );
		//else 
		j++;
	}
	//Residence time is greater than 50% of total time
	// printf("\npourcentage=%d, duree=%d , %.2lf\n",(int)(50.00*molecule->nbMolecule/100), resTime,(molecule->tshVal.pourtMinResT*molecule->nbMolecule/100) );
	if(!verif && resTime>(int)(50.00*molecule->nbMolecule/100)) verif=true; 
	return verif;
}

/*========
Check if the conformation appears at least a percentage of pourt along the trajectory
========*/
bool verif_timeRes_conf(struct confList* conf,double pourt,int size){
	// return true;
	int i=0;
	double resdTime=0.00;
	 //printf("pour=%d\n", (int)(pourt*size/100));
	//printf("\nconf=%d ; nbperiod=%d pour=%.2lf\n",conf->name, conf->imgList[0], pourt*size/100 );
	while(i<conf->imgList[0]){
		resdTime +=(conf->imgList[i*2+2]-conf->imgList[i*2+1]+1);
		i++;
		//printf("period %d = %d\n",i,conf->imgList[i*2+2]-conf->imgList[i*2+1]+1);
		//if((conf->imgList[i*2+2]-conf->imgList[i*2+1]+1)>(int)(pourt*size/100)){
			// printf("name=%d , periode = %d \n",conf->name,conf->imgList[i*2+2]-conf->imgList[i*2+1]+1 );
		//	return true;
		//}	
		//else 
		//	i++;	
	}
	resdTime= (resdTime*100.00)/(double)(size);
	conf->totalPerc = resdTime;
	printf("%d : %.2lf\n",conf->name,resdTime );
	if(resdTime >= pourt) 
		return true;
	else
		return false;
	//return false;

}

/*========
Get the index of the donor from table of donor-hydrogen-acceptor that contains all H-bonds (Update: March 2017)
========*/
int get_donHB(int **donHacc,int index,int size){
	int i=0;
	while(i<size){
		if(donHacc[index][i]==-4)
			return i;
		else
			i++;
	}
	printf("ERROR : donor not found for %d \n",index);
	exit(10);
	return -1;
}

/*========
Get the index of the acceptor from matrix of donor-hydrogen-acceptor (donHacc) that contains all H-bonds (Update: March 2017)
========*/
int get_accHB(int **donHacc,int index,int size){
	int i=0;
	while(i<size){
		if(donHacc[index][i]==-3)
			return i;
		else
			i++;
	}
	printf("ERROR : acceptor not found for %d \n",index);
	exit(10);
	return -1;
}

/*========
Check if the index of the donor is either index1 or index2
========*/
int get_HHB(int **donHacc,int index1, int index2, int size){
	int i=0;
	while(i<size){
	if(i!=index1 && i!=index2){
		if((donHacc[index1][i]==-4 && donHacc[index2][i]!=-1 )|| donHacc[i][i]==index1) //
			return index1;
		if((donHacc[index2][i]==-4 && donHacc[index1][i]!=-1 )|| donHacc[i][i]==index2) // 
			return index2;		
	}	
		i++;
	}
	return -1;
}

/*========
Get the index of the hydrogen atom involved in the HB formed between atoms index1 and index2
========*/
int get_HHHB(struct confList *conf,int index1, int index2, int size){
	int i=0;

	while(i<size){
		if(conf->img[i].atomType=='H'){
			if(conf->HB[index2][i]==-2 /*&& conf->CB[index1][i]>0*/){
				return i;
			}
			if(conf->HB[i][i]==index1){
				return i;
			}
		}
		i++;
	}

	// while(i<size){ //DH = 4 and AH = state
	// 	if(conf->img[i].atomType=='H'){
	// 		if(conf->HB[index1][i]==index1 && (conf->HB[index2][i]==0 || conf->HB[index2][i]==1))
	// 			return i;
	// 		if(conf->HB[index2][i]==index2 && (conf->HB[index1][i]==0 || conf->HB[index1][i]==1))
	// 			return i;
	// 	}
	// 	i++;
	// }
	return -1;
}

/*========
Get the index of the donor and check if there has been a proton transfer in the conformation conf
========*/
bool check_ptr(struct confList *conf,int index1,int index2,int size, int *D){
	int i=0;
	bool ptr=false; //We suppose that there was no proton transfer
	while(i<size){
		if(conf->img[i].atomType=='H'){
			if(conf->HB[index1][i]==index1 /*&& conf->CB[index1][i]>0*/ && conf->HB[index2][i]==-2){
				*D=index1;
			}
			if(conf->HB[index2][i]==index2 /*&& conf->CB[index2][i]>0*/ && conf->HB[index1][i]==-2){
				*D=index2;
			}
			if(conf->HB[i][i]==index1){
				*D=index1;
				ptr=true;
			}
			if(conf->HB[i][i]==index2){
				*D=index2;
				ptr=true;
			}
		}
		i++;
	}
	return ptr;
}

/*========
Get the isomer's name which corresponds to the conformation num_conf.
This function is used in the multiple analysis; when the isomerization analysis is activated.
========*/
int get_isom_num(struct ModelIso *isomList, int confName){
	int i=0;
	while(i<isomList->nbConfF){
		if(isomList->confIsoList[i].confName==confName)	return isomList->confIsoList[i].isoName;
		i++;
	}
	return -1;
}

/*========
Check if the conformation index is in tabState in at least one trajectory
========*/
int verif_type_state(int *tabState,int size,int index){
	int i;
	int cpt=0;
	for (i = 0; i < size; i++) if(tabState[i]==index) cpt++;
	return cpt;
}

/*========
Draw the distribution of the connected components according to their size (number of atoms or water molecules)
========*/
void draw_interface_distribution(char inputFN[], char outputFN[], int size, char *title, char *xlabel, char *ylabel, int type) {
//Draw results using gnuplot file
    FILE* outputD=fopen("distribution.gplt","w");

	fprintf(outputD,"set title '%s'\n",title);
	fprintf(outputD,"\n");

	fprintf(outputD,"set xlabel '%s'\n",xlabel);
	fprintf(outputD,"\n");

	fprintf(outputD,"set ylabel '%s'\n",ylabel);
	fprintf(outputD,"\n");

	fprintf(outputD,"set term png\n");
	fprintf(outputD,"set output '%s'\n",outputFN); // difference between distance's files
	fprintf(outputD,"\n");

	// fprintf(outputD,"set grid\n");

	fprintf(outputD,"set auto x \n");
	fprintf(outputD,"set yrange [0:%d]\n",size);
	// fprintf(outputD,"set xrange [0:%d]\n",size+1);
	fprintf(outputD,"\n");

	// fprintf(outputD, "set boxwidth 0.4 relative\n");
	// fprintf(outputD, "set style data histogram\n");
	// fprintf(outputD, "set style histogram cluster gap 1\n");
	if(type==1 || type==3){
		fprintf(outputD, "set style fill solid border 0 \n");
		fprintf(outputD, "set boxwidth 0.5\n");		
	}
	fprintf(outputD, "set xtic \n");// rotate by -45 scale 0

	// fprintf(outputD, "xoffset=1\n");
	// fprintf(outputD, "set key default\n");
	// fprintf(outputD, "set key font \",8\" \n"); 
	// fprintf(outputD, "set key width -8 \n");


	fprintf(outputD,"plot \\\n");
	// fprintf(outputD,"  '%s' using 1:2 with boxes notitle ,\\\n",inputFN);
	//fprintf(outputD,"  '%s' using ($1+0.2):3  with boxes title '#conf. found',\\\n",inputFN);
	if(type==1){/*xticlabels(3)*/
		fprintf(outputD,"  '%s'    using 1:2:xticlabels(3) with  boxes notitle  lc rgb \"#20B2AA\"",inputFN); //1:($2 + 0.5):($2) , labels : :xticlabels(1)
	}
	if(type==3){
		fprintf(outputD,"  '%s'    using 1:2  with  boxes notitle  lc rgb \"#20B2AA\"",inputFN); //1:($2 + 0.5):($2) , labels : :xticlabels(1)
	}
	if(type==0){
		fprintf(outputD,"  '%s'    using 1:2  with  points notitle pt 7 ps 0.5  lc rgb \"red\"",inputFN); //1:($2 + 0.5):($2) , labels : :xticlabels(1)		
	}
	if(type==2){
		fprintf(outputD,"  '%s'     using 1:2:3 with  points notitle  pt 7 ps 0.5   lc variable",inputFN); //1:($2 + 0.5):($2) , labels : :xticlabels(1)				
	}
	// fprintf(outputD,",  '%s'    using 1:2 with linespoints notitle lc rgb \"#C71585\"  pt 7 ps 0.2 \n",inputFN); //1:($2 + 0.5):($2)
	fprintf(outputD,"\n");
	fprintf(outputD,"\n");
	fclose(outputD);

//Generate the ".png" file
	char cmd[128];
	sprintf(cmd,"gnuplot %s","distribution.gplt");
	system(cmd);
}

/*========
Read energy for one snapshot
========*/
double get_snap_energy(FILE *traj,int size){
 	double energy=0.00 ;
	fscanf(traj,"%d",&size);
	skip_return(traj,1);		
	// char line[512];
	fscanf(traj,"%lf",&energy);
	// fgets(line,512,traj); //comment line
	// energy=get_real_conf_energy(line);//Get the energy of the current snapshot 
	skip_return(traj,size+1); 		
	return  energy; 		
}

/*========
Exctract energy from a line in the xyz file
========*/
double get_real_conf_energy(char *line){
	char enVal[256];
	sprintf(enVal,"%s",strrchr(line,'-'));
	if(strrchr(enVal,',')!=NULL)
		sprintf(enVal,"%s",str_sub(enVal,0,strlen(enVal)-strlen(strrchr(enVal,','))-1));
	
	return(atof(enVal));//*2625.5
	//return (atof(strrchr(line,'-')))*2625.5;
}

/*========
Test isomorphism with Mckay (dreadnaut). Check in the inputFN file if the two graphs are isomorphic or not.
========*/
bool verif_dreadnaut_isom( char inputFN[]){
	// sprintf(resFile,"dreadnaut %s > %s", inputFN,"graph.canonical");
	char resFile[512]="";
	sprintf(resFile,"dreadnaut < %s ", inputFN);
	system(resFile);
	//Read results from graph.isom
	char outputFN[512];
	sprintf(outputFN,"%s.%s",inputFN,"isom");
	FILE *outputF=fopen(outputFN,"r");
	//Read the isomorphism test result from the file .isom
	char line[512];
	fgets(line,512,outputF);
	if(strstr(line,"identical") != NULL){
		fclose(outputF);
		system("rm -f graph.dreadnaut*");
		return true;
	}
	else{
		fclose(outputF);
		system("rm -f graph.dreadnaut*");
		return false;
	}

}

/*========
Save partition of atoms according to their type (C, N, O...)
========*/
void save_partition(FILE *outputF,struct atom *img,int nbAtom){

//Save the partitions for the first graph
	fprintf(outputF, "%s","f=[" );
	int **part=allocate_matrix(NB_MaxAtomType,nbAtom,0,"partitionMatrix"); //Adjacency matrix
	int i;
	for(i=0;i<nbAtom;i++){
		//Oxygen in 0
		if(strcmp(img[i].atomName,"O")==0){
			part[0][0]++;
			part[0][part[0][0]]=i+1;
		}
		//Hydrogen in 1
		if(strcmp(img[i].atomName,"H")==0){
			part[1][0]++;
			part[1][part[1][0]]=i+1;
		}
		//Neitrogen in 2
		if(strcmp(img[i].atomName,"N")==0){
			part[2][0]++;
			part[2][part[2][0]]=i+1;
		}

		//Carbon in 3
		if(strcmp(img[i].atomName,"C")==0){
			part[3][0]++;
			part[3][part[3][0]]=i+1;
		}

		//Lithium 4 
		if(strcmp(img[i].atomName,"Li")==0){
			part[4][0]++;
			part[4][part[4][0]]=i+1;
		}

		//Chlorine 5		
		if(strcmp(img[i].atomName,"Cl")==0){
			part[5][0]++;
			part[5][part[5][0]]=i+1;
		}

		//Sulfur 6		
		if(strcmp(img[i].atomName,"S")==0){
			part[6][0]++;
			part[6][part[6][0]]=i+1;
		}

		//Manganese 7		
		if(strcmp(img[i].atomName,"Mn")==0){
			part[7][0]++;
			part[7][part[7][0]]=i+1;
		}

		//Argon 8		
		if(strcmp(img[i].atomName,"Ar")==0){
			part[8][0]++;
			part[8][part[8][0]]=i+1;
		}
		//Potassium 9		
		if(strcmp(img[i].atomName,"k")==0){
			part[9][0]++;
			part[9][part[9][0]]=i+1;
		}
		//Boron 10		
		if(strcmp(img[i].atomName,"B")==0){
			part[10][0]++;
			part[10][part[10][0]]=i+1;
		}
		//Brom 11		
		if(strcmp(img[i].atomName,"Br")==0){
			part[11][0]++;
			part[11][part[11][0]]=i+1;
		}
		//Silica 12		
		if(strcmp(img[i].atomName,"Si")==0){
			part[12][0]++;
			part[12][part[12][0]]=i+1;
		}
		//Phosphorus 13		
		if(strcmp(img[i].atomName,"P")==0){
			part[13][0]++;
			part[13][part[13][0]]=i+1;
		}
		//Flourine 14	
		if(strcmp(img[i].atomName,"F")==0){
			part[14][0]++;
			part[14][part[14][0]]=i+1;
		}
		//Ruthenium 15	
		if(strcmp(img[i].atomName,"Ru")==0){
			part[15][0]++;
			part[15][part[15][0]]=i+1;
		}
		//Sodium 16	
		if(strcmp(img[i].atomName,"Na")==0){
			part[16][0]++;
			part[16][part[16][0]]=i+1;
		}
		//Iodine 17	
		if(strcmp(img[i].atomName,"I")==0){
			part[17][0]++;
			part[17][part[17][0]]=i+1;
		}

		//Gold 18	
		if(strcmp(img[i].atomName,"Au")==0){
			part[18][0]++;
			part[18][part[18][0]]=i+1;
		}
		//Zinc 19	
		if(strcmp(img[i].atomName,"Zn")==0){
			part[19][0]++;
			part[19][part[19][0]]=i+1;
		}
	}
	for(i=0;i<NB_MaxAtomType;i++){
		if(part[i][0]>0){
			int j;
			for(j=1;j<=part[i][0];j++){
				fprintf(outputF, "%d",part[i][j]);
				if(j<part[i][0])
					fprintf(outputF, "%s",","); //coma between vertices of the same partition
			}
			fprintf(outputF, "%s","|"); //vertical bar between partitions 
		}
	}

	fprintf(outputF, "%s\n","]" );	
	free_matrix(part,NB_MaxAtomType);
}
