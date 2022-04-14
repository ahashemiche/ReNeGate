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
 * \file io.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief I/O treatments.
 * \details This file contains function used to save results using statistical and/or graphical representations. 
 * 	The graphical representations are generated using  the free softwares "graphViz" and "gnuplot".
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
#include "conformation.h"
#include "comfunct.h"
#include "memory.h"
#include "io.h"

/*========================Functions=========================*/

/*========
Save a change (appearance/disappearance of a bond) that has been occurred in the trajectory
========*/
void save_event(struct MyModel *molecule, int index, int index2, char type ,char mode){
	FILE* outputH=NULL;
	if(molecule->nbInputF==1){
		if (mode=='a')	{ //Add an event to the file	
			outputH =fopen(get_FN(molecule->inputDir,molecule->inputFN,".event\0"),"a");
			switch(type){
				//H-bond changes
				case '+' : 
				case '-' :
				case '0' :
					{
					    //Type : '+' if bond created , '-' if broken
						fprintf(outputH,"H-B\t%5d\t%c\t",molecule->num_img,type);
						if(molecule->bondHList[index].sens=='0'){//right sense
							fprintf(outputH, "%s-",format_str(molecule->img[molecule->bondHList[index].numDon].atomType, num_occ(molecule->img,molecule->bondHList[index].numDon , molecule->nbAtom)));  
							fprintf(outputH, "%s-",format_str('H',num_occ(molecule->img,molecule->bondHList[index].numH, molecule->nbAtom)));  
							fprintf(outputH, "%s ", format_str(molecule->img[molecule->bondHList[index].numAcc].atomType, num_occ(molecule->img,molecule->bondHList[index].numAcc, molecule->nbAtom)));  						
						}
						else{
							fprintf(outputH, "%s-", format_str(molecule->img[molecule->bondHList[index].numAcc].atomType, num_occ(molecule->img,molecule->bondHList[index].numAcc, molecule->nbAtom)));  						
							fprintf(outputH, "%s-",format_str('H',num_occ(molecule->img,molecule->bondHList[index].numH, molecule->nbAtom)));  
							fprintf(outputH, "%s ",format_str(molecule->img[molecule->bondHList[index].numDon].atomType, num_occ(molecule->img,molecule->bondHList[index].numDon , molecule->nbAtom)));  
						}
						fprintf(outputH, "\t%c ",molecule->bondHList[index].sens);
						fprintf(outputH, "\t%.2lf\t%.2lf\n",distance(molecule->img[molecule->bondHList[index].numH].x,molecule->img[molecule->bondHList[index].numH].y,molecule->img[molecule->bondHList[index].numH].z,molecule->img[molecule->bondHList[index].numAcc].x,molecule->img[molecule->bondHList[index].numAcc].y,molecule->img[molecule->bondHList[index].numAcc].z),(angle(molecule->img[molecule->bondHList[index].numDon].x,molecule->img[molecule->bondHList[index].numDon].y,molecule->img[molecule->bondHList[index].numDon].z,molecule->img[molecule->bondHList[index].numH].x,molecule->img[molecule->bondHList[index].numH].y,molecule->img[molecule->bondHList[index].numH].z,molecule->img[molecule->bondHList[index].numAcc].x,molecule->img[molecule->bondHList[index].numAcc].y,molecule->img[molecule->bondHList[index].numAcc].z))*180/PI);
					}
				break;
				//Reference snapshot changed
				case 'r' : 
						fprintf(outputH,"IRF\t%5d\t%c\n",molecule->num_img,type);
				break;
				//Proton transfer
				case 'p' : 
					{
						fprintf(outputH,"PTR\t%5d\t%c\t",molecule->num_img,type);
						fprintf(outputH, "%s-",format_str(molecule->img[molecule->bondHList[index].numAcc].atomType, num_occ(molecule->img,molecule->bondHList[index].numAcc , molecule->nbAtom)));  
						fprintf(outputH, "%s-",format_str('H',num_occ(molecule->img,molecule->bondHList[index].numH, molecule->nbAtom)));  					
						fprintf(outputH, "%s ", format_str(molecule->img[molecule->bondHList[index].numDon].atomType, num_occ(molecule->img,molecule->bondHList[index].numDon, molecule->nbAtom)));  
						fprintf(outputH, "\n" );
						//Put the the H-bond is broken in the other sense
						fprintf(outputH,"H-B\t%5d\t%c\t",molecule->num_img,'-');
						fprintf(outputH, "%s-",format_str(molecule->img[molecule->bondHList[index].numAcc].atomType, num_occ(molecule->img,molecule->bondHList[index].numAcc , molecule->nbAtom)));  
						fprintf(outputH, "%s-",format_str('H',num_occ(molecule->img,molecule->bondHList[index].numH, molecule->nbAtom)));  					
						fprintf(outputH, "%s ", format_str(molecule->img[molecule->bondHList[index].numDon].atomType, num_occ(molecule->img,molecule->bondHList[index].numDon, molecule->nbAtom)));  
						if(molecule->bondHList[index].sens=='0')	
							fprintf(outputH, "\t%c ",'1');
						else
							fprintf(outputH, "\t%c ",'0');
						fprintf(outputH, "\t%.2lf\t%.2lf\n",distance(molecule->img[molecule->bondHList[index].numH].x,molecule->img[molecule->bondHList[index].numH].y,molecule->img[molecule->bondHList[index].numH].z,molecule->img[molecule->bondHList[index].numAcc].x,molecule->img[molecule->bondHList[index].numAcc].y,molecule->img[molecule->bondHList[index].numAcc].z),(angle(molecule->img[molecule->bondHList[index].numDon].x,molecule->img[molecule->bondHList[index].numDon].y,molecule->img[molecule->bondHList[index].numDon].z,molecule->img[molecule->bondHList[index].numH].x,molecule->img[molecule->bondHList[index].numH].y,molecule->img[molecule->bondHList[index].numH].z,molecule->img[molecule->bondHList[index].numAcc].x,molecule->img[molecule->bondHList[index].numAcc].y,molecule->img[molecule->bondHList[index].numAcc].z))*180/PI);

					}	
				break;
				//Covalent bond changes
				case '2' : 
				case '3' :
					{
					    //Type : '+' if bond created , '-' if broken
						fprintf(outputH,"COV\t%5d\t",molecule->num_img);
						if(type=='2') fprintf(outputH, "%c\t",'+'); else  fprintf(outputH, "%c\t",'-'); 
						fprintf(outputH, "%s-",format_str(molecule->img[index].atomType, num_occ(molecule->img,index , molecule->nbAtom)));  
						fprintf(outputH, "%s ", format_str(molecule->img[index2].atomType, num_occ(molecule->img,index2, molecule->nbAtom)));  						
						fprintf(outputH, "\t%.2lf\t%d\t%d\n",distance(molecule->img[index].x,molecule->img[index].y,molecule->img[index].z,molecule->img[index2].x,molecule->img[index2].y,molecule->img[index2].z),molecule->img[index].bondMax,molecule->img[index2].bondMax);
					}
				break;
				//Intermolecular bonds changes
				case '4' : 
				case '5' :  
					{
					    //Type : '+' if bond created , '-' if broken
						fprintf(outputH,"ION\t%5d\t",molecule->num_img);
						if(type=='4') fprintf(outputH, "%c\t",'+'); else  fprintf(outputH, "%c\t",'-'); 
						fprintf(outputH, "%s..",format_str(molecule->img[index].atomType, num_occ(molecule->img,index , molecule->nbAtom)));  
						fprintf(outputH, "%s ", format_str(molecule->img[index2].atomType, num_occ(molecule->img,index2, molecule->nbAtom)));  						
						fprintf(outputH, "\t%.2lf\n",distance(molecule->img[index].x,molecule->img[index].y,molecule->img[index].z,molecule->img[index2].x,molecule->img[index2].y,molecule->img[index2].z));
					}
				break;
				//Organometallic bonds changes
				case '6' : 
				case '7' :  
					{
					    //Type : '+' if bond created , '-' if broken
						fprintf(outputH,"METAL\t%5d\t",molecule->num_img);
						if(type=='6') fprintf(outputH, "%c\t",'+'); else  fprintf(outputH, "%c\t",'-'); 
						fprintf(outputH, "%s..",format_str(molecule->img[index].atomType, num_occ(molecule->img,index , molecule->nbAtom)));  
						fprintf(outputH, "%s ", format_str(molecule->img[index2].atomType, num_occ(molecule->img,index2, molecule->nbAtom)));  						
						fprintf(outputH, "\t%.2lf\n",distance(molecule->img[index].x,molecule->img[index].y,molecule->img[index].z,molecule->img[index2].x,molecule->img[index2].y,molecule->img[index2].z));
					}
				break;
			}
		}
		else { outputH=fopen(get_FN(molecule->inputDir,molecule->inputFN,".event\0"),"w"); }
		fclose(outputH);		
	}
}

/*========
Save one snapshot for conformation conf into "file.xyz" file
========*/
void save_snapshot(struct MyModel *molecule, struct confList* conf){
	int i;
	FILE* outputF=NULL;
	char *resFile=NULL;
	resFile = malloc (sizeof (*resFile) * (500));
	if(molecule->nbInputF==1){
		sprintf(resFile,"%s",get_FN(molecule->inputDir,molecule->inputFN,"_xyz/coordinates_conf\0"));
	}
	else{
  		sprintf(resFile,"%s/%s",molecule->inputDir,"conf_xyz/coordinates_conf");
	}
	int cptatom=0;
	for(i=0;i<conf->nbAtom;i++) if(conf->img[i].atomIn) cptatom++; 
	char a[5];
	sprintf(a,"%d",conf->name);
	strcat(resFile,a); 
	strcat(resFile,".xyz\0");	
	outputF =fopen(resFile,"w");
	if (outputF!=NULL){

		fprintf(outputF, "%d\n", cptatom);
		fprintf(outputF, "conformation = %d",conf->name ); 
		//fprintf(outputF, "\t snapshot = %d", molecule->num_img);
		fprintf(outputF, "\n" );
		//fprintf(outputF,"Atom \t x \t\t\t y \t\t\t z \t\t\t  Valence\n");
		for (i = 0; i < conf->nbAtom; i++){
			if(conf->img[i].atomIn)
				fprintf(outputF,"%s \t %lf \t %lf \t %lf \n",conf->img[i].atomName, conf->img[i].x, conf->img[i].y, conf->img[i].z);
		}
	}
	fclose(outputF);
	if(molecule->nbInputF==1 && molecule->sysType[0]!='w'){ //Save the xyz with extra columns 
		sprintf(resFile,"%s",get_FN(molecule->inputDir,molecule->inputFN,"_2_xyz/coordinates_conf\0"));
		char a[5];
		sprintf(a,"%d",conf->name);
		strcat(resFile,a); 
		strcat(resFile,".xyz\0");	
		outputF =fopen(resFile,"w");
		if (outputF!=NULL){
			fprintf(outputF, "%d\n", conf->nbAtom);
			fprintf(outputF, "conformer = %d",conf->name ); 
			fprintf(outputF, "\n" );
			for (i = 0; i < conf->nbAtom; i++)
				fprintf(outputF,"%s \t %lf \t %lf \t %lf \t %d\n",conf->img[i].atomName, conf->img[i].x, conf->img[i].y, conf->img[i].z,check_atom_involved(molecule,conf,i));
		}
		fclose(outputF);
	}
	free(resFile);

}
/*========
Check if atom with index i is involved in dynamic bond or not, 0 if the atom is not involved in a dynamic bond, 1 else.
========*/
int check_atom_involved(struct MyModel *molecule, struct confList* conf, int index){
	int j;
	//Covalent bonds
	if(bit_1(molecule->changType,1)){
		j=0;
		//Search by lines
		while (j<index){
			if(conf->CB[j][index]>0 && molecule->bondDyn[index][j]==1)	
				return 1;
			else 
				j++;
		}
		//Search by columns
		while (j<molecule->nbAtom){
			if(conf->CB[index][j]>0 && molecule->bondDyn[index][j]==1)
				return 1;
			else
				j++;			
		}
	}
	
	//Hydrogen bonds
	if(bit_1(molecule->changType,0)){
		j=0;
		//Search by lines
		while (j<index){
			if(conf->HB[j][index]==-5 && molecule->bondDyn[index][j]==2)	
				return 2;
			else 
				j++;
		}
		//Search by columns
		while (j<molecule->nbAtom){
			if(conf->HB[index][j]==-5 && molecule->bondDyn[index][j]==2)
				return 2;
			else
				j++;			
		}
	}

	//Intermolecular interaction
	if(bit_1(molecule->changType,2)){
		j=0;
		if(strcmp(conf->img[index].atomName,"Ar")==0 || strcmp(conf->img[index].atomName,"K")==0 || strcmp(conf->img[index].atomName,"Li")==0 || strcmp(conf->img[index].atomName,"Na")==0 || strcmp(conf->img[index].atomName,"I")==0 || strcmp(conf->img[index].atomName,"Cl")==0 || strcmp(conf->img[index].atomName,"Br")==0){
			int I = get_index_ion(conf->IB,conf->nbIon,conf->nbAtom,index);
			//Search by columns only 
			while (j<molecule->nbAtom){
				if(conf->IB[I][j]>0 && molecule->bondDyn[index][j]==3)
					return 3;
				else
					j++;			
			}			
		}
		else{
			while(j<conf->nbIon){
				if(conf->IB[j][index]>0 && molecule->bondDyn[index][conf->IB[j][conf->nbAtom]]==3)
					return 3;
				else
					j++;
			}
		}
	}
	//Organometallic bonds 
	if(bit_1(molecule->changType,3)){
		j=0;
		if(strcmp(conf->img[index].atomName,"Mn")==0 || strcmp(conf->img[index].atomName,"Ru")==0 || strcmp(conf->img[index].atomName,"Au")==0){
			int I = get_index_ion(conf->MB,conf->nbMetal,conf->nbAtom,index);
			//Search by columns only 
			while (j<molecule->nbAtom){
				if(conf->MB[I][j]>0 && molecule->bondDyn[index][j]==4)
					return 4;
				else
					j++;			
			}			
		}
		else{
			while(j<conf->nbMetal){
				if(conf->MB[j][index]>0 && molecule->bondDyn[index][conf->MB[j][conf->nbAtom]]==4)
					return 4;
				else
					j++;
			}
		}
	}
	return 0;
}
/*========
Save the covalent bonds molecule->covBond into "file.cov" file
========*/
void save_cov_bonds(struct MyModel *molecule){
	int i= 1 ;  int k=1;
	FILE* outputF=NULL;
	char resFile[1043];
	char rad[256];
	strcpy(rad,molecule->inputFN);
	rad[strlen(molecule->inputFN)-strlen(strrchr(molecule->inputFN,'.'))]='\0';
	molecule->nbCovF++;
	sprintf(resFile,"%s/%s/%s_cov/%s%d.cov",molecule->inputDir,rad,rad,rad,molecule->nbCovF);
	outputF =fopen(resFile,"w");
	if (outputF!=NULL){
		fprintf(outputF,"Bond \t AtomA-atomB \t First atom \t Last atom \t Bond type \t distance\n");
		for ( i = 0; i < molecule->nbAtom; i++){ 
			int j = 0 ;
			for ( j = i+1; j < molecule->nbAtom; j++){
				if (molecule->covBond[i][j]>0){
					fprintf(outputF,"%d \t\t %c-%c \t\t\t %c%d \t\t\t %c%d \t\t %d \t\t\t %lf \n",k,molecule->img[i].atomType, molecule->img[j].atomType, molecule->img[i].atomType , num_occ(molecule->img, i, molecule->nbAtom), molecule->img[j].atomType,num_occ(molecule->img, j, molecule->nbAtom) , molecule->covBond[i][j] ,  distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z));
					k++;
				}
			}
		}
		fclose(outputF);
	}
}

/*========
Save H-bond dynamics into "file.Hbonds" and draw results using Gnuplot
========*/
void save_Hbond_dynamics(struct MyModel *molecule){
	FILE* outputH=NULL;
    int i=0;
	outputH =fopen(get_FN(molecule->inputDir,molecule->inputFN,".Hbonds\0"),"w");
	while(i<molecule->nbHbond){
		int j=0;
		//if(check_formed_Hbond(molecule,i)){
			for(j=0;j<molecule->bondHdyn[i].nChange;j++){
				if(molecule->bondHdyn[i].imgList[j*3]!=-1){
					fprintf(outputH,"%5d\t",molecule->bondHdyn[i].imgList[(j*3)+1]);
					fprintf(outputH,"%2d\t%2d\t%s%d-%s%d-%s%d\t%2d",i+1,molecule->bondHdyn[i].imgList[j*3],molecule->imgRef[molecule->bondHdyn[i].bondH[0]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[0], molecule->nbAtom),molecule->imgRef[molecule->bondHdyn[i].bondH[1]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[1], molecule->nbAtom),molecule->imgRef[molecule->bondHdyn[i].bondH[2]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[2], molecule->nbAtom),molecule->bondHdyn[i].nChange);
					fprintf(outputH,"\n");
					if(j==molecule->bondHdyn[i].nChange-1 && molecule->bondHdyn[i].imgList[j*3+2]==-1){
						fprintf(outputH,"%5d\t",molecule->nbMolecule-1);
						fprintf(outputH,"%2d\t%2d\t%s%d-%s%d-%s%d\t%2d",i+1,molecule->bondHdyn[i].imgList[j*3],molecule->imgRef[molecule->bondHdyn[i].bondH[0]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[0], molecule->nbAtom),molecule->imgRef[molecule->bondHdyn[i].bondH[1]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[1], molecule->nbAtom),molecule->imgRef[molecule->bondHdyn[i].bondH[2]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[2], molecule->nbAtom),molecule->bondHdyn[i].nChange);
						fprintf(outputH,"\n");
					}
					else{
						fprintf(outputH,"%5d\t",molecule->bondHdyn[i].imgList[j*3+2]-1);
						fprintf(outputH,"%2d\t%2d\t%s%d-%s%d-%s%d\t%2d",i+1,molecule->bondHdyn[i].imgList[j*3],molecule->imgRef[molecule->bondHdyn[i].bondH[0]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[0], molecule->nbAtom),molecule->imgRef[molecule->bondHdyn[i].bondH[1]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[1], molecule->nbAtom),molecule->imgRef[molecule->bondHdyn[i].bondH[2]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[2], molecule->nbAtom),molecule->bondHdyn[i].nChange);
						fprintf(outputH,"\n\n");
					}					
				}
			}			
		//	i++;
			fprintf(outputH,"\n");
		//}
		// if(!check_formed_Hbond(molecule,i)){
		// 	update_HbondList(molecule,i);
		// }
		i++;
	}
	fclose(outputH);
//Update the conformations list
	// update_ConfList(molecule);

//Draw the H-bond dynamics using Gnuplot
	plot_Hbonds(molecule);
}

/*========
Save the identified conformations into "file.csv". This file will contain all details about each conformation :
	- Description about the whole trajectory : number of conformations, rotational axes, etc.
	- Dynamics of bonds : which appears.
	- Periods of appearances. 
	- Fragments if a collision. 
========*/
void save_conformations(struct MyModel *molecule){
    //Update conformations (conformation Vs intermediate state)
   	if(molecule->nbConf>1){
    	update_conf_states(molecule);    
   	}
    else{ //if we have only one conformation we consider it as a conformation and not transitional state
    	molecule->conformations->state ='C';
    	molecule->nbTransConf=0;
    }
    printf("save the conformations:\n");
	struct confList* head= molecule->conformations;
	FILE* outputC= NULL;		
	outputC =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_conf.txt\0"),"w");
	//Trajectory description
	fprintf(outputC,"Empirical formula;%s\n", molecule->name);
	fprintf(outputC, "System Type;");
	switch(molecule->sysType[0]){
		case 'i': fprintf(outputC, "MD of an isolated peptide in gas phase\n"); break;
		case 'c': fprintf(outputC, "MD of clusters \n"); break;
		case 'w': fprintf(outputC, "MD of of droplets\n"); break;
		case 'f': fprintf(outputC, "Follow a collision leading to the fragmentation of a peptide\n"); break;
		default : fprintf(outputC, "-\n");
	}
	fprintf(outputC,"Size of trajectory;%d snapshots\n", molecule->nbMolecule);
	fprintf(outputC,"Size of molecule; %d atoms\n",molecule->nbAtom);
	fprintf(outputC, "Number of conformations found;%d conformations\n",molecule->nbConf- molecule->nbTransConf);
	fprintf(outputC, "Number of transitional states found;%d states\n",molecule->nbTransConf);
	//Save the dynamics : 
	int cptCB =0;
	int cptHB =0;
	int cptIB =0;
	int cptMB =0;
	int i; 
	int j;
	if(molecule->changType==0){
		fprintf(outputC, "Dynamics : None\n");
	}
	else{
		if(bit_1(molecule->changType,1) && !molecule->partAnal){ //We don't put the covalent bond if there is an analysis with fragment
			fprintf(outputC, "Covalent bonds:\n");
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==1){
						cptCB++;
						fprintf(outputC,"%s%d-%s%d \t %3d \t %3d\t",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtomDyn),molecule->img[j].atomName, num_occ(molecule->img, j, molecule->nbAtomDyn),i,j);
						add_conf_list(outputC,i,j,1,molecule->conformations);
						fprintf(outputC,"\n");
						
					}
				}
			}
		}
		if(bit_1(molecule->changType,0)){
			fprintf(outputC, "Hydrogen bonds:\n");
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==2){
						cptHB++;
						fprintf(outputC,"%s%d-%s%d \t %3d \t %3d\n",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtomDyn),molecule->img[j].atomName, num_occ(molecule->img, j, molecule->nbAtomDyn),i,j);
					}
				}
			}
		}
		if(bit_1(molecule->changType,2)){
			fprintf(outputC, "Intermolecular bonds:\n");
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==3){
						cptIB++;
						fprintf(outputC,"%s%d-%s%d \t %3d \t %3d\t",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtomDyn),molecule->img[j].atomName, num_occ(molecule->img, j, molecule->nbAtomDyn),i,j);
						add_conf_list(outputC,i,j,3,molecule->conformations);
						fprintf(outputC,"\n");
					}
				}
			}
		}
		if(bit_1(molecule->changType,3)){
			fprintf(outputC, "Organometallic bonds:\n");
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==4){
						cptMB++;
						fprintf(outputC,"%s%d-%s%d \t %3d \t %3d\t",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtomDyn),molecule->img[j].atomName, num_occ(molecule->img, j, molecule->nbAtomDyn), i,j);
						add_conf_list(outputC,i,j,4,molecule->conformations);
						fprintf(outputC,"\n");
					}
				}
			}
		}
		fprintf(outputC,"\n");
	}
	printf("%3d \t %3d \t %3d \t %3d\n",cptCB,cptHB,cptIB,cptMB );
	fprintf(outputC,"\n");
	//Rotational motion 
	// if(molecule->isthList[0]>0){ //there is at least one isthmus (bridge)
	// 	bool confRot = update_isthList(molecule);
	// 	int cpt=0;
	// 	if(confRot){
	// 		fprintf(outputC, "Conformational rotation;");
	// 		int i=0;
	// 		for(i=0;i<molecule->isthList[0];i++){
	// 			if(molecule->isthList[i*3+3]==2){
	// 				fprintf(outputC, "%s%d-%s%d, ", molecule->imgRef[molecule->isthList[i*3+1]].atomName, num_occ(molecule->imgRef, molecule->isthList[i*3+1], molecule->nbAtom),molecule->imgRef[molecule->isthList[i*3+2]].atomName, num_occ(molecule->imgRef, molecule->isthList[i*3+2], molecule->nbAtom) );
	// 				cpt++;
	// 			}
	// 		}
	// 		fprintf(outputC, "\n");
	// 	}
	// 	else{
	// 		fprintf(outputC, "Conformational rotation;-\n");
	// 	}

	// 	if(cpt<molecule->isthList[0]){
	// 		fprintf(outputC, "Simple rotation;");
	// 		int i=0;
	// 		for(i=0;i<molecule->isthList[0];i++){
	// 			if(molecule->isthList[i*3+3]==1){
	// 				fprintf(outputC, "%s%d-%s%d, ", molecule->imgRef[molecule->isthList[i*3+1]].atomName, num_occ(molecule->imgRef, molecule->isthList[i*3+1], molecule->nbAtom),molecule->imgRef[molecule->isthList[i*3+2]].atomName, num_occ(molecule->imgRef, molecule->isthList[i*3+2], molecule->nbAtom) );
	// 				cpt++;
	// 			}
	// 		}
	// 		fprintf(outputC, "\n");		
	// 	}
	// 	else{
	// 		fprintf(outputC, "Simple rotation;-\n");
	// 	}
	// }
	// else{
	// 	fprintf(outputC, "Conformational rotation;-\n");
	// 	fprintf(outputC, "Simple rotation;-\n");		
	// }
	//Conformation description 
	fprintf(outputC, "Conformation");
	fprintf(outputC, "\tType");
	fprintf(outputC, "\t#snapshots");
	fprintf(outputC, "\n");
	// //H-bonds dynamics
	// if(bit_1(molecule->level,0)){
	// 	fprintf(outputC, ";#H-bonds;H-bonds");
	// }
	// //intermolecular bonds dynamics
	// if(bit_1(molecule->level,2)){
	// 	fprintf(outputC, ";#clusters;CN;Intermolecular bonds");
	// }
	// //covalent bonds dynamics
	// if(bit_1(molecule->level,1)){
	// 	fprintf(outputC, ";#fragments");//;covalent for fragments
	// }
	// fprintf(outputC, "\n");
	//FILE *perd=fopen(get_FN(molecule->inputDir,molecule->inputFN,"_perd.data\0"),"w");
	while(head){
			//Save on the file the informations
			fprintf(outputC, "%d", head->name);	
			if(head->state=='C')
				fprintf(outputC, "\t%s","conf.");
			else	
				fprintf(outputC, "\t%s","trans.");

			// for(i=0;i<=head->imgList[0];i++){
				// if(i==0)
					fprintf(outputC, "\t");
				// else 
					// fprintf(outputC, "\t\t");
				//Periods
				// if(i<head->imgList[0]){
					fprintf(outputC,"%d",head->imgList[0]);
				// 	//if(head->state=='C' ) //If one need to save only the stable conformations
				// 		fprintf(perd,"%d \t %d \t %d \n %d \t %d \t %d \n\n",head->imgList[i*2+1],head->name,head->name, head->imgList[i*2+2],head->name, head->name);
				// } 
				// else
				// 	fprintf(outputC,"\t");
				fprintf(outputC, "\n");
			// }
		//List of dynamics bonds 
		// =====	
		if(bit_1(molecule->changType,1) && !molecule->partAnal){ //We don't put the covalent bond if there is an analysis with fragment
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==1){
						if(head->CB[i][j]>0)
							fprintf(outputC,"\t\t CB \t %s%d-%s%d\n",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtomDyn),molecule->img[j].atomName, num_occ(molecule->img, j, molecule->nbAtomDyn));
					}
				}
			}
		}
		if(bit_1(molecule->changType,0)){
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==2){
						if(head->HB[i][j]==-5){
							int D=-1;
							int A=-1;
							bool ptr=check_ptr(head,i,j,head->nbAtom,&D);							
							if(D==i)
								A=j;
							if(D==j)
								A=i;
							if(D==-1){
								printf("Error on the donor (3)\n");
								exit(10);
							}	
							if(ptr){//proton transfer A-H-D
								fprintf(outputC,"\t\t THB \t %s%d-%s%d \n",head->img[A].atomName, num_occ(head->img, A, head->nbAtom),head->img[D].atomName, num_occ(head->img, D, head->nbAtom));								
							}
							else{
								fprintf(outputC,"\t\t HB \t %s%d-%s%d \n",head->img[D].atomName, num_occ(head->img, D, head->nbAtom),head->img[A].atomName, num_occ(head->img, A, head->nbAtom));								
							}
						}	
					}
				}
			}
		}
		if(bit_1(molecule->changType,2)){
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==3){
						int I =-1; 
						int A=-1;
						if(strcmp(head->img[i].atomName,"Ar")==0 || strcmp(head->img[i].atomName,"K")==0 || strcmp(head->img[i].atomName,"Li")==0 || strcmp(head->img[i].atomName,"F")==0 ||  strcmp(head->img[i].atomName,"Na")==0 ||  strcmp(head->img[i].atomName,"I")==0 ||  strcmp(head->img[i].atomName,"Cl")==0 ||  strcmp(head->img[i].atomName,"Br")==0 ){
							I=get_index_ion(head->IB,head->nbIon,head->nbAtom,i);
							
							A=j;
						}
						else{
							I=get_index_ion(head->IB,head->nbIon,head->nbAtom,j);
							A=i;
						}						
						if(I !=-1 && head->IB[I][A]>0){
							fprintf(outputC,"\t\t IB \t %s%d-%s%d \n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
						}
					}
				}
			}
		}
		if(bit_1(molecule->changType,3)){
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==4){
						int I =-1; 
						int A=-1;
						if(strcmp(head->img[i].atomName,"Mn")==0 || strcmp(head->img[i].atomName,"Ru")==0 ||  strcmp(head->img[i].atomName,"Au")==0){
							I=get_index_ion(head->MB,head->nbMetal,head->nbAtom,i);
							A=j;
						}
						else{
							I=get_index_ion(head->MB,head->nbMetal,head->nbAtom,j);
							A=i;
						}
						if(I !=-1 && head->MB[I][A]>0)
							fprintf(outputC,"\t\t OB \t %s%d-%s%d \n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));					
					}
				}
			}
		}
		// =====	
		//next conformation
		head = head->suiv;			
	}
	// fclose(perd);
	fclose(outputC);
//Save time evolution of conformations
	// save_evol_conf(molecule);
//Draw each conformation
	draw_mixed_graph(molecule,'s');
//Draw evolution of conformations
	// plot_conf_periods(molecule);	
//Draw the graph of conformations
	draw_trans_graph(molecule);
}
/*=========Part of the save conformation function
while(head){
			int i=0;
			int j=0;
			int f=0; //Index for covalent fixed
			int h=0;//Index for H-bond present
			int io=0; //Index for Intermolecular bonds
			char **covbond=allocate_char_matrix(molecule->nbCov,15,"covbond");
			//char *cov_Nfix[molecule->nbCov];
			char **hbond=allocate_char_matrix(molecule->nbHbond,30,"hbond");  //molecule->nbHbond
			char **iobond=allocate_char_matrix(NB_MaxInterBond,15,"iobond"); 
			//Description of H-bonds

		if(bit_1(molecule->level,0)){
			for(i=0;i<head->nbAtom;i++){
				for(j=i+1;j<head->nbAtom;j++){ //Update
					if(head->HB[i][j]==-5 && head->img[i].atomType!='H' && head->img[j].atomType!='H'){
						//check if there is proton transfer or not
						int D=-1;
						int A=-1;
						bool ptr=check_ptr(head,i,j,head->nbAtom,&D);							
							if(D==i)
								A=j;
							if(D==j)
								A=i;
							if(D==-1){
								printf("Error on the donor (1) : i=%d , j=%d \n",i,j);
								exit(10);
							}	
						if(ptr){//No prt
							sprintf(hbond[h],"%s%d-%s%d",molecule->imgRef[A].atomName, num_occ(molecule->imgRef, A, molecule->nbAtom),molecule->imgRef[D].atomName, num_occ(molecule->imgRef, D, molecule->nbAtom));
							h++;								
						}
						else{
							sprintf(hbond[h],"%s%d-%s%d",molecule->imgRef[D].atomName, num_occ(molecule->imgRef, D, molecule->nbAtom),molecule->imgRef[A].atomName, num_occ(molecule->imgRef, A, molecule->nbAtom));
							h++;								
						}
					}
				}					
			}						
				
		}
		//Description of intermolecular bonds
		if(bit_1(molecule->level,2)){
			for(i=0;i<molecule->nbIon;i++)
				for(j=0;j<head->nbAtom;j++){
					if(head->IB[i][j]>0 && !((strcmp(head->img[j].atomName,"Ar")==0 || strcmp(head->img[j].atomName,"Li")==0) && get_index_ion(head->IB,molecule->nbIon,head->nbAtom,j) < i)){ //Intermolecular bond
						sprintf(iobond[io],"%s%d-%s%d;", molecule->imgRef[head->IB[i][head->nbAtom]].atomName, num_occ(molecule->imgRef, head->IB[i][head->nbAtom], head->nbAtom),molecule->imgRef[j].atomName, num_occ(molecule->imgRef, j, head->nbAtom));
						io++;
					}
				}
		}
		//Description of covalent bonds	
		if(bit_1(molecule->level,1)){
			for(i=0;i<head->nbAtom;i++)
				for(j=i+1;j<head->nbAtom;j++)
					if(head->CB[i][j]>0 ){ //Check this condition later
						sprintf(covbond[f],"%s%d-%s%d", molecule->imgRef[i].atomName, num_occ(molecule->imgRef, i, head->nbAtom),molecule->imgRef[j].atomName, num_occ(molecule->imgRef, j, head->nbAtom));					
						f++;
					} 
		}			
			
			//Save on the file the informations
			fprintf(outputC, "%d", head->name );	
			if(head->state=='C')
				fprintf(outputC, ";%s","conf.");
			else	
				fprintf(outputC, ";%s","trans.");

			for(i=0;i<=Max(Max(h,head->imgList[0]),io);i++){
				if(i==0)
					fprintf(outputC, ";");
				else 
					fprintf(outputC, ";;");
				//Periods
				if(i<head->imgList[0]){
					fprintf(outputC,"%d-%d;",head->imgList[i*2+1],head->imgList[i*2+2]);
					//if(head->state=='C' ) //If one need to save only the stable conformations
						fprintf(perd,"%d \t %d \t %d \n %d \t %d \t %d \n\n",head->imgList[i*2+1],head->name,head->name, head->imgList[i*2+2],head->name, head->name);

				} 
				else
					fprintf(outputC,";");
				//H-bonds : #H-bonds and list
				if(bit_1(molecule->level,0)){
					if(i==0) 
						fprintf(outputC,"%d;",h);  //#H-bonds
					else
						fprintf(outputC,";");
					if(i<h) 
						fprintf(outputC,"%s;",hbond[i]);
					else
						fprintf(outputC,";");
				}
				// Intermolecular bonds : CN and list
				if(bit_1(molecule->level,2)){
					if(i==0) 
						fprintf(outputC,"%d;%d;",head->nbFrg, io); //#clusters, #CN
					else
						fprintf(outputC,";;");				
					if(i<io) 
						fprintf(outputC,"%s;",iobond[i]);
					else
						fprintf(outputC,";");
				
				}
				//Covalent bonds : #fragments and list for each fragment
				if(bit_1(molecule->level,1)){
					if(i==0) 
						fprintf(outputC,"%d;",head->nbFrg); //#fragments
					else
						fprintf(outputC,";");

					// if(i<f) 
					// 	fprintf(outputC,"%s;",covbond[i]);
					// else
					// 	fprintf(outputC,";");
				}			
				fprintf(outputC, "\n");
			}	
		//next conformation
		head = head->suiv;	
		//Release memory for hydrogen and ionic bonds
		   	free_char_matrix(hbond, molecule->nbHbond);
			free_char_matrix(iobond, NB_MaxInterBond);
			free_char_matrix(covbond,molecule->nbCov);

	}
=========*/

/*========
Add the list of conformers where the bond index1-index2 , of type bondType is present  
========*/
void add_conf_list(FILE *outputF,int index1,int index2, int bondType, struct confList  *conf){
	struct confList *head = conf;
	while(head){
		switch (bondType){
			case 1 : {
				if(head->CB[index1][index2]>0 ||head->CB[index2][index1]>0)
					fprintf(outputF,"%d,",head->name);						
			}
			break;
			case 2 : {
				if(head->HB[index1][index2]==-5 || head->HB[index2][index1]==-5 )
					fprintf(outputF,"%d,",head->name);						
			}
			break;
			case 3 :{
				int I =-1; 
				int A=-1;
					if(strcmp(head->img[index1].atomName,"Ar")==0 || strcmp(head->img[index1].atomName,"K")==0 || strcmp(head->img[index1].atomName,"Li")==0 || strcmp(head->img[index1].atomName,"F")==0 ||  strcmp(head->img[index1].atomName,"Na")==0 ||  strcmp(head->img[index1].atomName,"I")==0 ||  strcmp(head->img[index1].atomName,"Cl")==0 ||  strcmp(head->img[index1].atomName,"Br")==0 ){
					I=get_index_ion(head->IB,head->nbIon,head->nbAtom,index1);
					A=index2;
				}
				else{
					I=get_index_ion(head->IB,head->nbIon,head->nbAtom,index2);
					A=index1;
				}
				if(I !=-1 && head->IB[I][A]>0)
					fprintf(outputF,"%d,",head->name);					
			}
			break;
			case 4:{
				int I =-1; 
				int A=-1;
				if(strcmp(head->img[index1].atomName,"Mn")==0 || strcmp(head->img[index1].atomName,"Ru")==0 || strcmp(head->img[index1].atomName,"Au")==0){
					I=get_index_ion(head->MB,head->nbMetal,head->nbAtom,index1);
					A=index2;
				}
				else{
					I=get_index_ion(head->MB,head->nbMetal,head->nbAtom,index2);
					A=index1;
				}
				if(I !=-1 && head->MB[I][A]>0)
					fprintf(outputF,"%d,",head->name);	
			}
			break;
		
			default:
			break;
		}
		head = head->suiv;
	}
}
/*========
Save the time evolution of conformations (order by periods)
========*/
void save_evol_conf(struct MyModel *molecule){
	struct confList* head = molecule->conformations;
	FILE *output=fopen(get_FN(molecule->inputDir,molecule->inputFN,"_evol.txt\0"),"w");
	//Update links between conformations
	head = molecule->conformations;
	int i=0;
	int j=-1;
	int b = head->imgList[2];	
	while(b < molecule->nbMolecule-1){ //While not the end of the trajectory
		fprintf(output, "%d \t %d \t %d \n", head->imgList[i*2+1],b, head->name);
		j=-1;
		head= get_suiv(molecule,b+1,&j);							
		i=j ;			
		b = head->imgList[i*2+2];
	} 
	fclose(output);
}

/*========
Save list of conformations found in the whole trajectories. 
This function is used in the analysis of multiple trajectories. 
========*/
void save_dyn_desc(struct MyModel* molecule){
  char outputFN[384]; 
  sprintf(outputFN,"%s/%s%s",molecule->inputDir,strrchr(molecule->inputDir,'/'),"_desc.txt");
  FILE *outputC = fopen(outputFN,"w");
  struct confList* head= molecule->conformations;
  //File description
	fprintf(outputC,"Empirical formula; %s\n", molecule->name);
	fprintf(outputC, "System Type; ");
	switch(molecule->sysType[0]){
		case 'i': fprintf(outputC, "MD of an isolated peptide in gas phase\n"); break;
		case 'c': fprintf(outputC, "MD of clusters\n"); break;
		case 'w': fprintf(outputC, "MD of of droplets\n"); break;		
		case 'f': fprintf(outputC, "Follow a collision leading to the fragmentation of a peptide\n"); break;
		default : fprintf(outputC, "-\n");
	}
	//fprintf(outputC,"Size of trajectory;%d snapshots\n", molecule->nbMolecule);
	fprintf(outputC,"Number of trajectories ; %d trajectories\n",molecule->nbInputF);
	fprintf(outputC,"Size of molecule ; %d atoms\n",molecule->nbAtom);
	fprintf(outputC, "Number of conformations found ; %d conformations\n",molecule->nbConf);

	//Nature of dynamics
	fprintf(outputC, "Nature of dynamics;");
	if(molecule->changType==0){
		fprintf(outputC, "-\n");
	}
	else{
		if(bit_1(molecule->changType,0))
			fprintf(outputC, "H-bonds, ");
		if(bit_1(molecule->changType,1))
			fprintf(outputC, "covalent bonds, ");
		if(bit_1(molecule->changType,2))
			fprintf(outputC, "intermolecular bonds, ");
		if(bit_1(molecule->changType,3))
			fprintf(outputC, "organometallic bonds, ");
		fprintf(outputC,"\n");
	}
	//Rotational motion 
	if(molecule->isthList[0]>0){ //there is at least one isthme
		bool confRot = update_isthList(molecule);
		int cpt=0;
		if(confRot){
			fprintf(outputC, "Conformational rotation;");
			int i=0;
			for(i=0;i<molecule->isthList[0];i++){
				if(molecule->isthList[i*3+3]==2){
					fprintf(outputC, "%s%d-%s%d, ", molecule->imgRef[molecule->isthList[i*3+1]].atomName, num_occ(molecule->imgRef, molecule->isthList[i*3+1], molecule->nbAtom),molecule->imgRef[molecule->isthList[i*3+2]].atomName, num_occ(molecule->imgRef, molecule->isthList[i*3+2], molecule->nbAtom) );
					cpt++;
				}
			}
			fprintf(outputC, "\n");
		}
		else{
			fprintf(outputC, "Conformational rotation;-\n");
		}

		if(cpt<molecule->isthList[0]){
			fprintf(outputC, "Simple rotation;");
			int i=0;
			for(i=0;i<molecule->isthList[0];i++){
				if(molecule->isthList[i*3+3]==1){
					fprintf(outputC, "%s%d-%s%d, ", molecule->imgRef[molecule->isthList[i*3+1]].atomName, num_occ(molecule->imgRef, molecule->isthList[i*3+1], molecule->nbAtom),molecule->imgRef[molecule->isthList[i*3+2]].atomName, num_occ(molecule->imgRef, molecule->isthList[i*3+2], molecule->nbAtom) );
					cpt++;
				}
			}
			fprintf(outputC, "\n");		
		}
		else{
			fprintf(outputC, "Simple rotation;-\n");
		}
	}
	else{
		fprintf(outputC, "Conformational rotation;-\n");
		fprintf(outputC, "Simple rotation;-\n");		
	}
	while(head){
		//Save the cartesian coordinates of current conformation
		save_snapshot(molecule, head);
		//Save the conformation bonds 
		save_conf_bonds(molecule, head,'c');
		head = head->suiv;
	}
	fclose(outputC);
}

/*========
Save the bonds of a conformation conf into "conformation/bonds_conf/conf.txt" file. 
This is an intermediate file used in the multiple analysis (analysis of many trajectories simultaneously) 
========*/
void save_conf_bonds(struct MyModel* molecule,struct confList* conf, char type ){
	int i;
	int j;
	FILE* outputF=NULL;
	//char resFile[256]="";
	char *resFile=NULL;
	resFile = malloc (sizeof (*resFile) * (500));
	if(type=='c'){
		sprintf(resFile,"%s/%s",molecule->inputDir,"conformations/bonds_conf");
	}
	else{
		sprintf(resFile,"%s",get_FN(molecule->inputDir,molecule->inputFN,"_bonds/bonds_conf\0"));
	}
	char a[5];
	sprintf(a,"%d",conf->name);
	strcat(resFile,a); 
	strcat(resFile,".txt\0");	
	outputF =fopen(resFile,"w");
	if (outputF!=NULL){
		fprintf(outputF, "%d\n",conf->name ); 
		// fprintf(outputF, "%d\n", conf->nbFrg);
		//Covalent bonds (1 : single , 2: double, 3 : triple)
	    for(i=0;i<conf->nbAtom;i++){
 	        for(j=i+1;j<conf->nbAtom;j++){
	          if(conf->CB[i][j]>0){ //There is a covalent bond
	            fprintf(outputF,"%1d\t%3d\t%3d\n",1, i,j); // conf->CB[i][j]    
	          } 
	        }	    	
	    }
		//Hydorgen bonds (4: simple hydrogen bond, 5: proton transfer)
	    if(bit_1(molecule->level,0)){
	      for(i=0;i<conf->nbAtom;i++)
	        for(j=i+1;j<conf->nbAtom;j++){
				if(conf->HB[i][j]==-5 && conf->img[i].atomType!='H' && conf->img[j].atomType!='H'){
					//check if there is proton transfer or not
					int D=-1;
					int A=-1;
					bool ptr=check_ptr(conf,i,j,conf->nbAtom,&D);							
						if(D==i)
							A=j;
						if(D==j)
							A=i;
						if(D==-1){
							printf("Error on the donor (2)\n");
							exit(10);
						}	
					if(ptr){//proton transfer A-H-D
					    fprintf(outputF,"%1d\t%3d\t%3d\t%3d\n",5,D,A,get_HHHB(conf, D,A,conf->nbAtom));
					}
					else{//D-H-A
		                fprintf(outputF,"%1d\t%3d\t%3d\t%3d\n",4,D,A,get_HHHB(conf, D,A,conf->nbAtom));
					}
				}
	        }
	    }
		//Intermolecular bonds (7)
		if(bit_1(molecule->level,2)){
			for(i=0;i<molecule->nbIon;i++)		
		      for(j=0;j<molecule->nbAtom;j++)
		        if(conf->IB[i][j]>0){
		          fprintf(outputF,"%1d\t%3d\t%3d\n", 7, conf->IB[i][molecule->nbAtom],j);
		        }
	    }	
	    //Organometallic bonds (8)
	    if(bit_1(molecule->level,3)){
			for(i=0;i<molecule->nbMetal;i++)		
		      for(j=0;j<molecule->nbAtom;j++)
		        if(conf->MB[i][j]>0){
		          fprintf(outputF,"%1d\t%3d\t%3d\n", 8, conf->MB[i][molecule->nbAtom],j);
		        }	    	
	    }
		fclose(outputF);
	    free(resFile);
	}	
}

/*========
Save isomers list into "conf2iso.txt" file.
========*/
void save_isom_bonds(struct ModelIso *isomList){
	int i;
	int j;
	FILE* outputI=NULL;
	char resFileI[515]="";
	sprintf(resFileI,"%s%s",isomList->inputDir,"/conf2iso.txt");
	outputI =fopen(resFileI,"w");

	for(i=1;i<=isomList->nbIso;i++){
		fprintf(outputI, "%2d",i );
		int k=0;
		while(isomList->confIsoList[k].isoName!=i) k++;
		//Save bonds of the current isomer
		FILE* outputF=NULL;
		char resFile[520]="";
		sprintf(resFile,"%s%s",isomList->inputDir,"/isomers/bonds_isom");
		char a[5];
		sprintf(a,"%d",i);
		strcat(resFile,a); 
		strcat(resFile,".txt\0");	
		outputF =fopen(resFile,"w");
		fprintf(outputF, "%d\n",i ); 
		fprintf(outputF, "%d\n",isomList->confIsoList[k].nbBonds ); 
		fprintf(outputF, "%d\n", isomList->confIsoList[k].nbfrg);
	    for(j=0;j<isomList->confIsoList[k].nbBonds;j++){
	        fprintf(outputF,"%1d\t%3d\t%3d\n",1,isomList->confIsoList[k].tabBonds[j*2],isomList->confIsoList[k].tabBonds[j*2+1] );         
	    }
		fclose(outputF);
		//Save the Cartesian coordinates of the current isomer
		save_isom_xyz(isomList->inputDir,i,isomList->confIsoList[k].confName);
		//Get the list of conformations that belongs to the current isomer i 
		int cpt=0;
		while(k<isomList->nbConfF){
			if(isomList->confIsoList[k].isoName==i)	cpt++;
			k++;
		}
		fprintf(outputI, ";%d;",cpt );
		k=0;
		while(k<isomList->nbConfF){
			if(isomList->confIsoList[k].isoName==i){
				fprintf(outputI, "%d,",isomList->confIsoList[k].confName);
			}
			k++;
		}
		fprintf(outputI, "\n");
	}
	fclose(outputI);
}

/*========
Save Cartesian coordinates for isomer isom_num that corresponds to conformation conf_num. 
========*/
void save_isom_xyz(char inputDir[500], int isom_num, int conf_num){
	//Open files
	FILE* outputC=NULL;
	char resFile[256]="";
	sprintf(resFile,"%s%s",inputDir,"/conf_xyz/coordinates_conf");
	char a[5];
	sprintf(a,"%d",conf_num);
	strcat(resFile,a); 
	strcat(resFile,".xyz\0");	
	outputC =fopen(resFile,"r");

	FILE* outputI=NULL;
	sprintf(resFile,"%s%s",inputDir,"/isom_xyz/coordinates_isom");
	sprintf(a,"%d",isom_num);
	strcat(resFile,a); 
	strcat(resFile,".xyz\0");	
	outputI =fopen(resFile,"w");
	//Copy the content of outputC on outputI
	char line[512];
	//Number of atoms
	fgets(line,512,outputC);
	fputs(line,outputI);
	//Name of the isomer
	fgets(line,512,outputC);	
	fprintf(outputI, "isomer=%d\n",isom_num);
	while(fgets(line,512,outputC))	fputs(line,outputI);
	//Close files
	fclose(outputC);
	fclose(outputI);
}

/*========
Save the dynamics of isomers for each trajectory. 
The isomers are identified based upon change in covalent bonds.
========*/
void save_isom_traj(struct ModelIso *isomList){
  	DIR* rep = NULL;
  	char filPath[512];
  	struct dirent* TrajConFile = NULL; 
//Open dir
    sprintf(filPath,"%s",isomList->inputDir);  
	strcat(filPath,"/traj_conf");
  	rep =opendir(filPath);
//Get list of conformations and update table
  	while ((TrajConFile = readdir(rep)) != NULL){
  		if (strcmp(TrajConFile->d_name, ".") != 0 && strcmp(TrajConFile->d_name, "..") != 0){
			char confTrajFN[1024]="";
			sprintf(confTrajFN,"%s/%s",filPath,TrajConFile->d_name);			
  			//Get the isomer file name
			char resFile[512]="";
			sprintf(resFile,"%s%s",isomList->inputDir,"/traj_isom/");
			char isomTrajFN[256]="";
			sprintf(isomTrajFN,"%s",TrajConFile->d_name);
			isomTrajFN[strlen(isomTrajFN)-strlen(".conf")]='\0';
			strcat(isomTrajFN,".isom");
			strcat(resFile,isomTrajFN);
			//Open files
			FILE* outputC=fopen(confTrajFN,"r"); //Traj_conf
			FILE* outputI=fopen(resFile,"w"); //Traj_isom
			//Copy the content of outputC on outputI
			char line[512];
			//Name of trajectory
			fgets(line,512,outputC);
			fputs(line,outputI);
			//#snapshots
			fgets(line,512,outputC);
			fputs(line,outputI);
			//#atoms
			fgets(line,512,outputC);
			fputs(line,outputI);
			//Dynamics
			int c;
			int prev_isom =-1;
			int num_img =-1;
			int num_conf =-1;
		   	while((c=fgetc(outputC))!=EOF){
		   		ungetc(c,outputC);
				//Read the number of snapshot
				if(fscanf(outputC,"%d", &num_img)!=1){
			      	printf("can not read the number of the snapshot\n");
			      	exit(-1);
			   	}   		
				//Read the conformation
				if(fscanf(outputC,"%d", &num_conf)!=1){
			      	printf("can not read the name of the conformation\n");
			      	exit(-1);
			   	}   		
			   	//Save if the current snapshot has a different isomer
			   	int num_isom = get_isom_num(isomList,num_conf);
			   	if(num_isom != prev_isom){
			   		fprintf(outputI, "%d \t %d\n", num_img,num_isom);
			   		prev_isom = num_isom;
			   	}
			   	skip_return(outputC,1); //go to the next line 
		   	}

			//Close files
			fclose(outputC);
			fclose(outputI);

		}
	}
//Close dir
  	closedir(rep);
}

/*========
Save the graphs of transitions (different paths explored) for the trajectories analysed.
The conformations explored are divided to three types :
	- Initial state. 
	- Intermediate state.
	- Final state.
========*/
void draw_globalGraph(struct MyTrajGraph graphTraj){
	struct MyGraphList* head=graphTraj.graphList;
	char graphTrajFN[550];
	char fstateTrajFN[550];
	char lstateTrajFN[550];
	char istateTrajFN[550];

	FILE* outputG; //File to save global graphs
	FILE* outputFC; //File to save initial states with trajectories
	FILE* outputLC; //File to save final states with trajectories
	FILE* outputIC; //File to save intermediate states with trajectories

	if(head->stateG->dynType=='c'){
		sprintf(graphTrajFN,"%s/%s/%s",head->stateG->inputDir,"global_conf_graph","graphsTrajList.txt");
		sprintf(fstateTrajFN,"%s/%s/%s",head->stateG->inputDir,"conf_graph","fConfTrajList.txt");
		sprintf(lstateTrajFN,"%s/%s/%s",head->stateG->inputDir,"conf_graph","lConfTrajList.txt");
		sprintf(istateTrajFN,"%s/%s/%s",head->stateG->inputDir,"conf_graph","iConfTrajList.txt");

		outputFC=fopen(fstateTrajFN,"w");
		outputLC=fopen(lstateTrajFN,"w");
		outputIC=fopen(istateTrajFN,"w");

		fprintf(outputFC, "Initial isomer\t #occurrence\t Trajectories \n");
		fprintf(outputLC, "Final isomer\t #occurrence\t Trajectories \n");
		fprintf(outputIC, "Intermediate isomer\t Trajectories\t #occurrence \n");
	}
	else{
		sprintf(graphTrajFN,"%s/%s/%s",head->stateG->inputDir,"global_isom_graph","graphsTrajList.txt");
		sprintf(fstateTrajFN,"%s/%s/%s",head->stateG->inputDir,"isom_graph","fIsomTrajList.txt");
		sprintf(lstateTrajFN,"%s/%s/%s",head->stateG->inputDir,"isom_graph","lIsomTrajList.txt");
		sprintf(istateTrajFN,"%s/%s/%s",head->stateG->inputDir,"isom_graph","iIsomTrajList.txt");

		outputFC=fopen(fstateTrajFN,"w");
		outputLC=fopen(lstateTrajFN,"w");
		outputIC=fopen(istateTrajFN,"w");

		fprintf(outputFC, "Initial isomer\t #occurrence\t Trajectories \n");
		fprintf(outputLC, "Final isomer\t #occurrence\t Trajectories \n");
		fprintf(outputIC, "Intermediate isomer\t Trajectories\t #occurrence \n");

	}
	//Save type of states: initial, intermediate or final 
	int i;
	for(i=1;i<=head->stateG->nbdynm;i++){
		bool interState = true;
		if(verif_type_state(graphTraj.fState,graphTraj.nbInputF,i)>0){
			fprintf(outputFC, "%d\t%d\t",i,verif_type_state(graphTraj.fState,graphTraj.nbInputF,i) );
			for(int j=0;j<graphTraj.nbInputF; j++) if(graphTraj.fState[j]==i) fprintf(outputFC,"%s, ",graphTraj.files[j]);
			fprintf(outputFC, "\n");
			interState= false;
		}
		if(verif_type_state(graphTraj.lState,graphTraj.nbInputF,i)>0){
			fprintf(outputLC, "%d\t%d\t",i, verif_type_state(graphTraj.lState,graphTraj.nbInputF,i) );
			for(int j=0;j<graphTraj.nbInputF; j++) if(graphTraj.lState[j]==i) fprintf(outputLC,"%s, ",graphTraj.files[j]);
			fprintf(outputLC, "\n");
			interState= false;
		}
		if(interState){ //Current state is an intermediate conformation		
			int cpt=0;
			fprintf(outputIC, "%d\t",i);
			struct MyGraphList* currGraph=graphTraj.graphList;
			while(currGraph){
				//browse the graph 
				struct MyGraph *graph = currGraph->stateG;
				for(int j=0;j<=graph->nbdynm;j++){
					if(graph->adjMatrix[i][j]>0 ){
						for(int k=0;k<graphTraj.nbInputF; k++){//Browse the graphs
							if(graphTraj.filegraph[k]==currGraph->name){
								fprintf(outputIC,"%s, ",graphTraj.files[k]);	
								cpt++;			
							} 
						} 
					}
				}
				currGraph = currGraph->suiv;
			}
			fprintf(outputIC, "\t %d \n",cpt);
		}
	}
	//close file
	fclose(outputFC);
	fclose(outputLC);
	fclose(outputIC);

	outputG=fopen(graphTrajFN,"w");
	fprintf(outputG, "Global graph;Trajectories\n");
	//Draw global graphs
	while(head){
      //Draw the corresponding graph 
		fprintf(outputG, "%d;",head->name);
		//Save trajectories name file that contains head->name global graph
		int i; 	for(i=0;i<graphTraj.nbInputF; i++) if(graphTraj.filegraph[i]==head->name) fprintf(outputG,"%s, ",graphTraj.files[i]);
		fprintf(outputG, "\n");
		if(head->stateG->dynType=='c')		
			sprintf(graphTrajFN,"%s/%s/%s%d%s",head->stateG->inputDir,"global_conf_graph","globalGraph",head->name,".png");
      	else
			sprintf(graphTrajFN,"%s/%s/%s%d%s",head->stateG->inputDir,"global_isom_graph","globalGraph",head->name,".png");

      	draw_trajGraph(head->stateG,graphTrajFN, -1); //single      
		head = head->suiv; 
	}
	//close file
	fclose(outputG);

}

/*========
Draw the graph of transition file (GraphViz file) related to trajectory num_traj
========*/
void draw_trajGraph(struct MyGraph *graph, char graphTrajFN[256], int num_traj){
//Save the graph of transitions between conformations/isomers
	FILE* outputG=NULL;
	FILE* outputD=NULL;
	char trajgraph[535]="";
	char trajdynm[535]="";
	int i;
	int j;

	if(graph->dynType=='c'){
    	sprintf(trajgraph,"%s/%s",graph->inputDir,"conf_graph.gv");
  		sprintf(trajdynm,"%s/%s/%s",graph->inputDir,"traj_conf_graph","trajConfList.txt");
  		if(num_traj==1){
  			outputD=fopen(trajdynm,"w");
  			fprintf(outputD, "Trajectory;#snaphots;Conformations explored;#conformations\n");
  		}
  		else
  			outputD=fopen(trajdynm,"a");
	}
	else{
    	sprintf(trajgraph,"%s/%s",graph->inputDir,"isom_graph.gv");
  		sprintf(trajdynm,"%s/%s/%s",graph->inputDir,"traj_isom_graph","trajIsomList.txt");
  		if(num_traj==1){
  	  		outputD=fopen(trajdynm,"w");
	  		fprintf(outputD, "Trajectory;#snaphots;Isomers explored;#isomers\n");		
  		}
  		else
  			outputD=fopen(trajdynm,"a");
	}
	//Initialization of vector of conformations/isomers explored at 0;
	int dynExplored[graph->nbdynm]; //Save conformations/isomers explored along the trajectory
	for(i=0;i<=graph->nbdynm;i++) dynExplored[i]=0;
	outputG =fopen(trajgraph,"w");
	//Graph with frequencies
	fprintf(outputG, "digraph G {\n");
	fprintf(outputG, "label=\"Graph with frequencies\";\n");
	fprintf(outputG, "node [style=filled];\n");
	fprintf(outputG, "graph [bgcolor=transparent];\n");
	fprintf(outputG, "node [shape = circle];\n");
	//Browse the adjacency matrix
	for(i=1;i<=graph->nbdynm;i++){
		if(graph->adjMatrix[i][i]==-1 || graph->adjMatrix[i][i]==-3){ //Initial state
			fprintf(outputG, "%d[fillcolor=yellow]\n", i);//initial state, {t[shape=point]}->
			dynExplored[i]=1;
		}
		if(graph->adjMatrix[i][i]==-2 || graph->adjMatrix[i][i]==-3){ //Final state
			fprintf(outputG, "%d[shape = doublecircle];\n",i); //final state
			dynExplored[i]=1;
		}
		for(j=1;j<=graph->nbdynm;j++){
			if(graph->adjMatrix[i][j]>0){
				if(num_traj>=0)
					fprintf(outputG, "%d->%d[label=%d];\n", i,j,graph->adjMatrix[i][j]);
				else
					fprintf(outputG, "%d->%d;\n", i,j);
				dynExplored[i]=1;
				dynExplored[j]=1;
			}
		}
	}	
	fprintf(outputG, "}\n");
	fclose(outputG);
	if(num_traj>0){
		//Save name of the trajectory and its size
		fprintf(outputD, "%s;",graph->inputFN );
		fprintf(outputD, "%d;",graph->nbsnapshots);
		//Save conformations/isomers explored along the trajectory
		for(i=1;i<=graph->nbdynm;i++){
			if(dynExplored[i]==1){
				fprintf(outputD, "%d, ",i );
				dynExplored[0]++;
			}
		}
		fprintf(outputD, ";%d\n", dynExplored[0]);
	}
	fclose(outputD);

//Draw the graph of transition with frequencies using graphViz
	char resFile[640];
	sprintf(resFile,"%s %s -o %s ","dot -Tpng  ",trajgraph,graphTrajFN);
	system(resFile);
}

/*========
Save the connected components found and their size.
========*/
void save_cc(char *outputFN,struct nodCC* tabCC,int size, int *maxCC){
	int i=0;
	FILE* outputF=fopen(outputFN,"w");
	if(outputF == NULL){ printf("Can\'not open file %s \n",outputFN); exit(4);}  

	int ccLabel=1;
	for(i=1; i<=size; i++){
		if(tabCC[i].size!=0){
			//Save  label connected component (cpt number and real one) , size of the connected component
			fprintf(outputF, "%d \t %d \t %d \n",ccLabel, i, tabCC[i].size);
			ccLabel++;
		}
	}
	*maxCC = Max(*maxCC, ccLabel);
	fclose(outputF);
}

/*========
Save global values for graph of possible conformations
========*/
void save_global_values(char *outputFN,struct GraphModel* graphPoss){
	char outputVal[256];
	sprintf(outputVal,"%s_energy.val",outputFN);	
	FILE* outputF=fopen(outputVal,"w");
	if(outputF == NULL){ printf("Can\'not open file _energy.val \n"); exit(4);}  
	
	//#minEng : 
	fprintf(outputF,"%d \n",graphPoss->minEng);
	//#maxEng : 
	fprintf(outputF,"%d \n",graphPoss->maxEng);
	//#nbValEng : 
	fprintf(outputF,"%d \n",graphPoss->nbValEng);
	//#nbConf : 
	fprintf(outputF,"%d \n",graphPoss->nbConfT );
	//#maxCCLabel : 
	fprintf(outputF,"%d \n",graphPoss->ccLabel);
	//#maxCC : 	
	fprintf(outputF,"%d \n",graphPoss->maxCC);
	//#nbEdge : 	
	fprintf(outputF,"%d \n",graphPoss->nbEdge);
	//Close file
	fclose(outputF);

}

/*========
Save the conformations generated in the current level
========*/
void save_level(char *outputFN,struct conf_elt* level ){
	FILE *conFile=fopen(outputFN,"a");
	if(conFile == NULL){ printf("Can\'not open file  %s\n",outputFN); exit(4);}  

	struct conf_elt* head=level;
	while(head){
		fprintf(conFile, "%s \t %d\n",head->conf, head->index );
		head= head->suiv;
	}
	fclose(conFile);
}

/*========
Save a generated conformation
========*/
void save_conf(FILE *conFile,struct conf_elt* levelElt ){
	//FILE *conFile=fopen(outputFN,"a");
	if(conFile == NULL){ printf("Can\'not open file  \n"); exit(4);}  
	fprintf(conFile, "%s \t %d\n",levelElt->conf, levelElt->index );
	//fclose(conFile);
}

/*========
Save the real energy of conformations found in simulation, from the xyz file.
========*/
void save_real_energy(struct MyModel* molecule){
	
	printf("Get energy\n");
	char inputF[512]; 
	sprintf(inputF,"%s/%s",molecule->inputDir,molecule->inputFN);
	
	FILE* outputC =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_energy.txt\0"),"w");
	if(outputC== NULL){ printf("Can\'not open file _energy.txt\n"); exit(4);}

	FILE* outputP =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_periods.dat\0"),"r");
	if(outputP== NULL){ printf("Can\'not read file _periods.dat\n"); exit(4);}
	
	FILE *traj;
	if((traj=fopen(inputF,"r"))== NULL){ printf("Can\'not read file .xyz\n"); exit(3);}
	int confName=0;
	int snap_str = 0;
	int snap_end =0;
	double engMin[molecule->nbConf];
	int snapMin[molecule->nbConf];
	int confreq[molecule->nbConf];
	init_tab_double(engMin,molecule->nbConf,1000.00);
	init_tab(snapMin, molecule->nbConf,-1);
	init_tab(confreq, molecule->nbConf,0);
	bool finish = false; 
	while(!finish){
		//Read the conf name
		if(fscanf(outputP,"%d", &confName)!=1){
			printf("can not read the confNamex \n");
			exit(-1);
		}   		
		//Read the start of the period
		if(fscanf(outputP,"%d", &snap_str)!=1){
			printf("can not read the start of the period\n");
			exit(-1);
		}   	
		//Read the end of the period	
		if(fscanf(outputP,"%d", &snap_end)!=1){
			printf("can not read the end of the period\n");
			exit(-1);
		} 
		skip_return(outputP,1); 
		printf("%d %d %d\n",confName,snap_str,snap_end); 
		//read the energy for all the period 
		int i= snap_str;
		for (i = snap_str; i <= snap_end; i++){
			double eng= get_snap_energy(traj,molecule->nbAtom);
			if(engMin[confName-1]>eng){
				engMin[confName-1]=eng;
				snapMin[confName-1]=i+1;
			}
			confreq[confName-1]++;
			fprintf(outputC,"%d \t %d \t %lf\n",confName,i+1,eng); 		
		}
		if(snap_end==molecule->nbMolecule-1) finish = true;
	}
	fclose(outputC);
	fclose(outputP);
	struct confList* head=molecule->conformations;
	while(head){
		head->snapMin=snapMin[head->name-1];
		head->engMin=engMin[head->name-1];
		head->totalPerc=(double)(confreq[head->name-1])/(double)(molecule->nbMolecule);
		head = head->suiv;
	}
}

/*========
Compare the different conformers to the reference in terms of bonds/interactions
========*/
void save_diff_with_ref(struct MyModel* molecule){
	
	printf("Get the comparaison with the reference structure\n");
	
	int i;
	int j;

	struct confList* confRef= molecule->conformations;
	struct confList* head= molecule->conformations->suiv;
	// char *confLab;
	
	FILE* outputC =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_comp_to_ref.gv\0"),"w");
	if(outputC== NULL){ printf("Can\'not open file _comp_to_ref.gv\n"); exit(4);}

	fprintf(outputC, "digraph G {\n");
	fprintf(outputC, "rankdir=\"LR\";\n");
	fprintf(outputC, "node [style=filled,  fontcolor=black, fontname=\"blod, bold\"];\n");
	char rad[256]; 
	strcpy(rad,molecule->inputFN);
	rad[strlen(molecule->inputFN)-strlen(strrchr(molecule->inputFN,'.'))]='\0';
	fprintf(outputC, "\"%s-%d\" ", rad,confRef->name);
	fprintf(outputC, "[fillcolor=aquamarine1, fontcolor=black, fontname=\"blod, bold\", CA=%d, image=\"%s%d.png\", ",get_coord_num(confRef,0),get_FN(molecule->inputDir,molecule->inputFN,"_conf/conformation\0"),confRef->name);
	fprintf(outputC,"engMin=%lf, snapMin=%d, ",confRef->engMin,confRef->snapMin);
	fprintf(outputC,"freq=%.2lf, ",confRef->totalPerc);
	// printf("engMin=%lf, snapMin=%d, \n",confRef->engMin,confRef->snapMin);
	fprintf(outputC," label=\" %s \"]\n",rad);
	// get_label_conf(confRef->name,molecule->nbMolecule, get_FN(molecule->inputDir,molecule->inputFN,"_energy.txt\0"));
	//Compare the different conformers to the reference in terms of bonds/interactions
	while (head){
		fprintf(outputC, "\"%s-%d\"", rad,head->name);
		fprintf(outputC, "[fillcolor=grey92,fontcolor=black, fontname=\"blod, bold\",CA=%d, image=\"%s%d.png\"",get_coord_num(head,0), get_FN(molecule->inputDir,molecule->inputFN,"_conf/conformation\0"),head->name);
		fprintf(outputC,"engMin=%lf, snapMin=%d, ",head->engMin,head->snapMin);
		fprintf(outputC,"freq=%.2lf, ",head->totalPerc);
		// confLab= get_label_conf(head->name,molecule->nbMolecule, get_FN(molecule->inputDir,molecule->inputFN,"_energy.txt\0"));
		fprintf(outputC," label=\" %s \"]\n",get_label_conf(head->name,molecule->nbMolecule, get_FN(molecule->inputDir,molecule->inputFN,"_energy.txt\0")));
		// free(confLab);
		fprintf(outputC,"\"%s-%d\" -> \"%s-%d\" [color=black, label=\"",rad,confRef->name,rad,head->name);
		
		//Browse the different bonds and check the differences 
		//Covalent bonds
		if(bit_1(molecule->level,1)){
			for(i=0;i<confRef->nbAtom;i++){
				for(j=i+1;j<confRef->nbAtom;j++){
					if(head->CB[i][j]>0 &&  confRef->CB[i][j]<=0){ //(bond formed)
						fprintf(outputC,"+%s%d-%s%d, ",confRef->img[i].atomName, num_occ(confRef->img, i, confRef->nbAtom),confRef->img[j].atomName, num_occ(confRef->img, j, confRef->nbAtom));
					}
					if(head->CB[i][j]<=0 &&  confRef->CB[i][j]>0){ //(bond broken)
						fprintf(outputC,"-%s%d-%s%d, ",confRef->img[i].atomName, num_occ(confRef->img, i, confRef->nbAtom),confRef->img[j].atomName, num_occ(confRef->img, j, confRef->nbAtom));					
					}
				}
			}	
		}

		//Hydrogen bonds
		if(bit_1(molecule->level,0)){
			for (i = 0; i < confRef->nbAtom; i++){
				for (j = i+1; j < confRef->nbAtom; j++){
					if(confRef->HB[i][j]!=head->HB[i][j]) {
						if(confRef->HB[i][j]==-1 && head->HB[i][j]==-5 && confRef->img[i].atomType!='H' &&  confRef->img[j].atomType!='H'){
							fprintf(outputC,"+%s%d-%s%d, ",confRef->img[i].atomName, num_occ(confRef->img, i, confRef->nbAtom),confRef->img[j].atomName, num_occ(confRef->img, j, confRef->nbAtom));
						}
						if(confRef->HB[i][j]==-5 && head->HB[i][j]==-1){
							fprintf(outputC,"-%s%d-%s%d, ",confRef->img[i].atomName, num_occ(confRef->img, i, confRef->nbAtom),confRef->img[j].atomName, num_occ(confRef->img, j, confRef->nbAtom));
						}
					}			
				}
			}
	
		}

		//Intermolecular interactions
		if(bit_1(molecule->level,2)){	
			for(i=0; i<confRef->nbIon; i++){
				for(j=0;j<confRef->nbAtom;j++){
					if(confRef->IB[i][j]!=head->IB[i][j] && (confRef->IB[i][j]>0 || head->IB[i][j]>0) ){ //Repenser à comment presenter l'info !! si on a mm CN mais pas la mm distribution 
						if(confRef->IB[i][j]< head->IB[i][j]) //bond formed 
							fprintf(outputC,"+%s%d-%s%d, ",confRef->img[confRef->IB[i][confRef->nbAtom]].atomName, num_occ(confRef->img, i, confRef->nbAtom),confRef->img[j].atomName, num_occ(confRef->img, j, confRef->nbAtom));
						else	//bond broken 
							fprintf(outputC,"+%s%d-%s%d, ",confRef->img[confRef->IB[i][confRef->nbAtom]].atomName, num_occ(confRef->img, i, confRef->nbAtom),confRef->img[j].atomName, num_occ(confRef->img, j, confRef->nbAtom));
					}							
				}
			}
		}

		//Organomettalic interactions  
		if(bit_1(molecule->level,3)){
			for(i=0; i<confRef->nbMetal; i++){
				for(j=0;j<confRef->nbAtom;j++){
					if(confRef->MB[i][j]!=head->MB[i][j] && (confRef->MB[i][j]>0 || head->MB[i][j]>0) ){ //Repenser à comment presenter l'info !! si on a mm CN mais pas la mm distribution 
						if(confRef->MB[i][j]< head->MB[i][j]) //bond formed 
							fprintf(outputC,"+%s%d-%s%d, ",confRef->img[confRef->MB[i][confRef->nbAtom]].atomName, num_occ(confRef->img, i, confRef->nbAtom),confRef->img[j].atomName, num_occ(confRef->img, j, confRef->nbAtom));
						else	//bond broken 
							fprintf(outputC,"+%s%d-%s%d, ",confRef->img[confRef->MB[i][confRef->nbAtom]].atomName, num_occ(confRef->img, i, confRef->nbAtom),confRef->img[j].atomName, num_occ(confRef->img, j, confRef->nbAtom));
					}							
				}
			}
		}


		fprintf(outputC,"\"]");
		fprintf(outputC,"\n");

		head = head->suiv;
	}

	fprintf(outputC, "}\n");
	fclose(outputC);

	//Draw the graph  using graphViz
	// char resFile[512];
	// strcpy(resFile,"dot -Tpng  ");
	// strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_comp_to_ref.gv -o \0"));
	// strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_comp_to_ref.png \0"));
	// system(resFile);


}

/*========
Get the coordinate number of the atom index, we sum up all the interactions formed with this atom  
========*/
int get_coord_num(struct confList *conf, int index){
	int j=0;
	int index2 = conf->MB[index][conf->nbAtom];
	int cpt = conf->MB[index][conf->nbAtom+1];
	//Research by lines
	while (j<index2){
		if(conf->CB[j][index2]>0){
			cpt++;
			j++;
		}
		else
			j++;
	}
	if(j==index2) j++;
	//Research by columns
	while (j<conf->nbAtom){
		if(conf->CB[index2][j]>0){
			cpt++;
			j++;
		}
		else
			j++;
	}
	return cpt;

}
/*========
Return string which represents the label of the conformer confName in the graph of comparaison with the reference
========*/
char *get_label_conf(int confName,int size, char *fileName){
	FILE* inputF =fopen(fileName,"r");
	if(inputF== NULL){ printf("Can\'not open file of energies %s\n",fileName); exit(4);}

	char *label = malloc (sizeof (*label) * (5000));
	sprintf(label,"%s","");
	int conf;
	int step=0;
	// int cpt=0;
	double eng;
	while(step<size){
		//Read the conformation
		if(fscanf(inputF,"%d", &conf)!=1){
		printf("can not read the name of the conformation ,%d %s\n",step,fileName);
		exit(-1);
		}   
		//Read the step
		if(fscanf(inputF,"%d", &step)!=1){
		printf("can not read the name of the conformation\n");
		exit(-1);
		}   	
		if(conf==confName){
			//Read the energy value
			if(fscanf(inputF,"%lf", &eng)!=1){
			printf("can not read the name of the conformation\n");
			exit(-1);
			}   
			char a[15];
			sprintf(a,"%d,%.2lf;",step,eng);
			// cpt +=strlen(a);
			// printf("%d : %d : %d \n",confName,strlen(a),cpt);

			// itoa(step,a,10);
			strcat(label,a);
			// printf("%d : %d \n",confName,strlen(label));

		}	
		skip_return(inputF,1);

	}
	fclose(inputF);
	return label;
}

// void save_real_energy(struct MyModel* molecule){
// 	printf("Get energy\n");
// 	char inputF[512]; 
// 	sprintf(inputF,"%s/%s",molecule->inputDir,molecule->inputFN);
// 	FILE* outputC =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_energy.txt\0"),"w");
// 	if(outputC== NULL){ printf("Can\'not open file _energy.txt\n"); exit(4);}
// 	FILE* outputMinMax =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_minmax.txt\0"),"w");
// 	if(outputMinMax== NULL){ printf("Can\'not open file _minmax.txt\n"); exit(4);}
// 	struct confList* head= molecule->conformations;
// 	while(head){
// 		int i;
// 		int nbSnap=0; //Number of snapshots for one conformation 
// 	 	double eng=0; //enery of snapshot
// 	 	int midSnap=0;
// 	 	double minEng=10000; //Sum(eng) 
// 	 	double maxEng=-10000; //Sum(eng^2)
// 		int minSnap;
// 		int maxSnap;
// 		// printf("conformation=%d\n",head->name );	
// 		FILE *traj;
// 		for(i=0;i<head->imgList[0];i++){
// 			if((traj=fopen(inputF,"r"))== NULL){ printf("Can\'not read file .xyz\n"); exit(3);}
// 			//Go to the snapshot head->imgList[i*2+1]
// 			skip_snapshot(traj, head->imgList[i*2+1], snapshot_size(traj)); 	
// 			int j=head->imgList[i*2+1];
// 			if(nbSnap<((head->imgList[i*2+2]-head->imgList[i*2+1])/2)){
// 				nbSnap = (head->imgList[i*2+2]-head->imgList[i*2+1])/2;
// 				midSnap = (head->imgList[i*2+2]+head->imgList[i*2+1])/2;
// 			} 
// 			while(j<= head->imgList[i*2+2]){
// 				eng= get_snap_energy(traj,molecule->nbAtom);
// 				fprintf(outputC, "%2d \t %2d \t", head->name, j);//conformation name
// 				fprintf(outputC, "%lf ",eng);//Energy average ,  eng*630.00
// 				fprintf(outputC, "\n");
// 				minEng =Min(minEng,eng);
// 				if(minEng==eng) minSnap = j;
// 				maxEng =Max(maxEng,eng);
// 				if(maxEng==eng) maxSnap = j;
// 				// nbSnap++;	
// 				j++;			
// 			}
// 			fclose(traj);
// 		}
// 		//Save  results 
// 		// double A=S1/nbSnap; //Average of energy for this conformation
// 		// fprintf(outputC, "%2d \t", head->name);//conformation name
// 		// fprintf(outputC, "%lf \t", A);//Energy average
// 		// fprintf(outputC, "%lf \n", sqrt((S2/nbSnap)-A*A)));//nergy standard deviation
// 		//Go to the next conformation
// 		fprintf(outputMinMax, "%s \t %2d \t %2d \t %lf \n","min",head->name, minSnap, minEng );
// 		fprintf(outputMinMax, "%s \t %2d \t %2d \t %lf \n","max",head->name, maxSnap, maxEng );
// 		fprintf(outputMinMax, "%s \t %2d \t %2d \n","mid",head->name, midSnap );
// 		head=head->suiv;
// 	}
// 	fclose(outputC);
// 	fclose(outputMinMax);
// }

/*========
Plot the time evolution of H-bonds identified along the analysed trajectory, using Gnuplot.
========*/
void plot_Hbonds(struct MyModel *molecule){
	FILE* plotF=NULL;
	char resFile[256];
	char rad[256];
	int i;
//Create the ".gplt" file	
	strcpy(rad,molecule->inputFN);
	rad[strlen(molecule->inputFN)-strlen(strrchr(molecule->inputFN,'.'))]='\0';
	plotF =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_Hbonds.gplt\0"),"w");

	fprintf(plotF, "set term 'png'\n");
	fprintf(plotF, "set output '%s/%s/%s%s'\n",molecule->inputDir,rad,rad,"_HbondsFin.png");
	fprintf(plotF, "set title 'Molecular dynamics: %s'\n",rad);
	fprintf(plotF, "set xlabel '#snapshot'\n");
	fprintf(plotF, "set ylabel 'H-bonds'\n");
	fprintf(plotF, "set xrange [0:%d]\n",molecule->nbMolecule);
	fprintf(plotF, "set yrange [0:%d]\n",molecule->nbHbond+1);
	fprintf(plotF, "unset ytics\n");
	fprintf(plotF, "set grid\n");
	for(i=0; i<molecule->nbHbond;i++)
		fprintf(plotF, "set label '%s%d-%s%d-%s%d' at 10,%d.2 left\n",molecule->imgRef[molecule->bondHdyn[i].bondH[0]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[0], molecule->nbAtom),molecule->imgRef[molecule->bondHdyn[i].bondH[1]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[1], molecule->nbAtom),molecule->imgRef[molecule->bondHdyn[i].bondH[2]].atomName, num_occ(molecule->imgRef, molecule->bondHdyn[i].bondH[2], molecule->nbAtom),i+1);
	fprintf(plotF, "plot ");
	fprintf(plotF, "'%s/%s/%s.Hbonds' using 1:2:3 with lines notitle lw 3 lc variable,  \x5C\n",molecule->inputDir,rad,rad);
	fclose(plotF);

//Generate the ".png" file 
	sprintf(resFile,"gnuplot %s",get_FN(molecule->inputDir,molecule->inputFN,"_Hbonds.gplt\0"));
	system(resFile);

}

/*========
Plot the time evolution of the identified conformations, using Gnuplot
========*/
void plot_conf_periods(struct MyModel *molecule){
	printf("plot periods\n");
	FILE* plotF=NULL;
	char resFile[256];
	char rad[256];
//Create the ".gplt" file	
	strcpy(rad,molecule->inputFN);
	rad[strlen(molecule->inputFN)-strlen(strrchr(molecule->inputFN,'.'))]='\0';
	plotF =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_perd.gplt\0"),"w");

	fprintf(plotF, "set term 'png'\n");
	fprintf(plotF, "set output '%s/%s/%s%s'\n",molecule->inputDir,rad,rad,"_perdFin.png");
	fprintf(plotF, "set title 'The conformational dynamics of %s'\n",rad);
	fprintf(plotF, "set xlabel '#snapshot'\n");
	fprintf(plotF, "set ylabel 'Conformations'\n");
	fprintf(plotF, "set xrange [0:%d]\n",molecule->nbMolecule);
	fprintf(plotF, "set yrange [0:%d]\n",molecule->nbConf+1);
	fprintf(plotF, "set ytics 0,1,%d\n",molecule->nbConf+1);
	fprintf(plotF, "set grid\n");

	fprintf(plotF, "plot ");
	fprintf(plotF, "'%s/%s/%s_perd.data' using 1:2:3 with lines notitle lw 5 lc variable,  \x5C\n",molecule->inputDir,rad,rad);
	fclose(plotF);

//Generate the ".png" file 
	sprintf(resFile,"gnuplot %s",get_FN(molecule->inputDir,molecule->inputFN,"_perd.gplt\0"));
	system(resFile);
}


/*========
Draw the mixed graphs for each conformation identified along the trajectory (ies), using GraphViz
========*/
void draw_mixed_graph(struct MyModel *molecule, char type){
//create the graph for every conformation using graphViz
	struct confList* head= molecule->conformations;
    	char graphFN[270]="";
	FILE* outputG=NULL; 
	//Save adjacency matrix if a single analysis 
	FILE* outputM=NULL;
	if(molecule->nbInputF==1) outputM= fopen(get_FN(molecule->inputDir,molecule->inputFN,"_confMatrix.txt\0"),"w");
	
	while(head){	
		head->state='C';
		int **adjMatrix= allocate_matrix(head->nbAtom,head->nbAtom,0,"adjMatrix");
		if(molecule->nbInputF==1)fprintf(outputM, "conformer %d : %d x %d \n",head->name, head->nbAtom, head->nbAtom);
		//Get the graph		
		switch(type){
			case 's':{ //Single trajectory
				outputG =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_conf.gv\0"),"w");
				//Save the cartesian coordinates :
				save_snapshot(molecule, head);
				save_conf_bonds(molecule, head,'s');
			}
			break;
			case 'c':{ //Multiple trajectories - conformations
				sprintf(graphFN,"%s/%s",molecule->inputDir,"confGraph.gv");
				outputG=fopen(graphFN,"w");
			}
			break;
			case 'i':{ //Multiple trajectories - isomers
				sprintf(graphFN,"%s/%s",molecule->inputDir,"isomGraph.gv");
				outputG=fopen(graphFN,"w");
			}
			break;

			default : printf("ERROR on output type\n");
		}
		fprintf(outputG, "digraph G {\n");
		// fprintf(outputG, "label=\"Conformation %d (%.3lf) \";\n",head->name, head->totalPerc);
		fprintf(outputG, "label=\"conf.%d \";\n",head->name);
		fprintf(outputG, "node [style=filled];\n");
		fprintf(outputG, "graph [bgcolor=transparent];\n");
		fprintf(outputG, "node [shape = circle, fontsize=12];\n");

		int i=0;
		int j=0;
		// printf("%c %s\n",type,molecule->sysType );
		//Draw atoms except hydrogen atoms
		for(i=0;i<head->nbAtom;i++){
			if(head->img[i].atomIn){
				//Oxygen in red
				if(strcmp(head->img[i].atomName,"O")==0){
					if(molecule->sysType[0]=='w'){
						if(get_nb_CB_type(head->img,i,head->CB,head->nbAtom,"H")==2)
							fprintf(outputG, "\"%s%d\"[fillcolor=red, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
						else 
							fprintf(outputG, "\"%s%d\"[fillcolor=goldenrod1, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
					}
					else
						fprintf(outputG, "\"%s%d\"[fillcolor=red, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));

				} 
				//Neitrogen in blue
				if(strcmp(head->img[i].atomName,"N")==0)
					fprintf(outputG, "\"%s%d\"[fillcolor=blue, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
				//Carbon in gray
				if(strcmp(head->img[i].atomName,"C")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=gray32, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
				//Lithium in blueviolet
				if(strcmp(head->img[i].atomName,"Li")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=blueviolet, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
				//Chlore in green
				if(strcmp(head->img[i].atomName,"Cl")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=green, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
				//Sulphur in cyan
				if(strcmp(head->img[i].atomName,"S")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=cyan, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
				//Sulphur in cornsilk3
				if(strcmp(head->img[i].atomName,"Si")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=cornsilk3, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
				//Manganese in deeppink1 
				if(strcmp(head->img[i].atomName,"Mn")==0){
					fprintf(outputG, "\"%s%d\"[fillcolor=deeppink1, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				} 
				//Ruthinium in deeppink3 
				if(strcmp(head->img[i].atomName,"Ru")==0){
					fprintf(outputG, "\"%s%d\"[fillcolor=deeppink3, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				} 
				//Potassium in midnightblue  
				if(strcmp(head->img[i].atomName,"K")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=midnightblue, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Sodium in  steelblue 
				if(strcmp(head->img[i].atomName,"Na")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=steelblue, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Boron in lightpink  
				if(strcmp(head->img[i].atomName,"B")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=lightpink, fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Brom in lightgreen  
				if(strcmp(head->img[i].atomName,"Br")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=seagreen1 , fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Phosphorus in dodgerblue1  
				if(strcmp(head->img[i].atomName,"P")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=dodgerblue1 , fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Flourine in darkgoldenrod1  
				if(strcmp(head->img[i].atomName,"F")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=darkgoldenrod1 , fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Iodine in chocolate1  
				if(strcmp(head->img[i].atomName,"I")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=chocolate1 , fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Gold in chocolate1  
				if(strcmp(head->img[i].atomName,"Au")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=chocolate1 , fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Zinc in azure4  
				if(strcmp(head->img[i].atomName,"Zn")==0) 
					fprintf(outputG, "\"%s%d\"[fillcolor=azure4 , fontcolor=white, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				//Hydrogen in gray
				if(type=='c' && molecule->sysType[0]=='f' && strcmp(head->img[i].atomName,"H")==0){
					fprintf(outputG, "\"%s%d\"[fillcolor=gray87 , fontcolor=black, fontname=\"blod, bold\"];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));					
				} 
			}
		}
		//Make the graph according the covalent bonds
		for(i=0;i<head->nbAtom;i++){
			if(head->img[i].atomIn){
				for(j=i+1;j<head->nbAtom;j++){
					if(head->CB[i][j]>0 ){ //covalent fixed or not and not leave
						adjMatrix[i][j]=1;
						adjMatrix[j][i]=1;
						if(head->img[i].atomType!='H' && head->img[j].atomType!='H'){
							int isthm = verif_exist_isthm(molecule->isthList,i,j);
							if(isthm==-1){ //simple covalent bond
								// if(head->CB[i][j]==1)
									fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );
								// else
								// 	fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = \"black:invis:black\" , nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );							
							}   
							else{//rotational axis 
								if(isthm==2){ //conformational rotation
									// if(head->CB[i][j]==1)
										fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color=green, nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );
									// else
									// 	fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color=\"green:invis:green\", nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );								
								}
								else{	
									// if(head->CB[i][j]==1)	
										fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color=firebrick, nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );
									// else
									// 	fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color=\"firebrick:invis:firebrick\", nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );								
								}				
							} 
						}
						else{//make change
							if(type=='i'){
								fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );					
							} //In isomers we present hydrogen
							if(head->img[i].atomType=='H' && (head->CB[i][i]==11 || (type=='c' && molecule->sysType[0]=='f'))) {
								fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );
							}
							if(head->img[j].atomType=='H' && (head->CB[j][j]==11 || (type=='c' && molecule->sysType[0]=='f'))){
								fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );
							}
						}
					}
				}
			}
		}
				
		//Make the graph according the rotational motion
		// for(i=0;i<molecule->isthList[0];i++){
		// 	if(molecule->isthList[i*3+3]==2)
		// 		fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color=forestgreen, nodesep=0.5];\n",head->img[molecule->isthList[i*3+1]].atomName, num_occ(head->img, molecule->isthList[i*3+1], head->nbAtom),head->img[molecule->isthList[i*3+2]].atomName, num_occ(head->img, molecule->isthList[i*3+2], head->nbAtom) );
		// 	else
		// 		fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color=firebrick, nodesep=0.5];\n",head->img[molecule->isthList[i*3+1]].atomName, num_occ(head->img, molecule->isthList[i*3+1], head->nbAtom),head->img[molecule->isthList[i*3+2]].atomName, num_occ(head->img, molecule->isthList[i*3+2], head->nbAtom) );
		// }
		//Make the graph according the H-bonds
		int h=1;
		if(bit_1(molecule->level,0)){
			for(i=0;i<head->nbAtom;i++){
				if(head->img[i].atomIn){
					for(j=i+1;j<head->nbAtom;j++){ //Update
						if(head->HB[i][j]==-5 && head->img[i].atomType!='H' && head->img[j].atomType!='H'){
							// printf("i=%s%d, j=%s%d \n", head->img[i].atomName, num_occ(head->img, i, head->nbAtom), head->img[j].atomName, num_occ(head->img, j, head->nbAtom));	
							adjMatrix[i][j]=4;
							adjMatrix[j][i]=4;
							//check if there is proton transfer or not
							int D=-1;
							int A=-1;
							bool ptr=check_ptr(head,i,j,head->nbAtom,&D);							
								if(D==i)
									A=j;
								if(D==j)
									A=i;
								if(D==-1){
									printf("Error on the donor (3)\n");
									exit(10);
								}	
							if(ptr){//proton transfer A-H-D
								fprintf(outputG, "\"%s%d\"[fillcolor=gray87];\n","H",h);
								fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[A].atomName, num_occ(head->img, A, head->nbAtom),"H",h );					
								fprintf(outputG, "\"%s%d\"->\"%s%d\"[fontcolor=blue, color=blue , style=dashed];\n","H",h,head->img[D].atomName, num_occ(head->img, D, head->nbAtom));
								h++;								
							}
							else{//D-H-A
								//printf("%d) %s%d -- %s%d \n",h, head->img[D].atomName, num_occ(head->img, D, head->nbAtom),head->img[A].atomName, num_occ(head->img, A, head->nbAtom) ); 
								//int H= get_donHB(molecule->donHacc,D,head->nbAtom);
								//fprintf(outputG, "\"%s%d\"[fillcolor=gray87];\n","H",h);
								//fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[D].atomName, num_occ(head->img, D, head->nbAtom),"H",h );					
								//fprintf(outputG, "\"%s%d\"->\"%s%d\"[fontcolor=red, color=red , style=dashed];\n","H",h,head->img[A].atomName, num_occ(head->img, A, head->nbAtom));							
								fprintf(outputG, "\"%s%d\"->\"%s%d\"[fontcolor=red, color=red , style=dashed];\n",head->img[D].atomName, num_occ(head->img, D, head->nbAtom),head->img[A].atomName, num_occ(head->img, A, head->nbAtom));							
								h++;								
							}
						}
					}
				}
				// else{
				// 	printf("not %d \n",i );
				// }
			}	
		}	
		//Make the graph according the H-bonds (2nd way with Hbond)
		// for(i=0;i<head->nbAtom;i++)
		// 	for(j=0;j<head->nbAtom;j++){
		// 		if(head->HB[i][j]==1 || head->HB[i][j]==0){ //Simple H-bond
		// 			int h=0;
		// 			while(h<head->nbAtom){	
		// 				//if(get_Hbond_sens(head,head->nbAtom,j,i)==2){
		// 				if(head->HB[j][h]==2 && (head->CB[i][h]==1 || head->CB[h][i]==1)){ //D-H-A 
		// 					fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[h].atomName, num_occ(head->img, h, head->nbAtom) );							
		// 					fprintf(outputG, "\"%s%d\"->\"%s%d\"[fontcolor=red, color=red , style=dashed];\n",head->img[h].atomName, num_occ(head->img, h, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
		// 				}
		// 				//else{ 
		// 				if(head->HB[j][h]==3){ //A-H-D head->HB[j][h]==3
		// 					fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[j].atomName, num_occ(head->img, j, head->nbAtom),head->img[h].atomName, num_occ(head->img, h, head->nbAtom) );		
		// 					fprintf(outputG, "\"%s%d\"->\"%s%d\"[fontcolor=blue, color=blue , style=dashed];\n",head->img[h].atomName, num_occ(head->img, h, head->nbAtom),head->img[i].atomName, num_occ(head->img, i, head->nbAtom));												
		// 				}
		// 				h++;
		// 			}
		// 		}
		// 	}
		//Make the graph according the Intermolecular bonds
		if(bit_1(molecule->level,2)){
			for(i=0;i<head->nbIon;i++){
				if(head->IB[i][head->nbAtom]!=-1 && head->IB[i][head->nbAtom+1]>0 ){
					for(j=0;j<head->nbAtom;j++)
						if(head->IB[i][j]>0 && !((strcmp(head->img[j].atomName,"Ar")==0 || strcmp(head->img[j].atomName,"K")==0 || strcmp(head->img[j].atomName,"Li")==0 ||  strcmp(head->img[j].atomName,"Na")==0 || strcmp(head->img[j].atomName,"F")==0 || strcmp(head->img[j].atomName,"I")==0 || strcmp(head->img[j].atomName,"Cl")==0 || strcmp(head->img[j].atomName,"Br")==0) && get_index_ion(head->IB,head->nbIon,head->nbAtom,j) < i)){ //Intermolecular bond
							adjMatrix[head->IB[i][head->nbAtom]][j]=6;
							adjMatrix[j][head->IB[i][head->nbAtom]]=6;
							fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color=blue, nodesep=0.5,style=dashed];\n",head->img[head->IB[i][head->nbAtom]].atomName, num_occ(head->img, head->IB[i][head->nbAtom], head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );			
							if(strcmp(head->img[j].atomName,"H")==0 ){
								int d=get_cov_index(head->CB,head->nbAtom,j);
								if(d!=-1)
									fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[d].atomName, num_occ(head->img, d, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );
								else{
									printf("Error in cov index with hydrogen %d\n",j); exit(-1);
								}
							}
						}		
				}
			}
		}
					
		//Make the graph according the  organometallic bonds
		if(bit_1(molecule->level,3)){
			for(i=0;i<head->nbMetal;i++){
				if(head->MB[i][head->nbAtom]!=-1 && head->MB[i][head->nbAtom+1]>0 ){
					for(j=0;j<head->nbAtom;j++)
						if(head->MB[i][j]>0 && !((strcmp(head->img[j].atomName,"Mn")==0 ||  strcmp(head->img[j].atomName,"Ru")==0 ||  strcmp(head->img[j].atomName,"Au")==0) && get_index_ion(head->MB,head->nbMetal,head->nbAtom,j) < i)){ //Intermolecular bond
							adjMatrix[head->MB[i][head->nbAtom]][j]=7;
							adjMatrix[j][head->MB[i][head->nbAtom]]=7;
							fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color=deeppink1, nodesep=0.5,style=dashed];\n",head->img[head->MB[i][head->nbAtom]].atomName, num_occ(head->img, head->MB[i][head->nbAtom], head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );			
							if(strcmp(head->img[j].atomName,"H")==0 ){
								int d=get_cov_index(head->CB,head->nbAtom,j);
								if(d!=-1)
									fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black, nodesep=0.5];\n",head->img[d].atomName, num_occ(head->img, d, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom) );
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
		//Draw the graph of transition using graphViz
		char resFile[512]="";
		strcpy(resFile,"neato -Tpng  ");
		switch(type){
			case 's':{ //Single trajectory
				strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_conf.gv -o \0"));
				if(head->state=='C' || molecule->nbMolecule==1)	
					strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_conf/conformation\0"));
				else
					strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_trans/tranState\0"));
			}
			break;
			case 'c':{ //Multiple trajectories - conformations
				strcat(resFile,molecule->inputDir);
				strcat(resFile,"/confGraph.gv -o ");
				strcat(resFile,molecule->inputDir);
				strcat(resFile,"/conf_graph/conformation");
			}
			break;
			case 'i':{ //Multiple trajectories - isomers
				strcat(resFile,molecule->inputDir);
				strcat(resFile,"/isomGraph.gv -o ");
				strcat(resFile,molecule->inputDir);
				strcat(resFile,"/isom_graph/isomer");
			}
			break;

			default : printf("ERROR on output type\n");
		}		

		char a[5];
		sprintf(a,"%d",head->name);
		strcat(resFile,a); 
		strcat(resFile,".png\0");	
		system(resFile);
		//Save the adjacency matrix 
		if(molecule->nbInputF==1){
			for(i=0;i<head->nbAtom;i++){fprintf(outputM," %s%d",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));}
			fprintf(outputM,"\n");
			for(i=0;i<head->nbAtom;i++){
				fprintf(outputM,"%s%d ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom));
				for(j=0;j<head->nbAtom;j++){
					fprintf(outputM, "%d ",adjMatrix[i][j]);
				}
				fprintf(outputM,"\n");
			}
			fprintf(outputM, "\n\n");
			free(adjMatrix);
		}
		head = head->suiv;
	}
	if(molecule->nbInputF==1)fclose(outputM);
	char resFilegv[512]="";
	sprintf(resFilegv,"rm %s",get_FN(molecule->inputDir,molecule->inputFN,"_conf.gv\0"));
	system(resFilegv);

} 

/*========
Draw the graph of transitions using GraphViz
========*/
void draw_trans_graph(struct MyModel *molecule){
//Save the graph of transitions between conformations
	FILE* outputG=NULL;
	outputG =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_graphTrans.gv\0"),"w");
//Save the Adjacency matrix of transitions 	
	FILE* outputM=NULL;
	outputM =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_TransMatrix.txt\0"),"w");
	int **transMatrix = allocate_matrix(molecule->nbConf,molecule->nbConf,0,"transMatrix"); //[molecule->nbConf][molecule->nbConf];
//Graph with frequencies
	fprintf(outputG, "digraph G {\n");
	fprintf(outputG, "label=\"Graph of transitions\";\n");
	fprintf(outputG, "node [style=filled,  fontcolor=black, fontname=\"blod, bold\"];\n");

	struct confList* head= molecule->conformations;
	int i;
	int j;
	while(head){
		head->state='C';
	if(head->state=='C'){
		fprintf(outputG, "%d ", head->name);
		if(head->state=='C'){
			fprintf(outputG, "[fillcolor=indianred1, ");
		}
		else{ //Intermediate 
			fprintf(outputG, "[fillcolor=darkseagreen2, ");
		}
		if(head->name==molecule->lastConf->name){
			fprintf(outputG, "shape=box, ");
		}
		//Fill the bonds
		fprintf(outputG," label=\" %d :\n",head->name);
		//Covalent bonds
		// fprintf(outputG, "C : "); 
		if(bit_1(molecule->changType,1) && !(molecule->partAnal)){
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==1){
						if(head->CB[i][j]>0) 
							fprintf(outputG,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
					}
				}
			}
			fprintf(outputG, ".\n");	
		}
		
		//Hydrogen bonds
		// fprintf(outputG, "H : "); 
		if(bit_1(molecule->changType,0)){
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==2){
						if(head->HB[i][j]==-5){
							int D=-1;
							int A=-1;
							bool ptr=check_ptr(head,i,j,head->nbAtom,&D);							
							if(D==i)
								A=j;
							if(D==j)
								A=i;
							if(D==-1){
								printf("Error on the donor (3)\n");
								exit(10);
							}	
							if(ptr){//proton transfer A-H-D
								fprintf(outputG,"%s%d-%s%d, ",head->img[A].atomName, num_occ(head->img, A, head->nbAtom),head->img[D].atomName, num_occ(head->img, D, head->nbAtom));								
							}
							else{
								fprintf(outputG,"%s%d-%s%d, ",head->img[D].atomName, num_occ(head->img, D, head->nbAtom),head->img[A].atomName, num_occ(head->img, A, head->nbAtom));								
							}
						}	

					}
				}
			}
			fprintf(outputG, ".\n");
		}

		//Intermolecular interaction
		// fprintf(outputG, "I : "); 
		if(bit_1(molecule->changType,2)){
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==3){
						int I =-1; 
						int A=-1;
						if(strcmp(head->img[i].atomName,"Ar")==0 || strcmp(head->img[i].atomName,"K")==0 || strcmp(head->img[i].atomName,"Li")==0 || strcmp(head->img[i].atomName,"Na")==0 || strcmp(head->img[i].atomName,"F")==0 || strcmp(head->img[i].atomName,"I")==0 || strcmp(head->img[i].atomName,"Cl")==0){
							I=get_index_ion(head->IB,head->nbIon,head->nbAtom,i);
							
							A=j;
						}
						else{
							I=get_index_ion(head->IB,head->nbIon,head->nbAtom,j);
							A=i;
						}						
						if(I !=-1 && head->IB[I][A]>0){
							fprintf(outputG,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
						}
					}
				}
			}
			fprintf(outputG, ".\n");
		}
		//Organometallic bonds 
		// fprintf(outputG, "O : "); 
		if(bit_1(molecule->changType,3)){
			for(i=0;i<molecule->nbAtomDyn;i++){
				for(j=i+1;j<molecule->nbAtomDyn;j++){
					if(molecule->bondDyn[i][j]==4){
						int I =-1; 
						int A=-1;
						if(strcmp(head->img[i].atomName,"Mn")==0 || strcmp(head->img[i].atomName,"Ru")==0 ||  strcmp(head->img[i].atomName,"Au")==0){
							I=get_index_ion(head->MB,head->nbMetal,head->nbAtom,i);
							A=j;
						}
						else{
							I=get_index_ion(head->MB,head->nbMetal,head->nbAtom,j);
							A=i;
						}
						if(I !=-1 && head->MB[I][A]>0)
							fprintf(outputG,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
					}
				}
			}
			fprintf(outputG, ".\n");
		}
		fprintf(outputG, "%.2lf %% \n \"] \n",head->totalPerc);
	}
		head =head->suiv;
	}
	// Get the type of transitions
	head= molecule->conformations;
	while(head){
		struct successor *newSucc=head->succ;
		while(newSucc){
		if(newSucc->conf->state=='C' && head->state=='C'){
			char bondChg[256]="";
			char bond[15]="";
			transMatrix[head->name-1][newSucc->conf->name-1]=1;
			//Get type of change newSucc->transf
			fprintf(outputG, "%d -> %d [color=black, label=\" \n", head->name, newSucc->conf->name);
			//Bonds diff
			//Covalent bonds change
			if(bit_1(newSucc->transf,3) || bit_1(newSucc->transf,4)) {
				strcat(bondChg,",CB =\"\0");
				fprintf(outputG,"CB =\'");
				for(i=0;i<molecule->nbAtomDyn;i++){
					for(j=i+1;j<molecule->nbAtomDyn;j++){
						// if(head->name==40 && newSucc->conf->name == 41 && molecule->bondDyn[i][j]!=-1 ){
						// // if(head->name==29 && newSucc->conf->name == 30 && i==32 && j==36){
						// 	printf("xx %d : %d : %d : %s : %s : %d : %d : %d : %d \n",i,j,molecule->bondDyn[i][j], head->img[i].atomName, head->img[j].atomName,head->CB[i][j],newSucc->conf->CB[i][j],head->CB[j][i],newSucc->conf->CB[j][i]);
						// 	printf("%s", head->img[i].atomIn ? "Hi : true \n" : "Hi :false\n");
						// 	printf("%s", head->img[j].atomIn ? "Hj : true \n" : "Hj :false\n");
						// 	printf("%s", newSucc->conf->img[i].atomIn ? "Si : true \n" : "Si :false\n");
						// 	printf("%s", newSucc->conf->img[j].atomIn ? "Sj : true \n" : "Sj :false\n");
						// }
						if( (molecule->bondDyn[i][j]==1 || (molecule->bondDyn[i][j]==4 && head->CB[i][j]!= newSucc->conf->CB[i][j])) && ((head->img[i].atomIn || head->img[j].atomIn) &&  (newSucc->conf->img[i].atomIn || newSucc->conf->img[j].atomIn) ) ){//!(molecule->partAnal) //
						 	if(head->CB[i][j]<newSucc->conf->CB[i][j] && head->CB[i][j]!=-1 && newSucc->conf->CB[i][j]!=-1/*>0*/ ){
								fprintf(outputG,"+%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								sprintf(bond,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));							
								// sprintf(bond,"+%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));							
								strcat(bondChg, bond);														 		
						 	} 
						 	if(head->CB[i][j]>newSucc->conf->CB[i][j] && head->CB[i][j]!=-1 && newSucc->conf->CB[i][j]!=-1/*>0 && ==0*/ ){
								fprintf(outputG,"-%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								sprintf(bond,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));							
								// sprintf(bond,"-%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));							
								strcat(bondChg, bond);														 								 		
						 	} 
						}
					}
				}
				fprintf(outputG, ".\'\n");
				strcat(bondChg,"\"\0");				
			}
			//Hydrogen bonds change
			if(bit_1(newSucc->transf,0) ||bit_1(newSucc->transf,1) || bit_1(newSucc->transf,2) ) {
				strcat(bondChg,",HB =\"\0");
				fprintf(outputG,"HB =\'");				
				if( bit_1(newSucc->transf,2))
					fprintf(outputG, "(T), ");								
				for(i=0;i<molecule->nbAtomDyn;i++){
					for(j=i+1;j<molecule->nbAtomDyn;j++){
						if(molecule->bondDyn[i][j]==2){
						 	if(newSucc->conf->HB[i][j]!= head->HB[i][j] &&   newSucc->conf->HB[i][j]==-5){
								fprintf(outputG,"+%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								sprintf(bond,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								// sprintf(bond,"+%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								strcat(bondChg, bond);
						 	}
						 	if(head->HB[i][j]==-5 && newSucc->conf->HB[i][j]!= head->HB[i][j] ){
								fprintf(outputG,"-%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								sprintf(bond,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								// sprintf(bond,"-%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								strcat(bondChg, bond);
						 	} 
						}
					}
				}
				fprintf(outputG, ".\'\n");
				strcat(bondChg,"\"\0");				
			}
			//Intermolecular interactions change
			if(bit_1(newSucc->transf,5) || bit_1(newSucc->transf,6)) {
				strcat(bondChg,",IB =\"\0");
				fprintf(outputG,"IB =\'");				
				for(i=0;i<molecule->nbAtomDyn;i++){
					for(j=i+1;j<molecule->nbAtomDyn;j++){
						if(molecule->bondDyn[i][j]==3){
							int I =-1; 
							int A=-1;
							if(strcmp(head->img[i].atomName,"Ar")==0 || strcmp(head->img[i].atomName,"K")==0 || strcmp(head->img[i].atomName,"Li")==0 || strcmp(head->img[i].atomName,"Na")==0 || strcmp(head->img[i].atomName,"F")==0 || strcmp(head->img[i].atomName,"I")==0 || strcmp(head->img[i].atomName,"Cl")==0 || strcmp(head->img[i].atomName,"Br")==0){
								I=get_index_ion(head->IB,head->nbIon,head->nbAtom,i);
								A=j;
							}
							else{
								I=get_index_ion(head->IB,head->nbIon,head->nbAtom,j);
								A=i;
							}						
							if(I !=-1 && head->IB[I][A]>0 && newSucc->conf->IB[I][A]<=0 ){
								fprintf(outputG,"-%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								sprintf(bond,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								// sprintf(bond,"+%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								strcat(bondChg, bond);
						}
							if(I !=-1 && head->IB[I][A]<=0 && newSucc->conf->IB[I][A]>0){
								fprintf(outputG,"+%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								sprintf(bond,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								// sprintf(bond,"-%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								strcat(bondChg, bond);
							}
						}
					}
				}
				fprintf(outputG, ".\'\n");
				strcat(bondChg,"\"\0");				
			}
			//Organometallic bonds change
			if(bit_1(newSucc->transf,7) || bit_1(newSucc->transf,8) ) {
				strcat(bondChg,",OB =\"\0");
				fprintf(outputG,"OB =\'");				
				for(i=0;i<molecule->nbAtomDyn;i++){
					for(j=i+1;j<molecule->nbAtomDyn;j++){
						//if(molecule->bondDyn[i][j]==4){
						if( molecule->bondDyn[i][j]==4 || (molecule->bondDyn[i][j]==1 && head->CB[i][j]== newSucc->conf->CB[i][j]) ){
							int I =-1; 
							int A=-1;
							if(strcmp(head->img[i].atomName,"Mn")==0 || strcmp(head->img[i].atomName,"Ru")==0 ||  strcmp(head->img[i].atomName,"Au")==0){
								I=get_index_ion(head->MB,head->nbMetal,head->nbAtom,i);
								A=j;
							}
							else{
								I=get_index_ion(head->MB,head->nbMetal,head->nbAtom,j);
								A=i;
							}						
							if(I !=-1 && head->MB[I][A]>0 && newSucc->conf->MB[I][A]<=0 ){
								fprintf(outputG,"-%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								sprintf(bond,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								// sprintf(bond,"+%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								strcat(bondChg, bond);
							}
							if(I !=-1 && head->MB[I][A]<=0 && newSucc->conf->MB[I][A]>0){
								fprintf(outputG,"+%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								sprintf(bond,"%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								// sprintf(bond,"-%s%d-%s%d, ",head->img[i].atomName, num_occ(head->img, i, head->nbAtom),head->img[j].atomName, num_occ(head->img, j, head->nbAtom));
								strcat(bondChg, bond);
							}
						}
					}
				}
				fprintf(outputG, ".\'\n");
				strcat(bondChg,"\"\0");				
			}
			//Frequency
			fprintf(outputG, "freq= %d\n", newSucc->freq);
			fprintf(outputG, "\" \n%s\n] \n",bondChg);
			// fprintf(outputG, "//%s\n",bondChg);				
		}
			newSucc=newSucc->suiv;
		}
		head = head->suiv;
	}	
	fprintf(outputG, "}\n");
	fclose(outputG);
//Save the adjacency matrix (transitions between conformers)
	fprintf(outputM, "%d x %d \n", molecule->nbConf, molecule->nbConf);
	for(i=0;i<molecule->nbConf;i++){
		for(j=0;j<molecule->nbConf;j++){
			fprintf(outputM, "%d ",transMatrix[i][j]);
		}
		fprintf(outputM, "\n");
	}
	fclose(outputM);
//Draw the graph of transition with frequencies using graphViz
	char resFile[512];
	strcpy(resFile,"dot -Tpng  ");
	strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_graphTrans.gv -o \0"));
	strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_graphTrans.png \0"));
	system(resFile);
	sprintf(resFile,"rm %s",get_FN(molecule->inputDir,molecule->inputFN,"_graphTrans.gv\0"));
	system(resFile);
}

// void draw_trans_graph(struct MyModel *molecule){
// //Save the graph of transitions between conformations
// 	FILE* outputG=NULL;
// 	FILE* outputG2=NULL; 
// 	outputG =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_graphFreq.gv\0"),"w");
// 	outputG2 =fopen(get_FN(molecule->inputDir,molecule->inputFN,"_graphChang.gv\0"),"w");
// //Graph with frequencies
// 	fprintf(outputG, "digraph G {\n");
// 	fprintf(outputG, "label=\"Graph of conformations with frequencies\";\n");
// 	fprintf(outputG, "node [margin=0, fontsize=12, shape=point, style=filled];\n");
// 	fprintf(outputG, "graph [bgcolor=transparent];\n");
// //Graph with nature of changes
// 	fprintf(outputG2, "digraph G {\n");
// 	fprintf(outputG2, "label=\"Graph of conformations with nature of dynamics\";\n");
// 	fprintf(outputG2, "node [margin=0, fontsize=12, shape=point, style=filled];\n");
// 	fprintf(outputG2, "graph [bgcolor=transparent];\n");
// 	//First state
// 	//Add \n ; totalPerc
// 	if(molecule->conformations->state=='C'){
// 		fprintf(outputG, "\"conf-%d\"[shape=circle, fillcolor=lemonchiffon];\n", molecule->conformations->name);//first conformation
// 		fprintf(outputG2, "\"conf-%d\"[shape=circle, fillcolor=lemonchiffon];\n", molecule->conformations->name);//first conformation		
// 	}
// 	else{
// 		fprintf(outputG, "\"inter-%d\"[shape=circle, fillcolor=lemonchiffon];\n", molecule->conformations->name);//first conformation
// 		fprintf(outputG2, "\"inter-%d\"[shape=circle, fillcolor=lemonchiffon];\n", molecule->conformations->name);//first conformation
// 	}	
// 	//last state
// 	if(molecule->lastConf->state=='C'){
// 		fprintf(outputG, "\"conf-%d\"[shape = circle, fillcolor=lightblue1];\n", molecule->lastConf->name);//last conformation
// 		fprintf(outputG2, "\"conf-%d\"[shape = circle, fillcolor=lightblue1];\n", molecule->lastConf->name);//last conformation
		
// 	}
// 	else{
// 		fprintf(outputG, "\"inter-%d\"[shape = circle, fillcolor=lightblue1];\n", molecule->lastConf->name);//last conformation
// 		fprintf(outputG2, "\"inter-%d\"[shape = circle, fillcolor=lightblue1];\n", molecule->lastConf->name);//last conformation
// 	}	

// 	//Stable conformation (state)
// 	struct confList* head= molecule->conformations;
// 	while(head){
// 		if(head->state=='C'){
// 			fprintf(outputG, "\"conf-%d \"[shape = circle, fillcolor=indianred1, fontcolor=black, fontname=\"blod, bold\"];\n",head->name);
// 			fprintf(outputG2, "\"conf-%d\"[shape = circle, fillcolor=indianred1, fontcolor=black, fontname=\"blod, bold\"];\n",head->name);
// 		}
// 		else{
// 			fprintf(outputG, "\"inter-%d\"[shape = circle, fillcolor=darkseagreen2, fontcolor=black, fontname=\"blod, bold\"];\n",head->name);
// 			fprintf(outputG2, "\"inter-%d\"[shape = circle, fillcolor=darkseagreen2, fontcolor=black, fontname=\"blod, bold\"];\n",head->name);			
// 		}
// 		head =head->suiv;
// 	}

// 	head= molecule->conformations;
// 	while(head){
// 		struct successor *newSucc=head->succ;
// 		while(newSucc){
// 			//Get type of change newSucc->transf
// 			char change[256];
// 			int cpt =0;
// 			strcpy(change,"\"(");
// 			if(bit_1(newSucc->transf,0)) strcat(change, "H-A, \0");
// 			if(bit_1(newSucc->transf,1)) strcat(change, "H-D, \0");
// 			if(bit_1(newSucc->transf,2)) strcat(change, "H-T, \0");
// 			if(bit_1(newSucc->transf,3)) strcat(change, "C-A, \0");
// 			if(bit_1(newSucc->transf,4)) strcat(change, "C-D, \0");
// 			if(bit_1(newSucc->transf,5)) strcat(change, "I-A, \0");
// 			if(bit_1(newSucc->transf,6)) strcat(change, "I-D, \0");
// 			if(bit_1(newSucc->transf,7)) strcat(change, "M-A, \0");
// 			if(bit_1(newSucc->transf,8)) strcat(change, "M-D, \0");
// 			strcat(change, ")\"\0");
// 			if(head->state=='C'){
// 				fprintf(outputG, "\"conf-%d\"->",head->name);
// 				fprintf(outputG2, "\"conf-%d\"->",head->name);
// 				cpt++;
// 			}
// 			else{
// 				fprintf(outputG, "\"inter-%d\"->",head->name);				
// 				fprintf(outputG2, "\"inter-%d\"->",head->name);				
// 			}	

// 			if(newSucc->conf->state=='C'){
// 				fprintf(outputG, "\"conf-%d\"",newSucc->conf->name);
// 				fprintf(outputG2, "\"conf-%d\"",newSucc->conf->name);
// 				cpt++;
// 			}
// 			else{
// 				fprintf(outputG, "\"inter-%d\"",newSucc->conf->name);				
// 				fprintf(outputG2, "\"inter-%d\"",newSucc->conf->name);				
// 			}	
// 			if(cpt==2){
// 				fprintf(outputG, "[label=%d, color=red];\n", newSucc->freq);
// 				fprintf(outputG2, "[label=%s, color=red];\n", change);				
// 			}
// 			else{
// 				fprintf(outputG, "[label=%d, color=black];\n", newSucc->freq);
// 				fprintf(outputG2, "[label=%s, color=black];\n", change);								
// 			}
// 			newSucc=newSucc->suiv;
// 		}
// 		head = head->suiv;
// 	}	
// 	fprintf(outputG, "}\n");
// 	fclose(outputG);
// 	fprintf(outputG2, "}\n");
// 	fclose(outputG2);

// //Draw the graph of transition with frequencies using graphViz
// 	char resFile[512];
// 	strcpy(resFile,"dot -Tpng  ");
// 	strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_graphFreq.gv -o \0"));
// 	strcat(resFile,get_FN(molecule->inputDir,molecule->inputFN,"_graphFreq.png \0"));
// 	system(resFile);
// // //Draw the graph of transition with frequencies using graphViz
// 	char resFile2[512];
// 	strcpy(resFile2,"dot -Tpng  ");
// 	strcat(resFile2,get_FN(molecule->inputDir,molecule->inputFN,"_graphChang.gv -o \0"));
// 	strcat(resFile2,get_FN(molecule->inputDir,molecule->inputFN,"_graphChang.png \0"));
// 	system(resFile2);
// }

/*========
Draw the distribution of conformations according to number of H-bonds formed, using gnuplot 
========*/
void draw_distribution(char inputFN[], char outputFN[]) {

// //Save data
//  , int vector[] , int size   
// FILE* outputD=fopen(inputFN,"w");
//     int i=0;
//     for(i=0;i<=size;i++)
//     	fprintf(outputD, "%d \t %d\n", i, vector[i]);
//     fclose(outputD); 

//Draw results using gnuplot file
    FILE* outputD=fopen("distribution.gplt","w");

	fprintf(outputD,"set title 'Distribution according H-bonds present'\n");
	fprintf(outputD,"\n");
	fprintf(outputD,"set grid\n");
	fprintf(outputD,"set xlabel '#H-bonds'\n");
	//fprintf(outputD,"set xrange [0:%d]\n",molecule->nbMolecule-1);
	//fprintf(outputD,"set ylabel '#conformations'\n");
	//fprintf(outputD,"set xrange [0:%d]\n",size+1);
	fprintf(outputD,"\n");
	fprintf(outputD, "set style fill solid 0.25 border -1\n");
	fprintf(outputD, "set boxwidth 0.4 relative\n");
	fprintf(outputD,"set term pdf\n");
	fprintf(outputD,"set output '%s'\n",outputFN); // difference between distance's files
	fprintf(outputD,"\n");
	fprintf(outputD,"plot \\\n");
	fprintf(outputD,"  '%s' using ($1-0.2):2 with boxes notitle ,\\\n",inputFN);
	//fprintf(outputD,"  '%s' using ($1+0.2):3  with boxes title '#conf. found',\\\n",inputFN);
	fprintf(outputD,"\n");
	fprintf(outputD,"\n");
	fclose(outputD);

//Generate the ".pdf" file
	char cmd[128];
	sprintf(cmd,"gnuplot %s","distribution.gplt");
	system(cmd);
}

/*========
Draw the mixed graph of conformation beloging to a path between two conformations.
We suppose that the difference is in H-bonds. The difference is only one HB with index indexdiff
========*/
void draw_pathConf(char *outputDir,struct MyModel *molecule,char *conf, int indexF,int indexdiff){
	//File to save
	char outputFN [512];
	sprintf(outputFN,"%s/%06d.pdf",outputDir,indexF);
	//GraphViz file
	FILE *outputG = fopen(get_FN(molecule->inputDir,molecule->inputFN,"_pathConf.gv\0"),"w");
	if(outputG == NULL){ printf("Can\'not open graphviz file \n"); exit(4);}  
	fprintf(outputG, "digraph G {\n");
	fprintf(outputG, "label=\"Conformation (%s)\";\n",conf);
	fprintf(outputG, "node [style=filled];\n");
	fprintf(outputG, "graph [bgcolor=transparent];\n");
	fprintf(outputG, "node [shape = circle, fontsize=10];\n");

	int i=0;
	int j=0;
	//Draw atoms except hydrogen atoms
	for(i=0;i<molecule->nbAtom;i++){
		//Oxygen in red
		if(strcmp(molecule->img[i].atomName,"O")==0) 
			fprintf(outputG, "\"%s%d\"[fillcolor=red, fontcolor=white, fontname=\"blod, bold\"];\n",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtom));
		//Neitrogen in blue
		if(strcmp(molecule->img[i].atomName,"N")==0) 
			fprintf(outputG, "\"%s%d\"[fillcolor=blue, fontcolor=white, fontname=\"blod, bold\"];\n",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtom));
		//Carbon in gray
		if(strcmp(molecule->img[i].atomName,"C")==0) 
			fprintf(outputG, "\"%s%d\"[fillcolor=gray32, fontcolor=white, fontname=\"blod, bold\"];\n",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtom));
		//Lithium in 
		if(strcmp(molecule->img[i].atomName,"Li")==0) 
			fprintf(outputG, "\"%s%d\"[fillcolor=blueviolet, fontcolor=white, fontname=\"blod, bold\"];\n",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtom));

	}
	//Make the graph according the covalent bonds
	for(i=0;i<molecule->nbAtom;i++){
		for(j=i+1;j<molecule->nbAtom;j++){
			if(molecule->covBond[i][j]!=0 ){ //covalent fixed or not and not leave
				if(molecule->img[i].atomType!='H' && molecule->img[j].atomType!='H')
					fprintf(outputG, "\"%s%d\"->\"%s%d\"[dir=none, color = black];\n",molecule->img[i].atomName, num_occ(molecule->img, i, molecule->nbAtom),molecule->img[j].atomName, num_occ(molecule->img, j, molecule->nbAtom) );
			}
		}		
	}
	//Make the graph according the H-bonds
	for(i=0;i<molecule->nbpairs;i++){
		if(conf[i]=='1' && i==indexdiff) //New hydrogen bond
			fprintf(outputG, "\"%s%d\"->\"%s%d\"[color = blue, label=%d];\n",molecule->img[molecule->cp[i*2]].atomName, num_occ(molecule->img, molecule->cp[i*2], molecule->nbAtom),molecule->img[molecule->cp[i*2+1]].atomName, num_occ(molecule->img, molecule->cp[i*2+1], molecule->nbAtom) ,i+1);					
		//if(conf[i]=='0' && i==indexdiff) //New hydrogen bond
			//fprintf(outputG, "\"%s%d\"->\"%s%d\"[color = yellow,  style=dashed, nodesep=0.5];\n",molecule->img[molecule->cp[i*2]].atomName, num_occ(molecule->img, molecule->cp[i*2], molecule->nbAtom),molecule->img[molecule->cp[i*2+1]].atomName, num_occ(molecule->img, molecule->cp[i*2+1], molecule->nbAtom) );					
		if(conf[i]=='1' && i!=indexdiff) //New hydrogen bond
			fprintf(outputG, "\"%s%d\"->\"%s%d\"[color = red, label=%d];\n",molecule->img[molecule->cp[i*2]].atomName, num_occ(molecule->img, molecule->cp[i*2], molecule->nbAtom),molecule->img[molecule->cp[i*2+1]].atomName, num_occ(molecule->img, molecule->cp[i*2+1], molecule->nbAtom),i+1 );					
	}
	//Close file
	fprintf(outputG, "}\n");
	fclose(outputG);
	//Draw the graph of transition using graphViz
	char resFile[550]="";
	sprintf(resFile,"%s %s -o %s","neato -Tpdf  ",get_FN(molecule->inputDir,molecule->inputFN,"_pathConf.gv\0"),outputFN);	
	system(resFile);
}


/*========
Save the list of hydrogen bonds for each conformation 
========*/
void save_Hbonds_list(struct MyModel *molecule){
	FILE* outputH=NULL;
  	if((outputH = fopen(get_FN(molecule->inputDir,molecule->inputFN,"_confHB.txt\0"),"w")) == NULL){ printf("Can\'not open file to save HB \n"); exit(-2);}
	struct confList* head= molecule->conformations;
	fprintf(outputH, "#conformations=%d\n",molecule->nbConf );
	while(head){
		fprintf(outputH, "conf %d :\n",head->name );
		for(int i=0;i<head->nbAtom;i++){
			if(head->img[i].atomIn){
				for(int j=i+1;j<head->nbAtom;j++){ //Update
					if(head->HB[i][j]==-5 && head->img[i].atomType!='H' && head->img[j].atomType!='H'){
						//check if there is proton transfer or not
						int D=-1;
						int A=-1;
						bool ptr=check_ptr(head,i,j,head->nbAtom,&D);							
							if(D==i)
								A=j;
							if(D==j)
								A=i;
							if(D==-1){
								printf("Error on the donor (3)\n");
								exit(-10);
							}	
						if(ptr){//proton transfer A-H-D
							fprintf(outputH, "%s%d-%s%d | %d \t %d\n",head->img[A].atomName, num_occ(head->img, A, head->nbAtom),head->img[D].atomName, num_occ(head->img, D, head->nbAtom),A,D);							
						}
						else{//D-H-A
							fprintf(outputH, "%s%d-%s%d | %d \t %d\n",head->img[D].atomName, num_occ(head->img, D, head->nbAtom),head->img[A].atomName, num_occ(head->img, A, head->nbAtom),D,A);							
						}
					}
				}
			}
		}		
		fprintf(outputH, "\n");
		head = head->suiv;	
	}
	fclose(outputH);
}
