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
 * \file hbond.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief Hydrogen bonds treatment
 * \details 
 *
 * 	This file contains all functions used to compute H-bonds and proton transfer.
 *
 * 	An H-bond is created if two conditions are satisfied:		
 *
 *  distance(Hydrogen, Acceptor)< Dr , where Dr is set in the constants file. 
 *
 * 	angle(Donor,Hydrogen, Acceptor) in interval [pi-alpha,pi+alpha], where alpha is set in the constants file.
 *	
 *	The angle can be ignored in the condition. Check the manual for more details.
 *	
 *  Two steps in order to identify the H-BONDS: 				
 *
 *	Step (I): calculate the orbits : a set of atoms that are	
 * 	at a given distance from a Hydrogen atom and which can 			
 * 	potentially form a H-bond. This idea allows to optimize 	
 * 	the number of comparisons and hence the computation time.		
 *
 * 	Step (II) : browse all hydrogen atoms and check if there are 
 *  H-bonds, which already exist, broken or if new H-bonds are formed
 *	with the current orbit.  	
 *  Proton transfer is identified at the Step(II).			
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
#include "hbond.h"

/*========================Functions=========================*/
/*========
Calculate the orbit of hydrogen
========*/
void get_orbits(struct MyModel *molecule,double alpha){
	int i=0;  //index on molecule->bondHList
	int m=0 ; //index in orbit			

  //Step (0): initialize the file and the list of hydrogen bonds
    init_Hbonds(molecule);
  //Step (I): calculate the orbits
  //COMPLEXITY : O(nbAtom x nbAtom)
	for ( i = 0; i < molecule->nbHAtom; i++){ 
		if (molecule->bondHList[i].numDon!=-1){
			int j=0; //index in molecule->img list
			int oldSize=molecule->bondHList[i].nbOrb;
			int oldOrbit[oldSize];
			copy_tab(molecule->bondHList[i].numOrb, oldOrbit,oldSize);
			molecule->bondHList[i].nbOrb=0;
			while(j<molecule->nbAtom){
				if(pair_valid(molecule,molecule->bondHList[i].numDon,j) ){
					double rH =0.00;
					if(molecule->sysType[0]=='w'){ //Use the PBC 
						rH = distance_pbc(molecule->img[molecule->bondHList[i].numH],molecule->img[j],molecule->tshVal.boxX,molecule->tshVal.boxY,molecule->tshVal.boxZ);
					}
					else{
						rH = distance(molecule->img[molecule->bondHList[i].numH].x,molecule->img[molecule->bondHList[i].numH].y,molecule->img[molecule->bondHList[i].numH].z,molecule->img[j].x,molecule->img[j].y,molecule->img[j].z);
					}				
					//check if the atoms are in the orbit , remove alpha for the interfaces.
					if(m==NB_MaxOrb){
						printf("Maximum number of atoms in orbit allowed (NB_MaxOrb= %d) is reached at snapshot %d \n",NB_MaxOrb,molecule->num_img);
						exit(-1);
					}					
					if (rH <= /*3*molecule->tshVal.distAHMax*/ (alpha-1)*molecule->tshVal.distAHMax){ 
						molecule->bondHList[i].numOrb[m]=j;
						m++;
						// printf("m=%d\n",m);
					}
				}
				j++;
			}
			//Save the size of the current orbit. 
			molecule->bondHList[i].nbOrb=m;
			m=0;
			//Update the donHacc list if the number of atoms has changed. 
			if(!molecule->nbAtomChange) 
				update_donHacc_orbit(molecule,i,oldOrbit,oldSize);
		}
		//else printf("H%d atom of index %d doesn't have match at image %d \n",num_occ(molecule->img, i, molecule->nbAtom),i, molecule->num_img);	
	}
}

/*========
Update donHacc if orbit has changed
========*/
void update_donHacc_orbit(struct MyModel *molecule,int indexH, int *oldOrbit, int size){
	int i ;
	for(i=0; i<size; i++){
		if(!exit_atom_inOrbit(molecule->bondHList[indexH].numOrb , molecule->bondHList[indexH].nbOrb,oldOrbit[i])){ //The atom i doesn't belongs to the orbit any more
			if(molecule->donHacc[ molecule->bondHList[indexH].numH ][oldOrbit[i]]==1 || molecule->donHacc[ molecule->bondHList[indexH].numH ][oldOrbit[i]]==0 ){
				//Update donHacc , delete the current hydrogen bond because the atom is out of the orbit
				molecule->donHacc[ molecule->bondHList[indexH].numH ][oldOrbit[i]]=-1;
				molecule->donHacc[oldOrbit[i]][ molecule->bondHList[indexH].numH ]=-1;

				molecule->donHacc[ molecule->bondHList[indexH].numDon ][oldOrbit[i]]=-1;
				molecule->donHacc[oldOrbit[i]][ molecule->bondHList[indexH].numDon ]=-1;

				molecule->donHacc[ molecule->bondHList[indexH].numDon ][molecule->bondHList[indexH].numH]=-1;
				molecule->donHacc[molecule->bondHList[indexH].numH][ molecule->bondHList[indexH].numDon ]=-1;

			}
		}
	}

}

/*========
Check if the atom index belongs to the orbit
========*/
bool exit_atom_inOrbit(int *tabOrbit, int  size, int index){
	int i=0;

	while(i<size)
		if(tabOrbit[i]==index)
			return true;
		else
			i++;
	return false;
}

/*========
Initialize the H-bonds table
========*/
void init_Hbonds(struct MyModel *molecule){
	int i;  //index on molecule->img
	int k=0; //index on bondHList

	//COMPLEXITY : O(nbHAtom)
	for(i=0;i<molecule->nbAtom;i++){
		if (molecule->img[i].atomType=='H')
		{
			molecule->bondHList[k].numH = i;
			int oldDon = molecule->bondHList[k].numDon;
			molecule->bondHList[k].numDon = get_don(i,molecule); //Get the index of the atom with which 'H' is bonded
			if(molecule->num_img==0 || molecule->nbAtomChange) {
				molecule->bondHList[k].sens ='0';//Initialize the sense  just at first time	
				molecule->bondHList[k].numAcc = -1;
				molecule->bondHList[k].state=-1;
			}
			else
				if(oldDon != molecule->bondHList[k].numDon) //Covalent has changed
					molecule->bondHList[k].numAcc = -1; //Reset the acceptor
			k++;
		}	
	}
}

/*========
Get the index of donor related to the hydrogen atom with index H
========*/
int get_don(int H,struct MyModel *molecule){
	int j=0; 
	while(j<H ) 
	if(molecule->covBond[j][H] ==1) return j ; else j++;
	j++;
	while(j<molecule->nbAtom ) 
	if(molecule->covBond[H][j] ==1) return j ; else j++;
return -1;
}


/*========
Check if an hydrogen bond can be created
========*/
bool verify_Hbond(struct MyModel *molecule,int D, int H, int A,int *state, char *sens){
	//There is no donor 
	if(D==-1)  
		return false;
	//Distances between H, donor & acceptor
	double DH  = distance(molecule->img[D].x,molecule->img[D].y,molecule->img[D].z,molecule->img[H].x,molecule->img[H].y,molecule->img[H].z); //distance beteween hydrogen and donor
	double AH  = distance(molecule->img[A].x,molecule->img[A].y,molecule->img[A].z,molecule->img[H].x,molecule->img[H].y,molecule->img[H].z); //distance beteween hydrogen and acceptor
	double DA  = 0.00; //Distance between the three atoms
	double ang = 0.00;
	bool verif=true;
	
	if(molecule->sysType[0]!='w'){
		//Distance Donor-Acceptor
		DA=distance(molecule->img[D].x,molecule->img[D].y,molecule->img[D].z,molecule->img[A].x,molecule->img[A].y,molecule->img[A].z); //distance beteween acceptor and donor	
		//Angle between Donor, H and the acceptor
		ang = (angle(molecule->img[D].x,molecule->img[D].y,molecule->img[D].z,molecule->img[H].x,molecule->img[H].y,molecule->img[H].z,molecule->img[A].x,molecule->img[A].y,molecule->img[A].z))*180/PI;
	}
	else{
		//Distance Donor-Acceptor
		DA=distance_pbc(molecule->img[A], molecule->img[D],molecule->tshVal.boxX,molecule->tshVal.boxY,molecule->tshVal.boxZ);
	    //Angle between Donor, H and the acceptor
		ang = angle_pbc(molecule->img[D],molecule->img[H], molecule->img[A],molecule->tshVal.boxX,molecule->tshVal.boxY,molecule->tshVal.boxZ)*180/PI;
	}
	
	//Four states for the H-bond : No H-bond, Donor and acceptor are in position of an H-bond, H-bond created or a proton transfer
	if(molecule->anghbond=='1'){
		if(ang-molecule->tshVal.angDHAMin< EPSILON)//molecule->tshVal.angDHAMin
			verif=false;
	}
	if(molecule->sysType[0]=='w'){//DA<3.2 
		if(DA - molecule->tshVal.distDAMax > EPSILON)
			verif = false;
	}	
	else{ //AH<2.0  or 2.3 
		if(AH - molecule->tshVal.distAHMax > EPSILON)
			verif = false;
	}
	if(verif){ //H-bond created  
		*state =1;
		if(molecule->sysType[0]=='i' || molecule->sysType[0]=='f' ){ //Check proton transfer for systems different of droplet
			if(DH>molecule->tshVal.prtMarg*DA){ //Proton transfer
				if(*sens=='0')	*sens='1';	else *sens='0';
				if(DH>molecule->tshVal.distAHMax)
					*state =0; //The H-bond doesn't exist in the other sense
			}			
		}
		return true;			
	}
	else{ 
	    //Updated March 2017 : we forget atoms in position of hydrogen bonds
		// verif=true;
		// if(molecule->anghbond=='1'){
		// 	if(ang<molecule->parmtrs[3]) //angle < angle_min
		// 		verif=false;
		// }
		// if(AH<molecule->parmtrs[2] && verif ){ //We keep the H-bond (dist< dist_max and ang > ang_min)
		// 	*state=0; //In position of hydrogen bond
		// 	*sens='0';
		// 	return true;
		// }
		// else{ //No hydrogen bond can be created atoms are far and angles (if applied) are not potential
		*state=-1;
		*sens='0';
		return false;
		// } 
	}
	//return false;
}

/*========
Check if a new H-bond is created with atom of position index
========*/
int check_new_Hbond(struct MyModel *molecule, int index){
	int i=0; int k = -1 ;
	//COMPLEXITY: O(nbOrb)
	molecule->bondHList[index].state=-1;
	molecule->bondHList[index].sens='0';
	//We choose the closest atom to the hydrogen
	while(i<molecule->bondHList[index].nbOrb){
		int state=-1;	
		char sens='0';
		if(verify_Hbond(molecule,molecule->bondHList[index].numDon,molecule->bondHList[index].numH,molecule->bondHList[index].numOrb[i],&state,&sens)){ // an H bond is created
			if(k==-1){
				molecule->bondHList[index].state=state;
				molecule->bondHList[index].sens=sens;
				k= molecule->bondHList[index].numOrb[i];
			}
			else
				if(distance(molecule->img[molecule->bondHList[index].numOrb[i]].x,molecule->img[molecule->bondHList[index].numOrb[i]].y,molecule->img[molecule->bondHList[index].numOrb[i]].z,molecule->img[molecule->bondHList[index].numH].x,molecule->img[molecule->bondHList[index].numH].y,molecule->img[molecule->bondHList[index].numH].z)<distance(molecule->img[k].x,molecule->img[k].y,molecule->img[k].z,molecule->img[molecule->bondHList[index].numH].x,molecule->img[molecule->bondHList[index].numH].y,molecule->img[molecule->bondHList[index].numH].z)){
					molecule->bondHList[index].state=state;
					molecule->bondHList[index].sens=sens;
					k=molecule->bondHList[index].numOrb[i];
				}
			i++;
		}					
		else i++;
	}
	return k;
}

/*========
Proton transfer leads to a change on the covalent bond
========*/
void change_cov(struct MyModel *molecule, int index){
	if(molecule->bondHList[index].numH<molecule->bondHList[index].numAcc){
		molecule->covBond[molecule->bondHList[index].numH][molecule->bondHList[index].numAcc]=0 ;
	}
	else {
		molecule->covBond[molecule->bondHList[index].numAcc][molecule->bondHList[index].numH]=0 ;
	}
	
	if(molecule->bondHList[index].numH <molecule->bondHList[index].numDon){
		molecule->covBond[molecule->bondHList[index].numH][molecule->bondHList[index].numDon]=1 ;
	}
	else{
		molecule->covBond[molecule->bondHList[index].numDon][molecule->bondHList[index].numH]=1 ;		
	} 
}

/*========
Update the H-bonds states (Donor,Acceptor) of hydrogen atom of position index 
========*/
void update_AH(struct MyModel *molecule, int index){
	int D = molecule->bondHList[index].numDon;
	int A = molecule->bondHList[index].numAcc;
	int H = molecule->bondHList[index].numH;
	int sens = molecule->bondHList[index].sens ;
	int state = molecule->bondHList[index].state;
	if(sens=='0'){
		//imp_hbond(molecule,D,A,sens);
		molecule->donHacc[D][A]=state;
		if(state==-1)
			molecule->donHacc[A][H]=-1;
		else
			molecule->donHacc[A][H]=2;			
	}
	else{
		//imp_hbond(molecule,A,D,sens);
		molecule->donHacc[A][D]=state;
		if(state==-1)
			molecule->donHacc[D][H]=-1;
		else
			molecule->donHacc[D][H]=3;			
	}
}

/*========
Get index of atom with index2 in the list of orbit of hydrogen atom index
========*/
int get_indexP(struct MyModel *molecule, int index, int index2){
	int i=0; 
	while(i<molecule->bondHList[index].nbOrb &&  molecule->bondHList[index].numOrb[i]!=index2) i++;
	if(i<molecule->bondHList[index].nbOrb)  return i; else return -1;
}

/*========
Mark the atoms involved in a H-bond
========*/
void imp_hbond(struct MyModel *molecule, int i, int j, char sens){
// i  : donor, j: hydrogen , k: acceptor 
	switch(sens){
		case '0' : //donor
		if(j> i)
			molecule->covBond[j][i]=1;
		else 	
			molecule->covBond[i][j]=1;
		break;

		case '1' : //acceptor
		if(j> i)
			molecule->covBond[j][i]=-1;
		else 	
			molecule->covBond[i][j]=-1;
		break;

		case '2' : //bond broken
		if(j> i)
			molecule->covBond[j][i]=0;
		else 	
			molecule->covBond[i][j]=0;
		break;

		default:printf("ERROR on the sens of the H-bond\n");
	}
}

/*========
Dynamic analysis of H-bonds (proton transfer included)
========*/
bool get_Hbonds_dynamics(struct MyModel *molecule){
	int **copy_donHacc= allocate_matrix(molecule->nbAtom,molecule->nbAtom,-1,"copy_donHacc");
	if(!molecule->nbAtomChange){
		// copy_donHacc= allocate_matrix(molecule->nbAtom,molecule->nbAtom,-1,"copy_donHacc");
		copy_matrix(molecule->donHacc,copy_donHacc,molecule->nbAtom);	
	}
	else{
		molecule->nbHbond =0;
	}
	int i=0 ;

//Browse the list of Hydrogen atoms
	for (i = 0; i < molecule->nbHAtom; i++){
		if(molecule->bondHList[i].numDon != -1){
			int hbond[3];
			hbond[0]=molecule->bondHList[i].numDon; //Donor
			hbond[1]=molecule->bondHList[i].numH; //Hydrogen
			//No acceptor at the moment
			int curAcc=-1;
			int cpt=0;
			double curDist=0.00;
			double curAng=0.00;
			char curSens='0';
			//Browse atoms in the orbit of the current hydrogen atom
			int j=0;
			while(j<molecule->bondHList[i].nbOrb && curSens=='0'){
			//Check if there is hydrogen bond with current orbit
				hbond[2]=molecule->bondHList[i].numOrb[j];
				molecule->bondHList[i].numAcc=hbond[2];
				int state=-1;
				char sens='0';
				if(verify_Hbond(molecule,hbond[0],hbond[1],hbond[2],&state,&sens)){
					if(sens=='1'){//Check if there is proton transfer
						if(molecule->donHacc[hbond[1]][hbond[1]]==hbond[2]){//Hbonds becomes to nature sense
							sens='0';
							update_donHacc(molecule->donHacc,hbond[2],hbond[1],hbond[0],state,sens,molecule->nbAtom); //HH=-1, DH/HD = 4, AH/HA=state 
						} 
						else //New proton transfer 
							update_donHacc(molecule->donHacc,hbond[0],hbond[1],hbond[2],state,sens,molecule->nbAtom); //HH=D, DH/HD = state, AH/HA=3  
						curSens='1';
						curAcc = hbond[2]; //PTR stop 		
					}
					else{ //No PTR, so check if a H-Bond is formed
						if(molecule->donHacc[hbond[1]][hbond[1]]==-1){ //H-bond in nature sense
							double orbDist=distance(molecule->img[hbond[2]].x,molecule->img[hbond[2]].y,molecule->img[hbond[2]].z,molecule->img[hbond[1]].x,molecule->img[hbond[1]].y,molecule->img[hbond[1]].z); //distance between hydrogen and acceptor
							double orbAng=(angle(molecule->img[hbond[0]].x,molecule->img[hbond[0]].y,molecule->img[hbond[0]].z,molecule->img[hbond[1]].x,molecule->img[hbond[1]].y,molecule->img[hbond[1]].z,molecule->img[hbond[2]].x,molecule->img[hbond[2]].y,molecule->img[hbond[2]].z))*180/PI;
							if(curAcc==-1 || state==0){ //First Hydrogen bond created
								update_donHacc(molecule->donHacc,hbond[0],hbond[1],hbond[2],state,sens,molecule->nbAtom); //DH/HD = 4, AH/AH=state (0/1)
								cpt++;
								if(state==1){
									curAcc=hbond[2];
									curDist=orbDist;
									curAng=orbAng;
								}
							}
							else{ //There is already Hydorgen bond && state=1
								if((curDist- orbDist>0.10) && (curAng - orbAng <10.00)){ //the new hydrogen bond is better than the current one
									update_donHacc(molecule->donHacc,hbond[0],hbond[1],curAcc,-1,sens,molecule->nbAtom); //Becomes potential HB (ignore 0 at the moment, jan 2018)
									update_donHacc(molecule->donHacc,hbond[0],hbond[1],hbond[2],state,sens,molecule->nbAtom); //Becomes HB
									//Change the current acceptor
									curAcc=hbond[2]; 
									curDist=orbDist;
									curAng=orbAng;
									cpt++;
								}
								else{//Keep the H-bond that already exists and put the new one as a potential H-bond.
									update_donHacc(molecule->donHacc,hbond[0],hbond[1],hbond[2],-1,sens,molecule->nbAtom); //DH/HD = 4, AH/HA=0 (ignore 0 at the moment, jan 2018)						
								}
							}
						}
						else{ //Check if the orbit is the old donor
							if(molecule->donHacc[hbond[1]][hbond[1]]==hbond[2]){
								update_donHacc(molecule->donHacc,hbond[2],hbond[1],hbond[0],state,'1',molecule->nbAtom); //DH/HD =state
								if(state==1) curAcc=hbond[2]; 				
							}
						}							
					}
				}
				else{ //H-bonds broken or doesn't exist 
					if(molecule->donHacc[hbond[1]][hbond[2]]!=-1){ //Hbond broken
						// if(molecule->sysType[0]!='w' && verif_pot_Hbond(molecule,hbond[0],hbond[1],hbond[2])){ //check if we can potentially form a hydrogen bond
						// 	if(molecule->donHacc[hbond[1]][hbond[1]]==hbond[2]){
						// 		update_donHacc(molecule->donHacc,hbond[2],hbond[1],hbond[0],0,'1',molecule->nbAtom);//state =-1 -> DH/HD = -1 : check sense
						// 	}
						// 	else{
						// 		update_donHacc(molecule->donHacc,hbond[0],hbond[1],hbond[2],0,sens,molecule->nbAtom);//state =-1 -> AH/HA = -1 : check sense
						// 	}

						// }
						// else{
							if(molecule->donHacc[hbond[1]][hbond[1]]==hbond[2]){
									update_donHacc(molecule->donHacc,hbond[2],hbond[1],hbond[0],-1,'1',molecule->nbAtom);//state =-1 -> DH/HD = -1 : check sense
								}
							else{
								update_donHacc(molecule->donHacc,hbond[0],hbond[1],hbond[2],-1,sens,molecule->nbAtom);//state =-1 -> AH/HA = -1 : check sense
							}
							save_event(molecule,i,0,'-','a');	
						// }
					}
				}

				j++;

			}
			
			if(!molecule->nbAtomChange){			
				if(curAcc!=-1 && (copy_donHacc[molecule->bondHList[i].numH][curAcc]==-1 || copy_donHacc[molecule->bondHList[i].numH][molecule->bondHList[i].numH]!=molecule->donHacc[molecule->bondHList[i].numH][molecule->bondHList[i].numH])){//New hydrogen bond
					if(curSens=='0'){ //New hydrogen bond without PTR
						molecule->bondHList[i].numAcc=curAcc;
						molecule->bondHList[i].state=1;
						//molecule->bondHList[i].sens='0';
						//Save that a new H-bond has been created
						save_event(molecule,i,0,'+','a');
					}
					else{ //There is a proton transfer
						//Update donor/acceptor & state of the other side of the HB
						int indexOrb=get_indexP(molecule,i,curAcc);
						molecule->bondHList[i].numOrb[indexOrb]=molecule->bondHList[i].numDon; //donor becomes an atom of orbit
						molecule->bondHList[i].numAcc=molecule->bondHList[i].numDon; 
						molecule->bondHList[i].numDon=curAcc;//We change donor
						if(molecule->bondHList[i].sens=='1') molecule->bondHList[i].sens='0'; else molecule->bondHList[i].sens='1';
						//Change cov
						change_cov(molecule,i);
						//Save that a proton transfer has been created
						save_event(molecule,i,0,'p','a');
						//Put all states with this hydrogen atom at -1
						update_donHacc_ptr(molecule->donHacc,molecule->bondHList[i].numOrb,molecule->bondHList[i].numH,indexOrb,molecule->bondHList[i].nbOrb);
					}
				}
			}	
		}
	}

	//Check the number of acceptors per HB
	check_nb_acc_per_hbond(molecule);
	//Check changes in hydrogen bonds that occur 
	if(!molecule->nbAtomChange){
		if(verif_change_Hbonds(molecule, molecule->donHacc,copy_donHacc,molecule->nbAtom)){ //There is a change in HB
			free_matrix(copy_donHacc,molecule->nbAtom);
			return true;				
		}
		else{//No change we keep the current donHacc
			free_matrix(copy_donHacc,molecule->nbAtom);
			return false;		
		} 		
	}
	else{
		free_matrix(copy_donHacc,molecule->nbAtom);
		return true;
	}
}

/*========
Verify if three atoms D,H and A can potentially form a hydrogen bond
========*/
bool verif_pot_Hbond(struct MyModel *molecule, int D,int H, int A ){
	double DH  = distance(molecule->img[D].x,molecule->img[D].y,molecule->img[D].z,molecule->img[H].x,molecule->img[H].y,molecule->img[H].z); //distance beteween hydrogen and donor
	double AH  = distance(molecule->img[A].x,molecule->img[A].y,molecule->img[A].z,molecule->img[H].x,molecule->img[H].y,molecule->img[H].z); //distance beteween hydrogen and acceptor
	double DA= DH+AH;

	if(DA - molecule->tshVal.distDHAMax < EPSILON){
		return true;
	}
	else{
		return false;
	}
}

/*========
Check that each hydrogen atom has at most 2 hydrogen bond 
========*/
void check_nb_acc_per_hbond(struct MyModel *molecule){
	int i;
	int j;
	for(i=0;i<molecule->nbAtom;i++){
		int cptD=0;
		int cptA=0;
		int dIndex[4];
		init_tab(dIndex,4,-1);

		if(molecule->img[i].atomIn){
			j=0;
			//Research by lines
			while (j<i){
				if(molecule->donHacc[i][j]==-5  && molecule->img[i].atomType!='H' && molecule->img[j].atomType!='H'){
					if(get_HHB(molecule->donHacc, i,j,molecule->nbAtom)==i)
						cptD++;
					else{
						dIndex[cptA]=j;
						cptA++;
					}
				
				}
				j++;
			}
			if(j==i) j++;
			//Research by columns
			while (j<molecule->nbAtom){
				if(molecule->donHacc[i][j]==-5  && molecule->img[i].atomType!='H' && molecule->img[j].atomType!='H'){				
					if(get_HHB(molecule->donHacc, i,j,molecule->nbAtom)==i)
						cptD++;
					else{
						dIndex[cptA]=j;
						cptA++;
					}
				}
				j++;
			}
			while(cptA>2){
				printf("atom %d has formed %d Hydrogen bonds at snapshot %d ",i,cptA,molecule->num_img);
				//Update the donHacc matrix :
				int maxIndex=-1;
				double maxDist=0.00;
				for(j=0;j<cptA;j++){
					if(dIndex[j]!=-1){
						double dist=distance(molecule->img[i].x,molecule->img[i].y,molecule->img[i].z,molecule->img[dIndex[j]].x,molecule->img[dIndex[j]].y,molecule->img[dIndex[j]].z);
						if(maxDist<dist){
							maxDist=1;
							maxIndex=j;
						}
					}
				} 
				molecule->donHacc[i][dIndex[maxIndex]]=-1;
				molecule->donHacc[dIndex[maxIndex]][i]=-1;
				dIndex[maxIndex]=-1;
				cptA--;

			}
		}
	} 
}
/*========
Update the matrix of donor-hydrogen-acceptor 
========*/
void update_donHacc(int **donHacc,int D,int H, int A, int state,char sens, int size){

	if(sens=='1'){ //Proton transfer
		//H-Donor
		donHacc[H][H]=D;
		donHacc[D][H]=state;
		donHacc[H][D]=state;	
		//H-Acceptor
		donHacc[A][H]=-3;
		donHacc[H][A]=-3;
	}
	else{
		//Donor
		donHacc[H][H]=-1;		
		donHacc[D][H]=-4;
		donHacc[H][D]=-4;
		//Acceptor
		donHacc[A][H]=state;
		donHacc[H][A]=state;
	}
	//Donor-Acceptor
	if(state>=0){
		donHacc[D][A]=-5;
		donHacc[A][D]=-5;
	}
	else{
		if(get_Hbond_sens(donHacc,size,A,D)==-1){
			donHacc[D][A]=-1;
			donHacc[A][D]=-1;						
			donHacc[D][H]=-4;
			donHacc[H][D]=-4;		
		}
	}
} 

/*========
Verify if there has been a change in the hydrogen bonds 
========*/
bool verif_change_Hbonds(struct MyModel *molecule,int **donHacc,int **copy_donHacc,int size){
	int i =0 ;
	int j=0 ;
	bool hbChang=false;
	for (i = 0; i < size; i++){
		for (j = i+1; j < size; j++){
			if(copy_donHacc[i][j]!=donHacc[i][j] && (molecule->img[i].atomType=='H' || molecule->img[j].atomType=='H') ) { //I or J is Hydrogen atom
				int hbond[3];
				if(molecule->img[i].atomType=='H'){
					hbond[1]=i;
					hbond[2]=j;					
				} 
				else{
					hbond[1]=j;
					hbond[2]=i;										
				}
				if(copy_donHacc[hbond[1]][hbond[2]]==-3 || donHacc[hbond[1]][hbond[2]]==-3){//Proton transfer
					if(donHacc[hbond[1]][hbond[2]]==-3){//From sens =0 to sens =1
						hbond[0]=donHacc[hbond[1]][hbond[1]];
						// add_Hbond(molecule,hbond,3); //3 for blue colour
					}
					else{//From sens =1 to sens =0
						hbond[0]=copy_donHacc[hbond[1]][hbond[1]];//Get donor
					}
					
					if(donHacc[hbond[1]][hbond[2]]==-3)
						add_Hbond(molecule,hbond,3);
					else								
						add_Hbond(molecule,hbond,-1);
					
					hbChang=true;									
				}
				else{ //Hydrogen bond formed or broken
					if(copy_donHacc[hbond[1]][hbond[2]]!=-4 && donHacc[hbond[1]][hbond[2]]!=-4){
						if(donHacc[hbond[1]][hbond[1]]==-1){ //Initial sens 0
							hbond[0]=get_donHB(donHacc,hbond[1],size);//Get donor index		
							if(copy_donHacc[hbond[1]][hbond[2]]==-1 || donHacc[hbond[1]][hbond[2]]==-1){//New hydrogen bond or hydrogen bond broken
								hbChang=true;
							}
							add_Hbond(molecule,hbond,donHacc[hbond[1]][hbond[2]]);															
						}
						else{ //PTR  : donHacc[hbond[1]][hbond[1]]==D
							hbond[0]=donHacc[hbond[1]][hbond[1]];
							if(hbond[0]!=hbond[2]){
								if(copy_donHacc[hbond[1]][hbond[2]]==-1 || donHacc[hbond[1]][hbond[2]]==-1){//New hydrogen bond or hydrogen bond broken
									// if(copy_donHacc[hbond[1]][hbond[2]]==-1){ 
									// 	add_Hbond(molecule,hbond,donHacc[hbond[1]][hbond[2]]);	
									// }
									// else{//donHacc[hbond[1]][hbond[2]]==-1 HB broken
									// 	add_Hbond(molecule,hbond,-1);	
									// }
									hbChang=true;
								}
								// else{ //From 0 to 1 or 1 to 0
									if(donHacc[hbond[1]][hbond[2]]==-3)
										add_Hbond(molecule,hbond,3);
									else								
										add_Hbond(molecule,hbond,-1);
								// }														
							}						
						}
					}
				}
			}				
		}
	}
	return hbChang;
}

/*========
If there has been a proton transfer, set all atoms in the orbit (different of indexOrb) of indexH to -1.
========*/
void update_donHacc_ptr(int **donHacc,int *numOrb, int indexH,int indexOrb, int size){
	int i=0;
	for(i=0;i<size;i++){
		//if(donHacc[indexH][numOrb[i]]!=4 && donHacc[indexH][numOrb[i]]!=-1){
		if(i!=indexOrb){
			donHacc[indexH][numOrb[i]]=-1;
			donHacc[numOrb[i]][indexH]=-1;			
		}	
		//}
	}
}

/*========
Add dynamics of H-bond bondH to molecule->bondHdyn
========*/
void add_Hbond(struct MyModel* molecule,int bondH[3], int state){
	int i = exist_Hbond(molecule,bondH);
	if (i != -1){ // The H-bond already exists
		int j=molecule->bondHdyn[i].nChange;
 		if(molecule->bondHdyn[i].imgList[j*3-1]==-1)
			molecule->bondHdyn[i].imgList[j*3-1]=molecule->num_img;
		if(state!=-1){
			//j++;// commented for ylene, to be changed 	
			if(j==NB_PERIOD){
				printf("Maximum number of periods allowed for HD (NB_PERIOD=%d) is reached at snapshot %d \n",NB_PERIOD,molecule->num_img);
				exit(-1);
			}
			if(molecule->sysType[0]!='w' ){
				molecule->bondHdyn[i].imgList[j*3-3]= state;
				molecule->bondHdyn[i].imgList[j*3-2]= molecule->num_img;
				molecule->bondHdyn[i].imgList[j*3-1]=-1; 
				molecule->bondHdyn[i].nChange=j;	
			}
		}
	}
	else{ //the H-bond doesn't exist yet, it's inserted at the end
		if(molecule->nbHbond==NB_MaxHBond){
			printf("Maximum number of Hydrogen bonds allowed (NB_MaxHBond=%d) is reached at snapshot %d \n",NB_MaxHBond,molecule->num_img);
			exit(-1);
		}
		add_new_Hbond(molecule,bondH,state);
	}
}

/*========
Add a new H-bond to the "bondHdyn" list
========*/
void add_new_Hbond(struct MyModel* molecule,int bondH[3], int state){
	molecule->bondHdyn[molecule->nbHbond].bondH[0]= bondH[0];//D
	molecule->bondHdyn[molecule->nbHbond].bondH[1]= bondH[1];//H
	molecule->bondHdyn[molecule->nbHbond].bondH[2]= bondH[2];//A
	molecule->bondHdyn[molecule->nbHbond].nChange=1;
	molecule->bondHdyn[molecule->nbHbond].imgList[0]= state;
	molecule->bondHdyn[molecule->nbHbond].imgList[1]= molecule->num_img;
	molecule->bondHdyn[molecule->nbHbond].imgList[2]=-1; 
	molecule->nbHbond++;
}

/*========
Check if the H-bond bondH already exists in bondHdyn list
========*/
int exist_Hbond(struct MyModel* molecule, int bondH[3]){
	int i=0;
	bool found=false;	
	while(i<molecule->nbHbond && !(found)){
		if (molecule->bondHdyn[i].bondH[0]==bondH[0] && molecule->bondHdyn[i].bondH[1]==bondH[1] && molecule->bondHdyn[i].bondH[2]==bondH[2]) found =true; else i++;
	}  
	if(found) return i; else return -1;
}

/*========
Clean Hbonds according percentage of appearance and disappearance.
Note : this function is not used yet in our program.
========*/
void clean_Hbond_dynamics(struct MyModel* molecule){
  	int i=0;
  	//
  	int presHB[molecule->nbHbond];
  	init_tab(presHB,molecule->nbHbond,0);

  	int absHB[molecule->nbHbond];
  	init_tab(absHB,molecule->nbHbond,0);

	while(i<molecule->nbHbond){
		int j=0;
		for(j=0;j<molecule->bondHdyn[i].nChange;j++){
			if(molecule->bondHdyn[i].imgList[j*3]!=-1){	
				if(j==0 && molecule->bondHdyn[i].imgList[j*3+1]!=0){ //The H-bond is not present from the beginning 
					absHB[i] += (int) pow(molecule->bondHdyn[i].imgList[j*3+1]-0,2);
				}
				if(j<molecule->bondHdyn[i].nChange-1){
					presHB[i] += (int) pow(molecule->bondHdyn[i].imgList[j*3+2]-molecule->bondHdyn[i].imgList[j*3+1],2);
					absHB[i] +=(int) pow(molecule->bondHdyn[i].imgList[(j+1)*3+1]-molecule->bondHdyn[i].imgList[j*3+2],2);
				}
				else{
					if(molecule->bondHdyn[i].imgList[j*3+2]==-1){
						presHB[i] += (int) pow(molecule->nbMolecule-molecule->bondHdyn[i].imgList[j*3+1],2);
					}
					else{
						presHB[i] += (int) pow(molecule->bondHdyn[i].imgList[j*3+2]-molecule->bondHdyn[i].imgList[j*3+1],2);
						absHB[i] +=(int) pow(molecule->nbMolecule-molecule->bondHdyn[i].imgList[j*3+2],2);
					}
				}

			}
		}	
		printf("HB %d is present %lf and absent %lf \n",i, (double)presHB[i]*100/pow(molecule->nbMolecule,2), (double)absHB[i]*100/pow(molecule->nbMolecule,2) );		//*100/pow(molecule->nbMolecule,2)
		i++;
	}
}
