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
 * \file usage.c
 * \author Bougueroua Sana (2014-2022), Ali Hashemi (2020-2022)
 * \version 1.0
 * \date  2022
 * \brief usage documentation.
 *
 * This file contains a mini help of usage of the program.
 */
/*========================Libraries=========================*/ 
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include "constants.h"
#include "struct.h"
#include "comfunct.h"
#include "usage.h"

/*========================Functions=========================*/

/*========
Display the usages of the different parts of the Getaway program
========*/
void display_usage( char *argv[]) {
            fprintf(stderr, "Syntax error %s",str_sub(strrchr(argv[0], '/'),1,strlen(strrchr(argv[0], '/')))); 
    // if(strcmp(str_sub(strrchr(argv[0], '/'),1,strlen(strrchr(argv[0], '/'))),"singanalysis")==0){
            fprintf(stderr,"Usage: singanalysis  -s file [-i|-f|-c] [-d val] [-h car] [-a] [-l] [-t] [-g] \n");
            fprintf(stderr, "\t -s : single analysis\n" );
            fprintf(stderr," \t file : trajectory file for single analysis. It must have the extension: .xyz or .xmol \n");
            fprintf(stderr, "\t [-i|-f|-c] : the type of the molecular system to analyse:\n");
            fprintf(stderr, "\t\t -i: MD of an isolated peptide in gas phase (by default: only hydrogen bonds dynamics are analysed)\n");
            fprintf(stderr, "\t\t -f: Follow a collision leading to the fragmentation of a peptide (by default: only covalent bonds dynamics are analysed)\n");
            fprintf(stderr, "\t\t -c: MD of an isolated cluster in gas phase (by default: only hydrogen bonds dynamics are analysed)\n");
            fprintf(stderr, "\t [-d val] : dynamics to take into account:\n");
            fprintf(stderr, "\t\t It is 4 bits number that let you choose which kinds of bonds do you would take into account for the conformational dynamics analysis. These bits 1 to 4 are assigned to variables as follow: \n");
            fprintf(stderr, "\t\t\t Bit     Variables                       Value \n");
            fprintf(stderr, "\t\t\t  1       H-Bonds                     \t\t   1 \n");
            fprintf(stderr, "\t\t\t  2       Covalent bonds               \t\t  2 \n");
            fprintf(stderr, "\t\t\t  3       Intermolecular interactions    \t\t 4 \n");
            fprintf(stderr, "\t\t\t  4       Organometallic interactions    \t\t 8 \n");
            fprintf(stderr, "\t\t To only analyze H-bonds, choose 1. To analyze H-bonds with intermolecular interactions, choose 1+4=5. To analyse H-bonds, covalent bonds and intermolecular interactions, choose 1+2+4=7, etc.\n" );
            fprintf(stderr, "\t [-h 0/1] : 1 to take into account angle in hydrogen bonds, 0 else (by default 1)\n");
            fprintf(stderr, "\t [-a] : save the real energy average computed from the simulation for each conformation identified\n");        
            fprintf(stderr, "\t [-l] : this field aims to perform a partial analysis of the molecular system. Only the fragments containing metals (_Mn_ in our case) are taking into account\n");        
            fprintf(stderr, "\t [-t] : save the atoms involved in the dynamics of the analyzed molecular system.\n");        
            fprintf(stderr, "\t [-g] : analyze the list of fragments that have been observed along the trajectory\n");        
    //}
    exit(1);    
}

/*========
Return the default values for options used in the Getaway program 
========*/
struct MyOpt default_val() {
    struct MyOpt p;

    p.pathFile =""; // Full path for the input file
    p.s_on=0; // Single analysis module
    p.isol_on=0; // Isolated peptide  
    p.cid_on=0;  // Collision CID
    p.clst_on=0; // Clusters 
    p.bonds_on=0; // Bonds dynamics
    p.bonds_val = 0; // Bonds taken into account 
    p.ang_HB_on=true; // Angle for HBs
    p.eng_anal=false;
    p.path_on=0;
    p.frg_on=0;
    p.total_on=0;
    return p;
}

/*========
Read arguments introduced by user and set the corresponding options 
========*/
struct MyOpt usage(int argc, char *argv[]) {
    char c=' ';
    struct MyOpt p;
    p = default_val();
    //Read and check arguments  introduced by the user
    while ( (c = getopt(argc, argv, "s:ifcd:h:altg") ) != -1)  {
        switch (c)  {
   			//Type of analysis 
            case 's': //Single analysis 
                p.s_on = 1;
                p.pathFile=optarg;
                //Check the extension of the input file
                if (strcmp(strrchr(p.pathFile,'.'),".xyz")!=0 && strcmp(strrchr(p.pathFile,'.'),".xmol")!=0 ) display_usage(argv);
            //Type of molecular system
            break;
            case 'i': //Isolated peptide
                p.isol_on =1;
            break;
            case 'f': //CID (collision induced dissociation)
                p.cid_on=1;
            break;
            case 'c': //Clusters
                p.clst_on=1;
            break;
            //Parameters of analysis
            case 'd': //Bonds to take into account
                p.bonds_on=1;
                p.bonds_val = atoi(optarg);
                if(p.bonds_val<0 || p.bonds_val>15) display_usage(argv);
            break;
            // case 'r': //Rotational motion
            //     p.rot_on=true;
            // break;
            case 'h': //Angle in H-bonds
                if(atoi(optarg)!=0 && atoi(optarg)!=1)display_usage(argv);
                if(atoi(optarg)==1) p.ang_HB_on=true; 
                if(atoi(optarg)==0) p.ang_HB_on=false; 
            break;
            case 'a': //Extraction of energy from the input file
                // p.isom_anal= true;
                p.eng_anal=true;
            break;
            case 'l': //Partial analysis of the system (focus on the most important fragments) 
                p.path_on=1;
            break;
            case 't': //Get the statistical analysis (cycles, connected components, etc.) for the interface 
                p.total_on=1;
            break;
            case 'g': //Get the fragments analysis (analyse the list of fragments that have been observed along the trajectory, g: grain)
                p.frg_on=1;
            break;

            //Other argument    
        	case '?':
            	display_usage(argv);
            break;
        }
    }

//At least one module (type of analysis) should be introduced
    if (p.s_on != 1) display_usage(argv);

//If a single or multiple analysis are performed, check the other options
    if(p.s_on==1){
        if(p.isol_on+p.cid_on+p.clst_on+p.bonds_on==0) display_usage(argv);
    }

    return p;
}

/*========
Set the different options and threshold values in the global variable \a molecule according to \a opt
========*/
void get_anal_opt(struct MyModel* molecule,struct MyOpt opt ){
//Get directory and trajectory file names
    char dirName[512];
    sprintf(dirName,"%s",opt.pathFile);
    molecule->partAnal = false;
    molecule->fragAnal = false;
    if(strrchr(dirName,'/')!=NULL){
        char *fileName=str_sub(strrchr(dirName, '/'),1,strlen(strrchr(dirName, '/')));
        sprintf(molecule->inputFN,"%s", fileName);
        dirName[strlen(dirName)-strlen(strrchr(dirName,'/'))]='\0';
        sprintf(molecule->inputDir,"%s",dirName);
        free(fileName);
    }
    else{ //Use the current directory '.'
        sprintf(molecule->inputFN,"%s",opt.pathFile);
        sprintf(molecule->inputDir,"%s",".");
    }   
    if(opt.s_on){
        //Some treatments for a single analysis
        molecule->num_F=0; 
        molecule->nbInputF=1; //Only one trajectory is analysed         

        //The system type to which the trajectory belongs to : 'i' (isolated peptide), 'c' (cluster), 'f' (collision leading to fragmentation), etc.
        if(opt.isol_on+opt.cid_on+opt.clst_on==0 || opt.isol_on==1){
            molecule->sysType[0]='i';
            molecule->level=1; //H-bonds dynamics
        }
        else{
            if(opt.cid_on==1){
                molecule->sysType[0]='f';
                molecule->level=2; //Covalent bonds dynamics
            }
            else{
                molecule->sysType[0]='c'; 
                if(opt.path_on){ // localize the conformational dynamics to the most important fragments(partial analysis) 
                    molecule->partAnal = true;
                }
                molecule->level=13;//Electrostatic interactions + hydrogen bonds + organometallic interactions        
            }
        }
        //Bonds introduced by user
        if(opt.bonds_on==1){
            molecule->level=opt.bonds_val;
        }
        //Angle for the H-bond condition
        if(opt.ang_HB_on)
            molecule->anghbond='1'; 
        else
            molecule->anghbond='0';
        //Analysis of fragments
        if(opt.frg_on)
            molecule->fragAnal=true;
    }
    else{
        molecule->sysType[0]='i'; 
        molecule->anghbond='1';
        molecule->level=1; 
        molecule->isomAnal=false;
    }
//Initialize the reference snapshot which is the first snapshot in the trajectory. 
    molecule->num_imgRef = 0;
    molecule->num_img = 0;
//#Conformations
    molecule->nbConf=0;
//#Intermediate states 
    molecule->nbTransConf=0;    
//Conformations list : this is related to the graph of transitions
    molecule->conformations=NULL; 
    molecule->lastConf=NULL;  
//Fragments list : this is related to the list of fragments
    molecule->fragments=NULL;
//#Covalent bonds
    molecule->nbCov=0;
    molecule->nbCovF=0;
//#H-bonds found
    molecule->nbHbond=0;
//Initialize the type of changes, from right to left: bit_0 for change in H-bonds, bit_1 for change in the covalent bonds and bit_2 for change in electrostatic interactions.
    molecule->changType=0;
//Initialize of #fragments
    molecule->nbFrg=1;
//CM coordinates
    molecule->cm[0]=0.00;
    molecule->cm[1]=0.00;
    molecule->cm[2]=0.00;
//Threshold values for bonds
    get_parameter(molecule);
}

/*========
Get the threshold values used to calculate the different types of bonds (covalent, hydrogen bonds and electrostatic interactions)
========*/
void get_parameter(struct MyModel *molecule){
    

    if(true){//See after to let it a condition
        //Covalent bonds
        molecule->tshVal.covrC=COVRC;// Covalent radius of a carbon atom (Angstroms)
        molecule->tshVal.covrH=COVRH;// Covalent radius of a Hydrogen atom (Angstroms) 
        molecule->tshVal.covrN=COVRN;// Covalent radius of a Nitrogen atom (Angstroms) 
        molecule->tshVal.covrO=COVRO;// Covalent radius of a Oxygen atom (Angstroms) 
        molecule->tshVal.covrSi=COVRSi;// Covalent radius of a Silica atom (Angstroms) 
        molecule->tshVal.covrS=COVRS;// Covalent radius of a Sulfur atom (Angstroms) 
        molecule->tshVal.covrMn=COVRMn;// Covalent radius of a Maganese atom (Angstroms) 
        molecule->tshVal.covrB=COVRB;// Covalent radius of a Boron atom (Angstroms) 
        molecule->tshVal.covrP=COVRP;// Covalent radius of a Phosphorus atom (Angstroms) 
        molecule->tshVal.covrF=COVRF;// Covalent radius of a Flourine atom (Angstroms) 
        molecule->tshVal.covrRu=COVRRu;// Covalent radius of a Ruthenium atom (Angstroms) 
        molecule->tshVal.covrCl=COVRCl;// Covalent radius of a Chlorine atom (Angstroms) 
        molecule->tshVal.covrZn=COVRZn;// Covalent radius of a Zinc atom (Angstroms) 
        molecule->tshVal.covMarg=COVMARG; // Margin for distance used before considering the broken of covalent bond (%) 
        //Hydrogen bonds
        molecule->tshVal.distAHMax=DIST_AH_MAX; //Distance maximum between the hydrogen atom and the acceptor to form hydrogen bond (Angstrom)
        molecule->tshVal.angDHAMin=ANG_DH_MIN; //Angle minimum between the heavy atom, the hydrogen atom and the acceptor atom  to form hydrogen bond (degree) 
        molecule->tshVal.prtMarg= PTR_MARG;//Percentage used to consider proton transfer (%)             
        molecule->tshVal.distDAMax=DIST_DA_MAX;//!<Distance maximum between the heavy atom (donor) and the acceptor to form hydrogen bond (Angstrom)
        molecule->tshVal.distDHAMax=DIST_DHA_MAX;//!<Distance maximum between the heavy atom (donor), hydrogen atom and the acceptor to potentially form hydrogen bond (Angstrom)
        //electrostatic interactions 
        molecule->tshVal.distLiAt=DIST_Li_ATOM; // Distance used for a interaction between Lithium atom and the other atoms (Angstroms)   
        molecule->tshVal.distLiAr=DIST_Li_Ar; // Distance used for interaction between Lithium and Argon (Angstrom) 
        molecule->tshVal.distArAt=DIST_Ar_ATOM; // Distance used for interaction between Argon and the other atoms (Angstrom) 
        molecule->tshVal.distKAt=DIST_K_ATOM; // Distance used for interaction between Potassium and the other atoms (Angstrom) 
        molecule->tshVal.distNaO=DIST_Na_O; // Distance used for interaction between Soduium and Oxygen (Angstrom) 
        molecule->tshVal.distNaN=DIST_Na_N; // Distance used for interaction between Sodium and Nitrogen Angstrom) 
        molecule->tshVal.distClCl=DIST_Cl_Cl; // Distance used for interaction between Chlorine and Chlorine Angstrom) 
        molecule->tshVal.distClO=DIST_Cl_O; // Distance used for interaction between Chlorine and Oxygen Angstrom) 
        molecule->tshVal.distBrO=DIST_Br_O; // Distance used for interaction between Brom and Oxygen Angstrom) 
        molecule->tshVal.distIO=DIST_I_O; // Distance used for interaction between Ionid and Oxygen Angstrom) 
        molecule->tshVal.distKO=DIST_K_O; // Distance used for interaction between Potassium and Oxygen Angstrom) 
        //organometallic interactions
        molecule->tshVal.distMnO=DIST_Mn_O; // Distance used for interaction between Maganese atom and oxygen atom (Angstroms)   
        molecule->tshVal.distMnH=DIST_Mn_H; // Distance used for interaction between Maganese and hydrogen atom (Angstrom) 
        molecule->tshVal.distMnBr=DIST_Mn_Br; // Distance used for interaction between Maganese and Brom atom (Angstrom) 
        molecule->tshVal.distMnF=DIST_Mn_F; // Distance used for interaction between Maganese and Flourine atom (Angstrom) 
        molecule->tshVal.distMnS=DIST_Mn_S; // Distance used for interaction between Maganese and Sulfer atom (Angstrom) 
        molecule->tshVal.distMnP=DIST_Mn_P; // Distance used for interaction between Maganese and Phosphorus atom (Angstrom) 
        molecule->tshVal.distRuCl=DIST_Ru_Cl; // Distance used for interaction between Ruthenium and Chlore atom (Angstrom) 
        molecule->tshVal.distRuF=DIST_Ru_F; // Distance used for interaction between Ruthenium and Flourine atom (Angstrom) 
        molecule->tshVal.distRuP=DIST_Ru_P; // Distance used for interaction between Ruthenium and Phosphorus atom (Angstrom) 
        molecule->tshVal.distAuAu=DIST_Au_Au; // Distance used for interaction between Gold and Gold atom (Angstrom) 
        //Conformations
        molecule->tshVal.pourtMinResT=POURT_MIN_RES_T; // Time residence to consider a conformation 
    }
    else{
        FILE *trashValF=fopen("constants.txt","r");
        //Covalent bonds
        if(fscanf(trashValF,"%lf", &molecule->tshVal.covrC)!=1){
            printf("can not read the covalent radius of carbone\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.covrH)!=1){
            printf("can not read the covalent radius of hydrogen\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.covrN)!=1){
            printf("can not read the covalent radius of nitrogen\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.covrO)!=1){
            printf("can not read the covalent radius of oxygen\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.covrS)!=1){
            printf("can not read the covalent radius of sulphur\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.covrMn)!=1){
            printf("can not read the covalent radius of maganese\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.covrB)!=1){
            printf("can not read the covalent radius of boron\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.covMarg)!=1){
            printf("can not read the margin percentage used for covalent bonds\n");
            exit(-1);
        }
        //Hydrogen bonds
        if(fscanf(trashValF,"%lf", &molecule->tshVal.distAHMax)!=1){
            printf("can not read the distance max used for hydrogen bonds\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.angDHAMin)!=1){
            printf("can not read the angle min used for hydrogen bonds\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.prtMarg)!=1){
            printf("can not read the percentage used for proton transfer\n");
            exit(-1);
        }
        //Electrostatic interactions
        if(fscanf(trashValF,"%lf", &molecule->tshVal.distLiAt)!=1){
            printf("can not read the distance max used between Lithium and an atom\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.distLiAr)!=1){
            printf("can not read the distance max used between Lithium and an Argon\n");
            exit(-1);
        }
        if(fscanf(trashValF,"%lf", &molecule->tshVal.distArAt)!=1){
            printf("can not read the distance max used between Argon and an atom\n");
            exit(-1);
        }        
        if(fscanf(trashValF,"%lf", &molecule->tshVal.distKAt)!=1){
            printf("can not read the distance max used between Potassium and an atom\n");
            exit(-1);
        }        
        if(fscanf(trashValF,"%lf", &molecule->tshVal.pourtMinResT)!=1){
            printf("can not read the minimum percentage of time residence of a conformation\n");
            exit(-1);
        }        
        fclose(trashValF);      
    }
}

