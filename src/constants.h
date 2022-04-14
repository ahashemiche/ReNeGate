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

/*! \file constants.h
    \brief This file contains all constants and macros used in the program.
*/

/*========
Macros
========*/

/*! \def Max(a,b)
    \brief A macro that returns the maximum of \a a and \a b.
*/
#define Max(a,b)((a)>(b)?(a):(b))

/*! \def Min(a,b)
    \brief A macro that returns the minimum of \a a and \a b.
*/
#define Min(a,b)((a)<(b)?(a):(b))

/*! \def bit_1(n,k)
    \brief A macro that returns TRUE if the bit \a k of \a n equal to 1.
*/
#define bit_1(n,k)((n)&(1<<(k))?1:0)

/*! \def bit_0(n,k)
    \brief A macro that returns TRUE if the bit \a k of \a n equal to 0.
*/
#define bit_0(n,k)((n)&(1<<(k))?0:1)
/*! \def add_1(n,k)
    \brief A macro that set 1 to the bit \a k of \a n.
*/
#define add_1(n,k)((n)|(1<<(k)))

/*! \def add_0(n,k)
    \brief A macro that set 0 to the bit \a k of \a n.
*/
#define add_0(n,k)((n)&(~(1<<(k))))

/*========
Constants
========*/

/*! \def NB_MaxAtom
    \brief Defines the maximum number of atoms by molecule.
*/
#define NB_MaxAtom 100 

/*! \def NB_MaxAtomType
    \brief Defines the maximum number of types of atoms : N, H, O, C, Li, Cl, S, Mn, Ar, B, K, Br, Si, P, F, Ru, Na, I, Au, Zn 
*/
#define NB_MaxAtomType 20

/*! \def NB_MaxHBond
    \brief Defines the maximum #H-bonds.
*/
#define NB_MaxHBond 5760 //2885 //600 //700

/*! \def NB_MaxInterBond
    \brief Defines the maximum #ionic-bonds.
*/
#define NB_MaxInterBond 200 

/*! \def NB_MaxOrb
    \brief Defines the maximum #atoms on the orbit of an hydrogen.
*/
#define NB_MaxOrb 300 //170 

/*! \def NB_MaxIsthm
    \brief Defines the maximum #isthmuses.
*/
#define NB_MaxIsthm 100 

/*! \def NB_PERIOD
    \brief Defines the maximum of periods for a hydrogen bond.
*/
#define NB_PERIOD 10 //850

/*! \def MAX_CYCLE_SIZE 
    \brief Defines the maximum size of a cycle in the interfaces trajectories  
*/
#define MAX_CYCLE_SIZE  50 

/*! \def ALPHA
    \brief Coefficient of optimization, used to change the reference snapshot and hence orbits
*/
#define ALPHA 1.00  

/*! \def NB_MAXPAIR
    \brief Defines maximum number of pairs generated.
*/
#define NB_MAXPAIR 200

/*! \def MAX_ENG
    \brief Parameter used only to fix the minimum and the maximum of energy value.
*/
#define MAX_ENG 100000  

/*! \def MAX_CONF
    \brief Maximum number of conformation to not relabel the connected component.
*/
#define MAX_CONF 1000 

/*! \def HASH_TAB_SIZE
    \brief Maximum of hash index in tabHash (the size of entries in the hash table)
*/
//#define HASH_TAB_SIZE 100000 

/*! \def MAX_PACK
    \brief Coefficient used for compression hash function.
*/
#define MAX_PACK 18  

/*! \def JUMP_SNAPS
    \brief Coefficient used to analyse just sub set of snapshots each JUM_SNAPS 
*/
#define JUMP_SNAPS 1900  

/*! \def PI
    \brief The mathematical constant Pi.
*/
#define PI 3.14159265358979323846 

/*! \def EPSILON
    \brief A constant to compare float/double values.
*/
#define EPSILON 0.001 

/*! \var typedef enum bool
    \brief A type definition for boolean variables.
*/
typedef enum { false, true } bool; 

/*========
Thrashold values 
========*/

/*! \var DIST_AH_MAX
    \brief Distance maximum between the hydrogen atom and the acceptor to form hydrogen bond (Angstroms)
*/
#define DIST_AH_MAX  2.3  //2.30

/*! \var DIST_DA_MAX
    \brief  Distance maximum between the heavy atom and the acceptor to form hydrogen bond (Angstroms)
*/
#define DIST_DA_MAX 3.2 //6 //3.2 

/*! \var ANG_DH_MIN
    \brief Angle minimum between the heavy atom, the hydrogen atom and the acceptor atom  to form hydrogen bond (degree) 
*/
#define ANG_DH_MIN 120 //100 //120.00

/*! \var DIST_DHA_MAX
    \brief  Distance maximum between the heavy atom, hydrogen atom and the acceptor to form hydrogen bond (Angstroms)
*/
#define DIST_DHA_MAX 3.5 

/*! \var PTR_MARG 
    \brief Percentage used to consider proton transfer (%)             
*/
#define PTR_MARG 0.7

/*! \var COVRC
    \brief Covalent radius of a carbon atom (Angstroms)
*/
#define COVRC 0.76
 
/*! \var COVRH
    \brief Covalent radius of a Hydrogen atom (Angstroms)
*/
#define COVRH 0.31 // to change to 0.35
  
/*! \var COVRN
    \brief Covalent radius of a Nitrogen atom (Angstroms)
*/
#define COVRN 0.66
  
/*! \var COVRO
    \brief Covalent radius of a Oxygen atom (Angstroms)
*/
#define COVRO 0.71
  
/*! \var COVRSi
    \brief Covalent radius of a silica atom (Angstroms)
*/
#define COVRSi 1.16

/*! \var COVRMn
    \brief Covalent radius of a Manganese atom (Angstroms)
*/
#define COVRMn 1.35
 
/*! \var COVRS
    \brief Covalent radius of a Sulfur atom (Angstroms)
*/
#define COVRS 1.02
 
/*! \var COVRB
    \brief Covalent radius of a Boron atom (Angstroms)
*/
#define COVRB 0.83

/*! \var COVRP
    \brief Covalent radius of a Phosphorus atom (Angstroms)
*/
#define COVRP 1.05  

/*! \var COVRCl
    \brief Covalent radius< of a Chlorine atom (Angstroms)
*/
#define COVRCl 0.71  


/*! \var COVRF
    \brief Covalent radius of a Florine atom (Angstroms)
*/
#define COVRF 0.64  

/*! \var COVRRu
    \brief Covalent radius of a Ruthenium atom (Angstroms)
*/
#define COVRRu 1.40  

/*! \var COVRZn
    \brief Covalent radius of a Zinc atom (Angstroms)
*/
#define COVRZn 1.45  

/*! \var COVMARG
    \brief Margin for distance used before considering the broken of covalent bond (%)
*/
#define COVMARG 1.3
  
/*! \var DIST_Li_ATOM
    \brief Distance used for a interaction between Lithium atom and the other atoms (Angstroms)   
*/

#define DIST_Li_ATOM 3.5 //1.82 + 1.20 (H) , 1.70 (C)
 
/*! \var DIST_Li_Ar
    \brief Distance used for interaction between Lithium and Argon (Angstroms)
*/
#define DIST_Li_Ar 3.60 // 1.82 + 1.88
  
/*! \var DIST_Ar_ATOM
    \brief Distance used for interaction between Argon and the other atoms (Angstroms) 
*/
#define DIST_Ar_ATOM 3.76 // 1.88 + 1.88

/*! \var DIST_F_H
    \brief Distance used for interaction between Flourine and hydrogen (Angstroms) 
*/
#define DIST_F_H  1.23 //0.95 //1.23 //  0.64  + 0.31

/*! \var DIST_K_ATOM
    \brief Distance used for a interaction between Potassium atom and the other atoms (Angstroms)   
*/
#define DIST_K_ATOM 3.0 

/*! \var DIST_K_O
    \brief Distance used for a interaction between Potassium atom and Oxygen atom (Angstroms)   
*/
#define DIST_K_O 3.54

/*! \var DIST_Cl_O 
    \brief Distance used for a interaction between Chlorine atom and Oxygen atom (Angstroms)   
*/
#define DIST_Cl_O 3.71

/*! \var DIST_Cl_Cl
    \brief Distance used for a interaction between Chlorine atom and Chlorine atom (Angstroms)   
*/
#define DIST_Cl_Cl 4.79 

/*! \var DIST_Br_O 
    \brief Distance used for a interaction between Brom atom and Oxygen atom (Angstroms)   
*/
#define DIST_Br_O 3.71

/*! \var DIST_I_O 
    \brief Distance used for a interaction between Iodine  atom and Oxygen atom (Angstroms)   
*/
#define DIST_I_O 4.2 

/*! \var DIST_Na_O
    \brief Distance used for a interaction between Sodium atom and Oxygen atom (Angstroms)   
*/
#define DIST_Na_O 2.5 

/*! \var DIST_Na_N
    \brief Distance used for a interaction between Sodium atom and the Nitrogen atom (Angstroms)   
*/
#define DIST_Na_N 2.4 

/*! \var DIST_Mn_H
    \brief Distance used for interaction between Manganese and the other hydrogen (Angstroms) 
*/
#define DIST_Mn_H 2.04

/*! \var DIST_Mn_F
    \brief Distance used for interaction between Manganese and the other Flourine (Angstroms) 
*/
#define DIST_Mn_F  2.50 // 0.64  + 1.35 = 1,94

/*! \var DIST_Mn_O
    \brief Distance used for interaction between Manganese and the other oxygen (Angstroms) 
*/
#define DIST_Mn_O 2.44

/*! \var DIST_Mn_S
    \brief Distance used for interaction between Manganese and the other Sulfer (Angstroms) 
*/
#define DIST_Mn_S  2.95 // 0.64  + 1.35 = 1,94

/*! \var DIST_Mn_P
    \brief Distance used for interaction between Manganese and the other Phosphorus (Angstroms) 
*/
#define DIST_Mn_P 2.97

/*! \var DIST_Mn_Br
    \brief Distance used for interaction between Manganese and  Brome  (Angstroms) 
*/
#define DIST_Mn_Br 2.60

/*! \var DIST_Ru_Cl
    \brief Distance used for interaction between Ruthenium  and Chlorine (Angstroms) 
*/
#define DIST_Ru_Cl 2.40


/*! \var DIST_Ru_F
    \brief Distance used for interaction between Ruthenium  and Florine (Angstroms) 
*/
#define DIST_Ru_F 2.05

/*! \var DIST_Ru_P 
    \brief Distance used for interaction between Ruthenium  and Phosphorus (Angstroms) 
*/
#define DIST_Ru_P  2.45

/*! \var DIST_Au_Au 
    \brief Distance used for interaction between Gold  and Gold (Angstroms) 
*/
#define DIST_Au_Au  3.2

/*! \var POURT_MIN_RES_T
    \brief Threshold Percentage to consider a conformation as stable else it will be a intermediate state 
*/
#define POURT_MIN_RES_T 4.00

/*! \var CC_SIZ_MIN
    \brief The minimum size to consider a connected component in the interfaces 
*/
#define CC_SIZ_MIN 1

/*! \var CYCLE_SIZE_MAX
    \brief The minimum size to consider a connected component in the interfaces 
*/
#define CYCLE_SIZE_MAX 10

/*! \var MIN_SIZE_CC_DIST
    \brief The minimum size to consider a connected component in the distribution  
*/
#define MIN_SIZE_CC_DIST 1

/*! \var MAX_SIZE_CC
    \brief The minimum size to consider a connected component in the distribution  
*/
#define MAX_SIZE_CC 300

/*========
Periodic Boundary Conditions (PBC)
========*/

/*! \var BOX_X
    \brief The X dimension of the box for the PBC (Angstroms)
*/
// #define BOX_X  13.386  //Water-Air interface wanlin 
#define BOX_X 19.734 //Water-Air interface 
// #define BOX_X 29.601 //Water-Air2 interface
// #define BOX_X 14.214 //Water-Alumina interface
// #define BOX_X 21.491 //Water-Graphene interface
// #define BOX_X 21.85 //Water-BN interface
// #define BOX_X 12.670 //Water-Silica interface

/*! \var BOX_Y
    \brief The Y dimension of the box for the PBC (Angstroms)
*/
// #define BOX_Y  13.386  //Water-Air interface wanlin 
#define BOX_Y 19.734 //Water-Air interface
// #define BOX_Y 29.601 //Water-Air2 interface
// #define BOX_Y 16.478 //Water-Alumina interface
// #define BOX_Y 22.224 //Water-Graphene interface
// #define BOX_Y 22.709 //Water-BN interface
// #define BOX_Y 13.270 //Water-Silica interface

/*! \var BOX_Z
    \brief The Z dimension of the box for the PBC (Angstroms)
*/
// #define BOX_Z  85.00  //Water-Air interface wanlin 
#define BOX_Z 35.00 //Water-Air interface
// #define BOX_Z 19.734 //Water-Air interface
// #define BOX_Z 45.00 //Water-Air2 interface
// #define BOX_Z 35.00 //Water-Alumina interface
// #define BOX_Z 35.00 //Water-Graphene interface
// #define BOX_Z 35.00 //Water-BN interface
// #define BOX_Y 37.00 //Water-Silica interface

