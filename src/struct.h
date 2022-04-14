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

/*! \file struct.h

    \brief This file contains all structures used by different parts of the proposed algorithms
*/

/*====================Structures============================*/

/*! \struct MyOpt
    \brief This class contains the program's options (all user's input parameters)
*/
struct MyOpt{
	char *pathFile;
	int s_on; //!< Single analysis
	int m_on; //!< Multiple analysis
	int p_on; //!< Conformational prediction 
	int w_on; //!< Interface analysis 
	int isol_on; //!< MD trajectory of an isolated peptide in gas phase
	int cid_on; //!<  MD trajectory a collision leading to the fragmentation of a peptide in gas phase
	int clst_on; //!< MD trajectory of an isolated cluster in gas phase
	int bonds_on; //!< Bonds to take into account for the conformational change
	int bonds_val; //!< Value of bonds to take into account
	bool rot_on; //!< Rotational motion around covalent bonds
	bool ang_HB_on; //!< Angle criterion to analyse the hydrogen bonds
	bool isom_anal; //!< Isomers analysis (used for the CID analysis)
	int eng_on; //!< Computation of energy
	int eng_coef; //!< Coefficient of H-bonds energy 
	int bar_on; //!< Construct the graph of possible conformations
	int bar_val; //!< Barrier of energy
	bool eng_anal; //!< Real energy average  (compute energy from the xyz file)
	int path_on; //!< Get the low cost path between two conformations  
 	int pbc_on; //!< Take into account the periodic boundary condition (PBC) 
    int total_on;//!< Add additional statistical analysis for the water interfaces 
    int frg_on; //!< Add additional analysis, the analysis of fragments 
    double pbcx_val; //!< The x size of the box for the PBC 
    double pbcy_val; //!< The y size of the box for the PBC 
    double pbcz_val; //!< The z size of the box for the PBC 	
};

/*! \struct threshVal
    \brief This class contains threshold values and constants used to compute different bonds 
*/
struct threshVal{
	//Covalent bonds
	double covrC; //!< Covalent radius of a Carbon atom (Angstrom)
	double covrH; //!< Covalent radius of a Hydrogen atom (Angstrom) 
	double covrN; //!< Covalent radius of a Nitrogen atom (Angstrom) 
	double covrO; //!< Covalent radius of a Oxygen atom (Angstrom) 
	double covrSi; //!< Covalent radius of a Silica atom (Angstrom) 
	double covrS; //!< Covalent radius of a sulfur atom (Angstrom) 
	double covrMn; //!< Covalent radius of a Manganese atom (Angstrom) 
	double covrB; //!< Covalent radius of a Boron atom (Angstrom) 
	double covrF; //!< Covalent radius of a Flourine atom (Angstrom) 
	double covrP; //!< Covalent radius of a Phosphorus atom (Angstrom) 
	double covrRu; //!< Covalent radius of a Ruthenium atom (Angstrom) 
	double covrCl; //!< Covalent radius of a Chlorine atom (Angstrom) 
	double covrZn; //!< Covalent radius of a Zinc atom (Angstrom) 
	double covMarg; //!< Margin for distance used before considering the broken of covalent bond (%) 
	//Hydrogen bonds
	double distAHMax; //!< The maximum distance between the hydrogen atom and the acceptor to form hydrogen bond (Angstroms)
	double distDAMax; //!< The maximum distance between the heavy atom (donor) and the acceptor to form hydrogen bond (Angstrom)
	double distDHAMax; //!< The maximum distance between the heavy atom (donor), hydrogen atom and the acceptor to potentially form hydrogen bond (Angstrom)
	double angDHAMin; //!< The minimum angle between the heavy atom, the hydrogen atom and the acceptor atom to form hydrogen bond (degree) 
	double prtMarg; //!< Percentage used to consider proton transfer (%)			 
	
	//Intermolecular bonds
	double distLiAt; //!< Distance used for an interaction between Lithium atom and other types of atoms (Angstrom) 	
	double distLiAr; //!< Distance used for an interaction between Lithium and Argon (Angstrom) 
	double distArAt; //!< Distance used for an interaction between Argon and other types of atoms (Angstrom) 
	double distKAt; //!< Distance used for an interaction between Potassium and other types of atoms (Angstrom) 
	double distFH; //!< Distance used for an interaction between Flourine and hydrogen (Angstrom) 
	double distNaN; //!< Distance used for an interaction between Sodium and Nitrogen (Angstrom) 
	double distNaO; //!< Distance used for an interaction between Sodium and Oxygen (Angstrom) 
	double distIO; //!< Distance used for an interaction between Iodine and Oxygen (Angstrom) 
	double distClO; //!< Distance used for an interaction between Chlorine and Oxygen (Angstrom) 
	double distClCl; //!< Distance used for an interaction between Chlorine and Chlorine (Angstrom) 
	double distBrO; //!< Distance used for an interaction between Brom and Oxygen (Angstrom) 
	double distKO; //!< Distance used for an interaction between Potassium and Oxygen (Angstrom) 
	
	//Organometallic bonds
	double distMnO; //!< Distance used for an interaction between Maganese atom and Oxygen (Angstrom)
	double distMnH; //!< Distance used for an interaction between Maganese atom and Hydrogen (Angstrom)
	double distMnBr; //!< Distance used for an interaction between Maganese atom and Brom (Angstrom)
	double distMnF; //!< Distance used for an interaction between Maganese atom and Flourine (Angstrom)
	double distMnS; //!< Distance used for an interaction between Maganese atom and Sulfer (Angstrom)
	double distMnP; //!< Distance used for an interaction between Maganese atom and Phosphorus (Angstrom)
	double distRuF; //!< Distance used for an interaction between Ruthenium atom and Flourine (Angstrom)
	double distRuP; //!< Distance used for an interaction between Ruthenium atom and Phosphorus (Angstrom)
	double distRuCl; //!< Distance used for an interaction between Ruthenium atom and Chlorine (Angstrom)
	double distAuAu; //!< Distance used for an interaction between Gold atom and Gold (Angstrom)

	//Conformation or intermediate state
	double pourtMinResT; //!< Minimum time residence percentage to consider a conformation as stable and not an intermediate state

	//Size of interface boxes (periodic boundary conditions - PBC)
	double boxX; //!< The x box size of the periodic boundary conditions (Angstrom)
	double boxY; //!< The y box size of the periodic boundary conditions (Angstrom)
	double boxZ; //!< The z box size of the periodic boundary conditions (Angstrom)

};

/*! \struct atom
    \brief This class describes the characteristics of an atom
*/
struct atom{
	char atomType; //!< One character that defines the atom type : 'C' (Carbon),'O' (Oxygen),'N' (Nitrogen), 'H' (Hydrogen) 'A' (Argon), 'L' for lithium etc.
	char atomName[3]; //!< Full chemical name of atom : 'C' (Carbon), 'Ar'(Argon), etc.
	double x ; //!< The x Cartesian coordinate
	double y ; //!< The y Cartesian coordinate
	double z ; //!< The z Cartesian coordinate
	int bondMax; //!< The maximum number of covalent bonds that an atom could form  
	bool atomIn; //< A boolean used as an indicator of the atom is taken into account or not to analyse one conformation (used for the clusters)
};

/*! \struct bondH
    \brief This class contains a hydrogen bond (HB) characteristic formed by one hydrogen atom of the molecular system 
*/
struct bondH{ 
	int numH; //!< Index of the Hydrogen atom
	int numDon; //!< Index of the donor atom (heavy atom)
	int numAcc; //!< Index of the acceptor atom
	char sens; //!< One character that defines sense of the HB '0' (without proton transfer) or '1'(with proton transfer)
	int state; //!< Get the state of the HB : '+1' HB present, '-1' HB absent, '0'  potential HB
	int nbOrb; //!< Number of atoms in the orbit of hydrogen atom (numH)
	int numOrb[NB_MaxOrb]; //!< List of atom indexes of hydrogen orbit
};

/*! \struct trajHBond
    \brief This class contains a HB dynamics (appearance/disappearance occurred along the trajectory)
*/
struct trajHBond{
	int bondH[3]; //!< The HB label Donor-Hydrogen-Acceptor
	int nChange; //!< Number of changes (appearance/disappearance) that have been occurred along the trajectory
	int imgList[NB_PERIOD*3]; //!< Periods of each appearance of the HB: for each period, three boxes of the table are used: 1st box for the HB sense (with/out proton transfer), 2nd box for instant where the HB appears, and a 3rd box for instant where it disappears
};

struct successor;

/*! \struct confList
    \brief This class contains all information related to one conformation 

    Each conformation represents a mixed graph where atoms are the vertices of the graph and bonds (covalent, hydrogen or electrostatics bonds) between them represent edges/arcs of the graph. The bonds are stocked as adjacency matrices.
	For each conformation, successors (neighbours) are also saved.   
	To that end, a graph of transitions is generated, where each conformation is a vertex in the graph and arcs (directed edges) represent the transitions explored along the dynamics 
*/
struct confList 
{
	int name ; //!< Label of the conformation, an unique number is assigned to each conformation 
	struct atom* img ; //!< Pointer to the list of atoms of the conformation (molecular system)
	int nbAtom; //!< Number of atoms in the conformation
	int imgList[3]; //!< Periods of each appearance of the conformation: for each period, two boxes are used: 1st box for instant where the conformation appears, and a 2nd box for instant it disappears
	int **HB; //!< Adjacency matrix of HB identified for this conformation
	int **CB; //!< Adjacency matrix of covalent bonds identified for this conformation, with additional information about rotational motion
	int **IB; //!< Adjacency matrix of intermolecular bonds identified for this conformation 
	int **MB; //!< Adjacency matrix of organometallic bonds identified for this conformation 
	int nbCov; //!< Number of covalent bonds
	int nbIon; //!< Number of ions in this conformation 
	int nbMetal; //!< Number of metal atoms in this conformation
	int nbFrg ; //!< Number of fragments in this conformation
	double totalPerc; //!< Total percentage of appearance
	double engMin; //!< The minimum energy value found for the conformer  
	int snapMin; //!< The snapshot/step related to the minimum energy vlue found for the conformer 
	char state; //!< One character that defines if a conformation is considered as stable or as an intermediate state
	struct successor* succ; //!< Pointer to the list of successors (neighbours) of this conformation 
	struct confList *suiv; //!<  Pointer to next conformation 
};

/*! \struct successor
    \brief This class contains characteristics of a conformation successor (transition)
*/
struct successor{ 
	struct confList* conf; //!< Pointer to a conformation (vertex)
	int freq;  //!< Frequency of going from a conformation to \a conf (can be seen as a edge weight in the graph of transitions)
	int transf; //!< Seven bits number, that gives details about changes occurring, from right to left: bit_0 for appearance of new HB, bit_1 for disappearance of HB, bit_2 for proton transfer, bit_3 for appearance of new covalent bond, bit_4 for disappearance of covalent bond, bit_5 for appearance of new intermolecular bond, bit_6 for disappearance of intermolecular bond.
	struct successor* suiv; //!< Pointer to next successor
};

/*! \struct frgList
    \brief This class contains all information related to one fragment 

    Each fragment represents a mixed graph where atoms are the vertices of the graph and bonds (covalent, hydrogen or electrostatics bonds) between them represent edges/arcs of the graph. The bonds are stocked as adjacency matrices.
*/
struct frgList {
	int name ; //!< Label of the fragment, an unique number is assigned to each fragment 
	struct atom* img ; //!< Pointer to the list of atoms of the conformation (molecular system)
	int nbAtom; //!< The number of atoms in the fragment (B matrix)
	int size;//!< Number of atoms in the fragment (the size)
	int **B; //!< Adjacency matrix of bonds involved in the fragment
	int nbocc; //!< Number of occurrence of the fragment
	int *conftab; //!< List of conformers where the fragment appears 
	struct frgList *suiv; //!<  Pointer to next conformation 
};

/*! \struct path_elt
    \brief This class contains one element of a queue (FIFO) or stack (LIFO)
*/
struct path_elt {
	int index; //!< Label of a node (element)
	struct path_elt* suiv; //!< Pointer to the next element of the queue or stack
};

/*! \struct path
    \brief This class contains the pointers of a queue or stack
*/
struct path  {
	struct path_elt *head; //!< Pointer to the first element inserted in the queue (or last element inserted for a stack)
	struct path_elt *queue; //!< Pointer to the last element inserted in the queue (or first element inserted for a stack)
};

/*! \struct MyModel
    \brief This is a global class that contains all variables used in the program 
	
	Note that a trajectory is a series of snapshots taken at a regular time scale. Hereafter the current snapshot means the current step in the trajectory analysed
*/
struct MyModel{
	char inputDir[256];	//!< Input directory name
	char **inputFNList; //!< Input files name list (multiple trajectories analysis)
	char inputFN[256]; //!< Input file name (one trajectory analysis)
	char periodFile[256];
	int nbInputF; //!< Number of input files
	int num_F; //!< The label of the current file
	char name[256]; //!< Empirical formula of the molecular system analysed
	char sysType[1]; //!< One character that defines which type of system the trajectory belongs to : 'i' (isolated peptide in gas phase), 'c' (cluster in gas phase), 'f' (collision leading to fragmentation (CID)), etc.
	bool partAnal; //!> A boolean variable, true if the user wants to concentrate the analysis on specific fragments (containing the metal, cation or anion atoms)
	bool fragAnal; //!> A boolean variable, true if the user wants to analyse the different analysis
	bool nbAtomChange; //!< A boolean variable. At a snapshot, this parameter is true if the atom number has been changed in comparison with the previous snapshot
	char anghbond; //!< One character used to consider or not the angle in the hydrogen bond conditions
	int level; //!< 4 bits number that let user set bond types (covalent, hydrogen or ionic bonds) would he take into account for the conformational dynamics analysis
	bool isomAnal; //!< A boolean variable, true if user would launch isomers analysis (for CID trajectories only covalent bonds are taken into account)
	struct threshVal tshVal; //!< A vector used to read and save threshold values, constants used to compute different bonds. 
	unsigned int changType; //!< 3 bits number used for changes occurring in the current snapshot, from right to left: bit_0 for change in H-bonds, bit_1 for change in the covalent bonds and bit_2 for change in intermolecular bonds.
	double max1; //!< The maximum displacement of atoms at the current snapshot according to a reference snapshot
	double max2; //!< The second maximum displacement of atoms at the current snapshot according to a reference snapshot
	double cm[3]; //!< Cartesian coordinates of center of mass of the molecular system
	int num_img; //!< The current snapshot label
	struct atom* img ; //!< Pointer to the list of atoms of the current snapshot
	int num_imgRef; //!< The reference snapshot label
	struct atom* imgRef ; //!< Pointer to the list of atoms of the reference snapshot	
	int nbAtom; //!< Size of a snapshot (number of atoms in the molecular system)
	double avrgAtom; //!< Average of #atoms on all trajectory (for trajectories where number of atoms changes from one step to another)
	int nbMolecule; //!< The trajectory size (number of snapshots)
	int nbCov; //!< Number of covalent bonds identified in the current snapshot  
	int nbCovF; //!< Number of covalent bonds lists found along all the trajectory 
	int nbHbond; //!< Number of HBs found along all the trajectory
	int nbIon; //!< Number of ion's atoms (cation, anion) in the molecular system 
	int nbMetal; //!< Number of ion's atoms (cation, anion) in the molecular system 
	int nbHAtom; //!< Get the number of hydrogen atoms on the molecular system
	struct bondH bondHList[NB_MaxHBond]; //!< List of HBs identified at the current snapshot 
	struct trajHBond bondHdyn[NB_MaxHBond]; //!< Hbonds dynamics analysis
	int **donHacc; //!< Dynamics of Donor-acceptor
	int **covBond; //!< List of covalent bonds (as an adjacency matrix)
	int **interBond; //!< List of internal covalent bonds (as an adjacency matrix)
	int **ionBond; //!< List of electrostatic interactions  (bonds with cation/anion , as an adjacency matrix)
	int **metalBond; //!< List of electrostatic interactions  (bonds with cation/anion , as an adjacency matrix)
	int **bondDyn; //!< List of bonds/interactions that have been changed at least one time along the trajectory : bondDyn[i][j]= 1 [covBond], 2 [HBond], 3 [interBond], 4 [metalBond]
	int nbAtomDyn; //!< Number of atoms for the bondDyn matrix
	int isthList[NB_MaxIsthm]; //!< List of isthmuses (bridges on the mixed graphs which are considered as a rotational axes)
	int nbConf; //!< Number of conformations found on the trajectory
	int nbTransConf; //!< Number of intermediate states (conformations that appears less than \a POURT_MIN_RES_T of time)
	struct confList* conformations; //!< Head of list of conformations explored along the trajectory
	struct confList* lastConf; //!< Tail list of conformations explored along the trajectory
	struct frgList* fragments; //!< Head of list of fragments explored along the trajectory
	int nbFrg ; //!< Number of fragments (especially used for CID trajectories or clusters)	
	int cp[NB_MAXPAIR*2]; //!< List of couples (heavy atom, acceptor) that could potentially form hydrogen bonds
	int nbpairs; //!< Number of Potential HB identified
	int dist[NB_MAXPAIR]; //!< List of euclidean distances between  couples (heavy atom, acceptor)
	int nbH[NB_MAXPAIR]; //!< Number of HBs that one heavy atom can potentially form (i.e. number of hydrogen atoms linked to the heavy atom)
	int OrbitHSize; //!< For performance , the average size of orbits around hydrogen atoms
	int possHB[16]; //!< Just for Ylene (2019)
};

/*! \struct confListIso
    \brief This class contains brief description of one conformation/isomer.
*/
struct confListIso  {
	int confName; //!< The label of the conformations
	int isoName; //!< The label of the isomer
	int nbBonds; //!< Number of covalent bonds
	int nbfrg; //!< Number of fragments
	int tabBonds[2*NB_MaxAtom]; //!< Covalent bonds list (AtomA-AtomB) , sorted by atom index
};

/*! \struct ModelIso
    \brief This class contains the list of conformations/isomers found for series of trajectories (conformations search).
*/
struct ModelIso{
  	char inputDir[500]; //!< Input directory name (the full path) 
  	char **inputconFNList; //!< Input conformation list 	
  	struct confListIso *confIsoList; //!< Conformations/isomers list
  	int nbConfF; //!< Number of input files for conformations
  	int nbIso; //!< Number of input files for isomers
};

/*! \struct MyGraph 
    \brief This class contains relationship between conformations/isomers found for  of trajectories (conformational dynamics).
*/
struct MyGraph{
  	char inputDir[500]; //!< Input directory name (the full path) 
  	char inputFN[256]; // !< The input file name
  	int nbsnapshots; //!< Get the size of the trajectory (number of snapshots)	
	int **adjMatrix; //!< Contains all transitions between conformations/isomers according all trajectories (global matrix) or for one trajectory
	int nbdynm; //!< #conformations or #isomers found for all the trajectories analysed
	char dynType; //!<  Type of the analysis performed : 'c' for conformational analysis and 'i' for isomer analysis
	int fstate; //!< first state : the first conformation that has been explored in the trajectory
	int lstate; //!< last state : the last conformation that has been explored in the trajectory
};

/*! \struct MyGraphList
    \brief This class is used to create a list of graphs of transitions between the conformations/isomers that have been identified in the different trajectories analysed.
*/
struct MyGraphList {
	int name; //!< A label for the graph of transitions 
	struct MyGraph* stateG; //!< Pointer to a MyGraph (description of the graph of transition) element
	struct MyGraphList* suiv; //!< Pointer to the next element of the list
};

/*! \struct MyTrajGraph
    \brief This class contains a description of the graph of transitions for each trajectory analysed 
*/
struct MyTrajGraph{
	char **files; //!< A Char vectors which contain the input files name (trajectories) that have been analysed 
  	int nbInputF; //!< Number of input files (trajectories) that have been analysed
	int *filegraph; //!< The filegraph[i] contains the label to the graph of transitions assigned to files[i]
	int *fState; //!< fState[i] contains the initial state for the trajectory files[i]
	int *lState; //!< lState[i] contains the final state for the trajectory files[i]
  	struct MyGraphList* graphList; //!< Pointer to the head of the graphs of transitions list
  	int nbgraph;
};


/*=========Structure used for the conformational prediction or interfaces analyses============*/

struct nodDijk;

/*! \struct conf_elt
    \brief This class contains a description of one conformation in terms of energy proprieties. 

    Note that a conformation is defined as n character where conf[i] is 
    the i_th hydrogen bond that can be formed in the molecular system. and conf[i] takes one of the two values '1' (H-bond present) or '0' (H-bond absent) 
*/
struct conf_elt {
	char *conf; //!< A label for a conformation
	int index; //!< An index for a conformation in a table (PS : it is also used to contains the number of Hydrogen bonds in this conformation)
	int eng; //!< The energy of the conformation
	int pen; //!< The penalty of conformation
	struct nodDijk* d; //!< Pointer to the list of successor or conformations already visited
	struct conf_elt* suiv; //!< Pointer to the next element of list
};


/*! \struct  confIndex
	\brief  This class describes a node used  to index the conformations
*/
struct confIndex{
	int index; //!< Index of conformation in confTab
	struct confIndex* suiv; //!< Pointer to the next element 
};

/*! \struct nodCC
	\brief  This class describes one element of the connected components list ((node of \a tabCC))
*/
struct nodCC{
	int size; // The size of the connected component
	struct confIndex* confH; //!< Pointer to the head of the conformations list
	struct confIndex* confQ; //!< Pointer to the queue of the conformations list

};

/*! \struct  nodHash
	\brief   This class contains a description of one element in the hash-age table (node of \a tabHash) which is used to accelerate computation   
*/
struct nodHash{
	struct conf_elt* confList; //!< Pointer to the conformations list 
};


/*! \struct nodEng 
	\brief This class describes one element of the conformations energies list  
*/
struct nodEng{
	int val; //!< Energy value
	int nbConf; //!< Number of conformations That have an energy value equal to \a val
	struct conf_elt* confList; //!< Pointer to the list of conformations That have an energy value equal to \a val
	struct nodEng* suiv; //!< Pointer to the next element of the list 
};

/*! \struct  nodConf
	\brief  This class describes a conformation in terms of connected components ((node of \a tabConf))

	Note that a conformation is defined as n character where conf[i] is 
    the i_th hydrogen bond that can be formed in the molecular system. and conf[i] takes one of the two values '1' (H-bond present) or '0' (H-bond absent) 
*/
struct nodConf{
	char *conf; //!< A label of the conformation  
	int nbNeigh; //!< Number of neighbours of conformation \a conf
	int cc; //!< The label of the connected component to which \a conf belongs 
};

/*! \struct GraphModel
    \brief This class contains the description of the graph of possibles conformations

    In the conformational prediction, a list of conformations are generated according to a set of local and global rules based on the change of hydrogen bonds. A theoric energy (indicator) value is
    assigned to each conformation. A graph of transitions is constructed. It is a weighted directed graphs such that conformations represent the vertices of the
    graph and arcs represent transition between them. Weights are assigned to the vertices, i.e. each conformation has a theoric energy (indicator). 
    According to the energy of the system, the graph of transitions changes and can be more or less connected (number of connected components increases or decreases)
*/
struct GraphModel{
	char inputDir[512]; //!< Input directory name
	struct nodEng* listEng;	 //!< Pointer to the list of energy values with conformations
	struct nodConf* tabConf; //!< Pointer to the list of conformations with their connected component labels 
	struct nodCC* tabCC; //!< Pointer to the list of connected components labels with their size and conformations
	struct nodHash* tabHash; //!<  Pointer to list of indexes, indexing process for the \a tabConf list using hash function
	int nbConfT; //!< Total number of conformations in the graph
	int nbValEng; //!< Number of energy values that are different 
	int minEng; //!< Minimum of energy value found
	int maxEng; //!< Maximum of energy value found
	int maxCC; //!< Maximum number of connected component
	int nbCC; //!< Number of connected components found in the graph (DFS algorithm)
	int ccLabel; //!< Last label used for connected components, this is used if the connected components of the graph should be merged.  
	int nbEdge; //!< Number of graph edges 
	int nbpairs; //!< Number of couples  (heavy atom, acceptor) that can form a hydrogen bonds in the molecular system 
};

/*! \struct weights
    \brief This class contains a description of the arc weights on the graph of possible conformations

    This class is used in order to evaluate a path cost between two conformations. In our case we use Dijkstra Algorithm
*/
typedef struct weights {
	double cost; //!< An arc is characterized by the cost of the change that has been occurred  between two conformations
	int dist; //!< Path size between two conformations in terms of number of arcs visited.
} weight;


/*! \struct nodDijk 
    \brief This class contains a description of a successor node or a node that has been visited and weights to arrive to this node

    This class is used in Dijkstra Algorithm to find paths with lower costs between conformations
*/
struct nodDijk{
	struct conf_elt *confP; //!< Pointer to conformation element
	double *tabPen; //!< Table of atoms penalties for conformation \a confP
	weight w; //!< Weight of going from conformation  pred->confP to conformation \a confP
	int visited; //!< An binary value which indicates if conformation \a confP has already been visited or not
	struct nodDijk *pred; //!< Pointer to the predecessor of conformation \confP
	struct nodDijk *prev; //!< Pointer to previous element in the list
	struct nodDijk *suiv; //!< Pointer to next element in the list
};

/*! \struct nodHB
    \brief This class contains the list of atoms belonging to one hydrogen bond and hence the size of the hydrogen bond in terms of the number of atoms
*/
struct nodHB{
	int *tabAtoms;
	int size;
};
/*! \struct DijkModel
    \brief This class contains fields used to browse paths between 2 conformations

    This class is used in Dijkstra Algorithm to find a set of paths of lower cost between two conformations.
    The Dijkstra algorithm constraints are defined in terms of path length (number of conformation browsed)
	and the maximum of changes costs that have been occurred in the browsed path  
*/
struct DijkModel{
	char outputDir[512]; // !< Output directory name
	struct nodDijk* visConfList; //!< Pointer to the list of conformations (nodes in the global graph) already visited 
	struct nodDijk* neighConfList; //!< Pointer to the list of neighbouring conformations (nodes that have not been visited yet)
	struct nodHash* tabHash; //!< Hash table to accelerate conformations search/access
	struct nodHB* tabHBAtoms; //!< A table to save lists of atoms belonging to each hydrogen bond 
	int nbHBonds; //!< Total number of potential hydrogen bonds that can be formed in the conformation
	int nbConf; //!< Number of conformations in the global graph 
	double maxEng; //! A barrier to limit paths taken into account between two conformations (only paths with a cost < \a maxEng are considered)
};

/*! \struct cycle
    \brief This class contains description of a cycle in a conformation (mixed graph).

    Cycles are obtained used Horton algorithm. 
*/
struct cycle{
	int *tabCycle;  //!< Pointer to the list of atoms (vertices) belonging to the cycle
	int *tabBonds;  //!< Pointer to the list of bonds (edges) of the cycle : tabBonds[i] = 1 if bondList[i] is present (see MyCycleGraph struct), 0 else
	int nbHB; //!< Number of hydrogen bonds in the cycle
	int size;  //!< Size of the cycle (number of edges)
	int atomIgnored; // !< A parameter used to choose a specific atoms (only used for interfaces)
	struct cycle *suiv;  //!< Pointer to the next element 
};

/*! \struct MyCycleGraph
    \brief This class contains cycles list that have been found using Horton algorithm
*/
struct MyCycleGraph{
	struct cycle *cyclesList; //!< Pointer to list of cycles 
	int nbCycles; //!< Number of cycles in the list \a cyclesList
	int *bondList; //!< Index of atoms for each covalent or hydrogen bond
	int nbBonds; //!< Total number of bonds in \a bondList
};

/*! \struct GaussModel
    \brief This class contains description of Gaussian model used to get minimum basis of cycles
	
	After application of Horton algorithm, a minimum basis can be extracted by Gaussian elimination
*/ 
struct GaussModel{
	int *tabBonds; //!< Gaussian vector,  tabBonds[i] = 1 if bondList[i] is present (see \a MyCycleGraph struct ), 0 else  
	struct cycle *cycl; //!< Pointer to a cycle in \a cyclesList
	int nbBonds; //!< Number of bonds belonging to the cycle
};


/*! \struct CCModel
    \brief This class contains description of a connected component for a graph.

    This class is mainly used to analyse the connectivity of molecules and ions in interfaces 
*/
struct CCModel{
	int *tabCC; //!< Statistic vector, tabCC[i] is the number of occurrence of the connected component of size i for the whole trajectory
	int *tabCCPour ;//!< Statistic vector, tabCCPour[i] is the pourcentage of occurrence of the connected component of size i for the whole trajectory (10%, 20%,30%...) 
	int bigCC; //!<  Size of the biggest connected component in the mixed graph in the whole trajectory 
	int smallCC; //!<  Size of the smallest connected component in the mixed graph at each 
	int *bigCCtab; //!< Pointer to the list of atoms belonging to the biggest connected component in the mixed graph at each step
	int nbCC; //!< Number of connected components found in the whole trajectory
	double averageCC; //!<  Average of the number of connected components in graphs for the whole trajectory 
	double averageSize; //!< Average of size of connected component found for the whole trajectory
};


//This structure is only used in the prediction part, for file processing. The "line" structure is used to read a line from .data file
struct line{
	char conf[256] ; // Vector of the conformation
	int label; // Label of the conformation
	int nb_hb; // Number of hydrogen bonds
	int pen; // Penalty of the conformation
	int eng; // Energy of the conformation
};