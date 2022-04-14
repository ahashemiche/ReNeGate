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

/*! \file comfunct.h
    \brief This file contains all headers related to common functions used by many parts of the algorithm
*/

#ifndef COMFUNCT_H
#define COMFUNCT_H 

/*! \fn skip_return (FILE *F, int n);
 *  \brief Skip n carriage return (\n) in a file.
 *  \param F Pointer to a file.
 *  \param n The number of return (\n) to skip in the file F.
 */
void skip_return (FILE *F, int n);

/*! \fn  skip_snapshot(FILE* traj, int num_snap, int nb_snaps);
 *  \brief Skip JUMP_SNAPS snapshots from num_snap
 *  \param F Pointer to a file.
 *  \param num_snap The start point.
 *  \param nb_snaps The tolal number of snapshots in the file traj    
 */
void skip_snapshot(FILE* traj, int num_snap, int nb_snaps);

/*! \fn create_res_folder(char path[500], char folder_name[256])
 *  \brief Create folder_name directory in path.
 *  \param path The full path where the folder is created.
 *	\param folder_name The directory name to be created.
 */
void create_res_folder(char path[500], char folder_name[256]);

/*! \fn get_FN(char pathFN[], char inputFN[],char* ext)
 *  \brief Create a file name from inputFN name, with 'ext' extension.
 *  \param pathFN The directory path to put the file created. 
 *  \param inputFN The file name used as prefix of the new one. 
 *	\param ext The new extension that replaces the extension of inputFN, in the new file name generated.
 *  \return The file name generated.
 */
char* get_FN(char pathFN[], char inputFN[],char* ext);

/*! \fn get_nbFiles(char filPath[500])
 *  \brief Get the number of files located in filePath.
 *  \param filPath The path of the location of files to be considered. 
 *  \return The number of files located in filePath.
 */
int get_nbFiles(char filPath[500]);

/*! \fn get_nbCoordFiles(char filPath[500]) 
 *  \brief Get the number of trajectories to be analysed. Only files with an extension ".xyz" or ".xmol" are taken into account.
 *  \param filPath The path of the location of files to be considered. 
 *  \return The number of files located in filePath with extensions ".xyz" and ".xmol".
 */
int get_nbCoordFiles(char filPath[500]);

/*! \fn get_FileList(struct MyModel *molecule)
 *  \brief Get the list of the trajectories files to be analysed.
 *  \param molecule A pointer to a structure of type MyModel. The list of files is saved in molecule->inputFNList list.
 */
void get_FileList(struct MyModel *molecule);

/*! \fn del_file(char *fileName)
 *  \brief Delete the file fileName if it exists.
 *  \param fileName The file name to be deleted.
 */
 void del_file(char *fileName);

/*! \fn del_dir(char *dirName)
 *  \brief Delete the directory dirName if it exists, with all its content.
 *  \param dirName The directory name to be deleted.
 */
 void del_dir(char *dirName);

/*! \fn snapshot_size(FILE* traj)
 *  \brief Compute the size of the first snapshot in traj file. We suppose that this size consists in the number of atoms plus 2 lines of comments.
 *  \param traj A pointer to a file.
 *  \return The size of the first snapshot in traj file.
 */
int snapshot_size(FILE* traj);

/*! \fn skip_snapshot(FILE* traj, int num_snap, int size)
 *  \brief Skip num_snap snapshots in trajectory  traj. 
 *  \param traj A pointer to a file.
 *	\param num_snap The number of snapshots to be skipped. 
 *	\param size The size of one snapshot. We suppose that all snapshots have the same size.
 */
// void skip_snapshot(FILE* traj, int num_snap, int size);

/*! \fn format_str(char c, int index)
 *  \brief Viewing format (AtomNumber).
 *  \param c The name of an atom.
 *	\param index The index of the atom.
 *  \return The name of the atom concatenated to its index. 
 */
char* format_str(char c, int i);

/*! \fn str_sub (const char *s, unsigned int start, unsigned int end)
 *  \brief Get a substring from string from start to end (inclusive). 
 *  \param s A string 
 *	\param start The position of the first character of substring
 *	\param end The position of the last character of substring
 *  \return The complete substring, from the first character to the last.
 */
char *str_sub (const char *s, unsigned int start, unsigned int end);

/*! \fn put_char(char *word, char c)
 *  \brief Replace the character in position index by the character c 
 *  \param word A string
 *	\param index A position in word 
 *	\param c A character to be used 
 *  \return A new string by replacing the the character in position index by the character c 
 */
char* put_char(char *word, int index, char c);

/*! \fn copy_matrix(int **matrix1,int **matrix2,int size)
 *  \brief Make a copy of the  matrix matrix1 into the matrix matrix2.
 *  \param[in] matrix1 A pointer to a 2 dimensional matrix of integers.
 *	\param[out] matrix2 A pointer to a 2 dimensional matrix of integers.
 *	\param size The size of the matrix (square matrix).
 */
void copy_matrix(int **matrix1,int **matrix2,int size);

/*! \fn copy_tab(int *tab1,int *tab2,int size)
 *  \brief Make a copy of the  vector tab1 into the vector tab2.
 *  \param[in] tab1 A pointer to a vector of integers.
 *	\param[out] tab2 A pointer to a vector of integers.
 *	\param size The size of the vector.
 */
void copy_tab(int *tab1,int *tab2,int size);

/*! \fn  init_tab(int *tab,int size, int val)
 *  \brief Initialise the vector tab of integers with the value val.
 *  \param[in, out] tab A point to a vector of integers.
 *	\param size The size of the vector.
 *	\param val An integer used to initialize the vector.
 */
void init_tab(int *tab,int size, int val);

/*! \fn  init_matrix(int **matrix,int sizeL, int sizeC, int val)
 *  \brief Initialise the matrix matrix of integers with the value val.
 *  \param matrix A pointer to a 2 dimensional matrix of integers.
 *  \param sizeL The number of rows in the matrix.
 *  \param sizeC The number of columns in the matrix.
 *	\param val An integer used to initialize the matrix.
 */
void init_matrix(int **matrix,int sizeL, int sizeC, int val);

/*! \fn  init_tab_double(double *tab,int size, double val)
 *  \brief Initialise the vector tab of doubles with the value val.
 *  \param[in, out] tab A point to a vector of doubles.
 *	\param size The size of the vector.
 *	\param val A double (float) used to initialize the vector.
 */
void init_tab_double(double *tab,int size, double val);

/*! \fn  init_matrix_double(double **matrix,int sizeL, int sizeC, double val)
 *  \brief Initialise the matrix matrix of doubles with the value val.
 *  \param matrix A pointer to a 2 dimensional matrix of doubles.
 *  \param sizeL The number of rows in the matrix.
 *  \param sizeC The number of columns in the matrix.
 *	\param val A double used to initialize the matrix.
 */
void init_matrix_double(double **matrix,int sizeL, int sizeC, double val);

/*! \fn distance (double xi,double yi,double zi,double xj,double yj,double zj)
 *  \brief Calculate the euclidean distance between two points (atoms).
 *  \param xi The x coordinate of the point i. 
 *	\param yi The y coordinate of the point i. 
 *	\param zi The z coordinate of the point i. 
 *  \param xj The x coordinate of the point j. 
 *	\param yj The y coordinate of the point j. 
 *	\param zj The z coordinate of the point j. 
 *  \return The euclidean distance between two points (atoms).
 */
double distance (double xi,double yi,double zi,double xj,double yj,double zj);

/*! \fn distance_pbc (struct atom atomA, struct atom atomB, double boxX,double boxY,double boxZ)
 *  \brief Calculate the euclidean distance between two atoms considering the Periodic Boundary Conditions (PBC).
 *  \param atomA A pointer to structure of type atom.
 *	\param atomB A pointer to structure of type atom.
 *	\param boxX The x coordinate of the box.
 *	\param boxY The y coordinate of the box.
 *	\param boxZ The z coordinate of the box.
 *  \return The euclidean distance between two points (atoms) with respect to the PBC.
 */
double distance_pbc (struct atom atomA, struct atom atomB, double boxX,double boxY,double boxZ);

/*! \fn angle (double xi,double yi,double zi,double xj,double yj,double zj,double xk,double yk,double zk)
 *  \brief Calculate the angle between three points (atoms).
 *  \param xi The x coordinate of the point i. 
 *	\param yi The y coordinate of the point i. 
 *	\param zi The z coordinate of the point i. 
 *  \param xj The x coordinate of the point j. 
 *	\param yj The y coordinate of the point j. 
 *	\param zj The z coordinate of the point j. 
 *  \param xk The x coordinate of the point k. 
 *	\param yk The y coordinate of the point k. 
 *	\param zk The z coordinate of the point k. 
 *  \return The angle between three points (atoms).
 */
double angle (double xi,double yi,double zi,double xj,double yj,double zj,double xk,double yk,double zk);

/*! \fn angle_pbc(struct atom atomD, struct atom atomH,struct atom atomA, double boxX,double boxY,double boxZ)
 *  \brief Calculate the angle between three atoms considering the Periodic Boundary Conditions (PBC).
 *  \param atomD A pointer to structure of type atom.
 *	\param atomH A pointer to structure of type atom.
 *	\param atomA A pointer to structure of type atom.
 *	\param boxX The x coordinate of the box.
 *	\param boxY The y coordinate of the box.
 *	\param boxZ The z coordinate of the box.
 *  \return The angle between three points (atoms) with respect to the PBC.
 */
double angle_pbc(struct atom atomD, struct atom atomH,struct atom atomA, double boxX,double boxY,double boxZ);

/*! \fn get_cc(int** CC, int size)
 *  \brief Get the number of connected components in a graph using the Breadth-First Search algorithm.
 *  \param CC The adjacency matrix of the graph.
 *	\param size The size of the adjacency matrix. 
 *	\param save A boolean, if it is set to "true", then the fragments are saved
 *  \return The number of connected components in a graph.
 */
int get_cc(int** CC, int size,bool save, int **frg);

/*! \fn init_adj_conf_matrix(int** CC, int nbAtom,struct confList* conf, int level,  char type)
 *  \brief Initialize the adjacency matrix with the different bonds according to the level used.
 *  \param [in,out]CC The global adjacency matrix of conformation conf.
 *  \param nbAtom The number of atoms in the molecular system.
 *	\param conf A pointer to a structure of type confList.
 *	\param level The level of comparison, a 4 bits number that indicates which bonds are taken into account in the comparison.  
 */
void init_adj_conf_matrix(int** CC, int nbAtom,struct confList* conf, int level,  char type);


/*! \fn get_cc1(int** CC, int size, char outputFN[],int *atomIgnored,struct CCModel* ccDistrib)
 *  \brief Get the number of connected components in a graph using the Breadth-First Search algorithm. 
 *	In this function more treatments are performed about the distribution of the identified connected components (cc).
 *  \param CC The adjacency matrix of the graph.
 *	\param size The size of the adjacency matrix. 
 *	\param[in] outputFN The name of the output file.  
 *	\param[in,out] atomIgnored the number of water molecules not taking into account
 *	\param[out] ccDistrib A pointer to a structure of type CCModel. This contains all the characteristics of the cc identified in the graph
 *  \param[in] nbHAtom number of atoms  
 *  \return The number of connected components in a graph.
 */
int get_cc1(int** CC, int size, char outputFN[],int *atomIgnored,struct CCModel* ccDistrib/*,int nbHAtom*/);

/*! \fn path* add_elt(struct path* F , int x)
 *  \brief Add element to de queue F (add at the end).
 *  \param F A pointer to a structure of type path (queue).
 *	\param x An integer to add in the queue.
 *  \return A pointer to the queue F with new element x.
 */
struct path* add_elt(struct path* F , int x);

/*! \fn get_elt(struct path* F)
 *  \brief Get the first element of the queue and delete this element from it.
 *  \param F A pointer to a structure of type path (queue).
 *  \return Return a pointer to the queue F without the first element.
 */
//int get_elt(struct path* F);

/*! \fn get_elt_notMarked(int** CC, int size, int index , int color[])
 *  \brief Get a covalent bond formed with the vertex index and which is not marked (not reached yet)
 *  \param CC The adjacency matrix of the graph.
 *	\param[in] outputFN The name of the output file.  
 *	\param index The index of the atom concerned. 
 *	\param colour A vector of colours of the vertices of the graph (atoms). 
 *  \return The position of the vertex (atom) that has a covalent bond with vertex index
 */
int get_elt_notMarked(int** CC, int size, int index , int color[]);

/*! \fn num_occ(struct atom* img, int atomA, int nbAtom)
 *  \brief Calculate the order of atom in the atomList, according to it type. 
 *  \param img A pointer to a structure of type atom. 
 *	\param atomA An index of atom. 
 *	\param nbAtom Number of atoms in the molecular system. 
 *  \return The order of the atom atomA in the atoms list according to it type.
 */
int num_occ(struct atom* img, int atomA, int nbAtom);

/*! \fn get_nb_frg(int **covBond,int size)
 *  \brief Get the number of fragments (connected components) according to the covalent bonds. 
 *  \param covBond A pointer to a 2 dimensional matrix of integers. The adjacency matrix of the covalent bonds.
 *	\param size The size of the covBond matrix.
 *  \return The number of fragments according to the covalent bonds.
 */
int get_nb_frg(int **covBond,int size);

/*! \fn  get_interface_frg(struct confList* conf,int size,char outputFN[], int *atomIgnored,struct CCModel* ccDistrib)
 *  \brief Get the number of fragments in conformation conf and save these fragments (connected components) in ccDistrib.
 *  \param conf A pointer to a structure of type confList.
 *	\param size The number of atoms in the molecular system.
 *	\param outputFN The output file name to save the size of each identified fragment.
 *	\param atomIgnored A pointer to a vector of integers. It contains the indexes of atoms that have been ignored.
 *	\param [out]ccDistrib A pointer to a structure of type CCModel. This contains all fragments identified in conformation conf.	
 *  \return The number of fragments in conformation conf.
 */
int get_interface_frg(struct confList* conf,int size,char outputFN[], int *atomIgnored,struct CCModel* ccDistrib);

/*! \fn get_nbTypeatoms(struct atom* img, char *atomName,int size)
 *  \brief Get the number of atoms with name atomName within the molecular system.
 *  \param img A pointer to a structure of type atom.
 *	\param atomName A chemical name of an atom (C for carbon, O for oxygen, Li for lithium, etc.)
 *  \param size Number of atoms in the molecular system
 *  \return The number of atoms with name atomName. 
 */
int get_nbTypeatoms(struct atom* img, char *atomName,int size);

/*! \fn get_index(struct MyModel *molecule,char *atomName, int occ)
 *  \brief Get index of the occurrence occ of an atom with name atomName.
 *  \param molecule A pointer to a structure of type MyModel.
 *	\param atomName A chemical name of an atom (C for carbon, O for oxygen, Li for lithium, etc.)
 *	\param occ An occurrence of an atom according to its type. 
 *  \return The index of the occurrence occ of an atom with name atomName.
 */
int get_index(struct MyModel *molecule,char *atomName, int occ);

/*! \fn get_index_ion(int **IB,int size1, int size2,int index)
 *  \brief Get the index of an ion in the IB matrix.
 *  \param IB A pointer to a 2 dimensional matrix of integers. The adjacency matrix of the electrostatic interactions. 
 *	\param size1 The number of rows in the IB matrix.
 *	\param size2 The number of columns  in the IB matrix.
 *	\param index The index of the ion in the atoms list. This index is already present in the IB matrix.
 *  \return The index of the ion in the IB matrix.
 */
int get_index_ion(int **IB,int size1, int size2,int index);

/*! \fn pair_valid(struct MyModel* molecule, int index1, int index2)
 *  \brief Check if atoms with indexes index1 and index2 can forms a hydrogen bond, according to a set of conditions. 
 *  \param molecule A pointer to a structure of type MyModel.
 *	\param index1 An index of an atom.
 *	\param index2 An index of an atom.
 *  \return True if the atoms index1 and index2 respect all conditions and can form a hydrogen bond. 
 */
bool pair_valid(struct MyModel* molecule, int index1, int index2);

/*! \fn atoms_near(struct MyModel* molecule, int index1, int index2)
 *  \brief Check if the atoms are at a distance less than 3 (in terms of separated bonds, i.e. edges in the graph)
 *  \param molecule A pointer to a structure of type MyModel.
 *	\param index1 An index of an atom.
 *	\param index2 An index of an atom.
 *  \return True if the atoms are at a distance less than 3. 
 */
bool atoms_near(struct MyModel* molecule, int index1, int index2);

/*! \fn exist_atom(struct MyModel* molecule, int index, char atomType)
 *  \brief Check if there is a covalent bond between index and atom with chemical type atomType.
 *  \param molecule A pointer to a structure of type MyModel. The function uses the molecule->covBond matrix.
 *	\param index An index of an atom.
 *	\param atomType A chemical type of an atom. 
 *  \return True if the atom with index index forms a covalent bonds with at least one atom with chemical type atomType.
 */
bool exist_atom(struct MyModel* molecule, int index, char atomType);

/*! \fn exist_H(struct MyModel* molecule, int index)
 *  \brief Check if there is a covalent bond between index and an hydrogen atom.
 *  \param molecule A pointer to a structure of type MyModel.
 *	\param index An index of an atom. 
 *  \return True if the atom with index forms a covalent bond with a hydrogen atom.
 */
bool exist_H(struct MyModel* molecule, int index);

/*! \fn get_H_num(struct MyModel* molecule, int index)
 *  \brief Get the number of hydrogen atoms related to the atom index.
 *  \param molecule A pointer to a structure of type MyModel.
 *	\param index An index of atom.
 *  \return The number of hydrogen atoms covalently bonded to the atom with index.
 */
int get_H_num(struct MyModel* molecule, int index);

/*! \fn get_Hbond_sens(int **HB,int size, int A, int D)
 *  \brief Find the orientation of (donor,acceptor). Simple H-bond or with proton transfer.
 *	\param HB A pointer to a 2 dimensional matrix of integers. The adjacency matrix of hydrogen bond.
 *  \param size The size of the HB matrix.
 *	\param A An index of an atom. 
 *	\param D An index of an atom.
 *  \return The 3 if there has been a proton transfer, 2 else.
 */
int get_Hbond_sens(int **HB,int size, int A, int D);

/*! \fn check_formed_Hbond(struct MyModel *molecule, int index)
 *  \brief Check if the H-bond is formed for a significant period of time along trajectory.
 *  \param molecule A pointer to a structure of type MyModel. This function uses molecule->bondHdyn field.
 *	\param index An index of a H-bond in molecule->bondHdyn.
 *  \return True if the H-bond has been appeared for a significant period of time.
 */
bool check_formed_Hbond(struct MyModel *molecule, int index);

/*! \fn verif_timeRes_conf(struct confList* conf,double pourt,int size)
 *  \brief Check if the conformation has been appeared at least a percentage of pourt along the trajectory. 
 *  \param conf A pointer to a structure of type confList.
 *	\param pourt The minimum percentage used to consider one conformation.  
 *	\param size The size of the trajectory.
 *  \return True if the conformation has been appeared at least a percentage of pourt along the trajectory. 
 */
bool verif_timeRes_conf(struct confList* conf,double pourt,int size);

/*! \fn get_donHB(int **donHacc,int index,int size)
 *  \brief Get the index of the donor (heavy atom) bonded to the hydrogen atom index.
 *  \param donHacc A pointer to a 2 dimensional matrix of integers.
 *	\param index An index of hydrogen atom. 
 *  \param size The size of the donHacc matrix.
 *  \return The index of the the donor (heavy atom) bonded to the hydrogen atom index.
 */
int get_donHB(int **donHacc,int index,int size);

/*! \fn get_accHB(int **donHacc,int index,int size)
 *  \brief Get the index of the acceptor from matrix of donor-hydrogen-acceptor (donHacc) that contains all H-bonds.
 *  \param  donHacc A pointer to a 2 dimensional matrix of integers.
 *	\param index An index of a hydrogen atom.
 *	\param size The size of the donHacc matrix.
 *  \return The index of an atom which forms a HB with a hydrogen index, -1 else. 
 */
 int get_accHB(int **donHacc,int index,int size);

/*! \fn get_HHB(int **donHacc,int index1, int index2, int size)
 *  \brief Check if the index of the donor is either index1 or index2.
 *  \param donHacc A pointer to a 2 dimensional matrix of integers.
 *	\param index1 An index of an atom.
 *	\param index2 An index of an atom. 
 *	\param size The size of the donHacc matrix.
 *  \return The index of the donor.
 */
int get_HHB(int **donHacc,int index1, int index2, int size);

/*! \fn check_ptr(struct confList *conf,int index1,int index2,int size, int *D)
 *  \brief Get the index of the donor and check if there has been a proton transfer in the conformation conf.
 *  \param conf A pointer to a structure of type confList.
 *	\param index1 An index of an atom.
 *	\param index2 An index of an atom. 
 *	\param size The number of atoms in the molecular system. This function uses conf->img list. 
 *	\param [out]D The index of the donor in the HB formed between index1 and index2.
 *  \return True if there has been a proton transfer.
 */
bool check_ptr(struct confList *conf,int index1,int index2,int size, int *D);

/*! \fn get_isom_num(struct ModelIso *isomList, int confName)
 *  \brief Get the isomer's name which corresponds to the conformation confName. This function is used in the isomerization analysis. 
 *  \param isomList A pointer to a structure of type isomList. 
 *	\param confName A name of a conformation.
 *  \return The name of the isomer which corresponds to the conformation confName.
 */
int get_isom_num(struct ModelIso *isomList, int confName);

/*! \fn verif_type_state(int *tabState,int size,int index)
 *  \brief Check if the conformation index is in tabState in at least one trajectory.
 *  \param tabState A pointer to a vector of integers. It contains conformations names. 
 *	\param size The size of the vector tabState.
 *	\param index An index of a conformation.
 *  \return The number of occurrence of conformation index in tabState.
 */
int verif_type_state(int *tabState,int size,int index);

/*! \fn get_HHHB(struct confList *conf,int index1, int index2, int size)
 *  \brief Get the index of the hydrogen atom involved in the HB formed between atoms index1 and index2.
 *  \param conf A pointer to a structure of type confList.
 *	\param index1 An index of an atom.
 *	\param index2 An index of an atom. 
 *	\param size The number of atoms in the molecular system. This function uses conf->img list. 
 *  \return The index of the hydrogen atom involved in the HB formed between atoms index1 and index2.
 */
int get_HHHB(struct confList *conf,int index1, int index2, int size);


/*! \fn draw_interface_distribution(char inputFN[], char outputFN[], int size, char *title, char *xlabel, char *ylabel, int type) 
 *  \brief Draw the distribution of the connected components (fragments) according to their size (number of atoms or water molecules).
 *  \param inputFN The input file name which contains the data to draw.
 *	\param outputFN The output file name.
 *	\param size The size maximum of a connected component (fragment)
 *	\param title The title of the plot. 
 *	\param xlabel The label of the x axis.
 *	\param ylabel The label of the y axis.
 *	\param type The type of drawing (points or histogram).
 */
void draw_interface_distribution(char inputFN[], char outputFN[], int size, char *title, char *xlabel, char *ylabel, int type) ;

/*! \fn get_snap_energy(FILE *traj,int size)
 *  \brief Read energy for one snapshot. We suppose that the energy appears in the second line after the number of atoms.
 *  \param traj The trajectory file.
 *	\param size The size of the snapshot.
 *  \return The energy related to the current snapshot of the molecular system.
 */
double get_snap_energy(FILE *traj,int size);

/*! \fn get_real_conf_energy(char *line)
 *  \brief Extract energy from a line in the xyz file
 *  \param line  A vector of characters. The content of a line in the xyz file.
 *  \return The value of energy contained in line.
 */
double get_real_conf_energy(char *line);

/*! \fn  get_atom_in_conf(struct MyModel *molecule,struct atom* img)
 *  \brief Get the list of atoms involved in the conformational search
 * 	\param molecule A pointer to a structure of type MyModel.
 *  \param img A pointer to a structure of type atom. 
 */
void get_atom_in_conf(struct MyModel *molecule,struct atom* img);

/*! \fn update_atom_in(int **CC, struct atom* img, int size)
 *  \brief Browse the fragments based on CC adjacency matrix and update the atomIn field in img structure 
 *  \param CC The adjacency matrix of the graph.
 *  \param img A pointer to a structure of type atom. 
 *	\param size The size of the snapshot.
 */
void update_atom_in(int **CC, struct atom* img, int size);

/*! \fn verif_dreadnaut_isom( char inputFN[])
 *  \brief Test isomorphism with Mckay (dreadnaut) algorithm. Check in the inputFN file if the two graphs are isomorphic or not.
 *  \param inputFN The name of the file containing the results of Mackay (nauty) algorithm.
 *  \return True if the graphs are isomorphic, false else.
 */
bool verif_dreadnaut_isom( char inputFN[]);

/*! \fn save_partition(FILE *outputF,struct atom *img,int nbAtom) 
 *  \brief  
 *  \param  
 *  \return 
 */
void save_partition(FILE *outputF,struct atom *img,int nbAtom);

/*! \fn get_cov_index (int **covbond,int size, int index) 
 *  \brief  
 *  \param  
 *  \return 
 */
int get_cov_index (int **covbond,int size, int index);


/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
void get_atom_in_interface(struct confList *conf);

/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
int get_nb_CB_type(struct atom *img,int index ,int **CB,int size,char *type);

/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
bool check_inv_HB(int **HB, int size, int index);

/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
void center_atoms(char *outputFN,struct atom* img, int size, double boxX, double boxY, double boxZ );


/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
double center_mass_axis(struct atom* img, int size, char axis);

/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
void get_stat(struct MyModel *molecule);

/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
int get_nb_atomIn(struct atom *img, int size);

/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
int get_nb_bonds(int **CC, struct atom *img, int size);

/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
int get_nb_CB_type(struct atom *img,int index ,int **CB,int size,char *type);


/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
int get_dir_path(int** CC, int size,  int **motifs);


/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */
int get_elt_notMarked2(int** CC, int size, int index , int color[]);


/*! \fn 
 *  \brief  
 *  \param  
 *  \return 
 */


#endif //COMFUNCT_H