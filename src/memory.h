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

/*! \file memory.h
    \brief This file contains all headers related to memory management.
*/

#ifndef MEM_H
#define MEM_H 

/*! \fn allocate_matrix( int sizeL,int sizeC, int init_val, char matrix_name[256]) 
 *  \brief Allocate memory for a 2 dimensional matrix of integer numbers with sizeL rows and sizeC columns.
 *  \param sizeL an integer to indicate number of rows in the matrix.
 *  \param sizeC an integer to indicate number of columns in the matrix.
 *	\param init_val an integer used to initialize the matrix.
 *	\param matrix_name a string to indicates label of the matrix.
 *  \return an integer pointer.
 */
int**  allocate_matrix( int sizeL,int sizeC, int init_val, char matrix_name[256]); 

/*! \fn free_matrix(int **matrix, int size) 
 *  \brief Release the memory allocated to a matrix of integer numbers.
 *	\param matrix A pointer to a 2 dimensional matrix.
 *  \param size an integer to indicate number of columns in the matrix.
 */
void free_matrix(int **matrix, int size);

/*! \fn allocate_char_matrix(int sizeL,int sizeC, char matrix_name[256]) 
 *  \brief Allocate memory for a 2 dimensional matrix of char with sizeL rows and sizeC columns.
 *  \param sizeL an integer to indicate number of rows in the matrix.
 *  \param sizeC an integer to indicate number of columns in the matrix.
 *	\param matrix_name a string to indicates label of the matrix.
 *  \return a character pointer.
 */
char**  allocate_char_matrix(int sizeL,int sizeC, char matrix_name[256]);

/*! \fn free_char_matrix(char **matrix, int size)
 *  \brief Release the memory allocated to a matrix of char.
 *	\param matrix A pointer to a 2 dimensional matrix.
 *  \param size an integer to indicate number of columns in the matrix.
 */
void free_char_matrix(char **matrix, int size);

/*! \fn  allocate_int_tab(int size,int init_val)
 *  \brief Allocate memory for a vector of type int 
 *  \param size an integer to indicate number of lines in the vector
 *	\param init_val a int value to initialize the elements of the vector
 *  \return a int pointer
 */
int* allocate_int_tab(int size,int init_val);

/*! \fn  void free_int_tab(int *tab)
 *  \brief Release the memory allocated to a vector of type int
 *  \param tab a vector of type int
 */
void free_int_tab(int *tab);

/*! \fn allocate_img(int size)
 *  \brief Allocate memory for a vector of structure atom (it represents one snapshot).
 *  \param size an integer to indicate number of elements of the vector.
 *  \return A pointer to a structure of type atom. 
 */
struct atom* allocate_img(int size);

/*! \fn free_img(struct atom* img)
 *  \brief Release the memory allocated to vector of type atom.
 *  \param img pointer to a vector of type atom.
 */
void free_img(struct atom* img);

/*! \fn allocate_confList( int sizeL, int sizeM, int sizeC,char* confList_name)
 *  \brief Allocate memory for a structure of type confList.
 *  \param sizeL an integer to indicate number of rows in the matrix.
 *  \param sizeM an integer to indicate number of rows in the matrix.
 *  \param sizeC an integer to indicate number of columns in the matrix.
 *	\param confList_name a string to indicates label of the matrix.
  *  \return A pointer to a structure of type confList.
 */
struct confList* allocate_confList( int sizeL, int sizeM, int sizeC,char* confList_name);

/*! \fn free_confList(struct confList* conf, int sizeL, int sizeM, int sizeC)
 *  \brief Release memory allocated for a struct confList.
 *	\param conf pointer to struct confList.
 *  \param sizeL an integer to indicate number of rows in the matrix.
 *  \param sizeM an integer to indicate number of rows in the matrix.
 *  \param sizeC an integer to indicate number of columns in the matrix.
 */
void free_confList(struct confList* conf, int sizeL, int sizeM, int sizeC);

/*! \fn allocate_frgList( int sizeL, int sizeM, int sizeC,char* frgList_name)
 *  \brief Allocate memory for a structure of type frgList.
 *  \param sizeL an integer to indicate number of rows in the matrix.
 *  \param sizeM an integer to indicate number of conformers in the vector.
 *  \param sizeC an integer to indicate number of columns in the matrix.
 *	\param frgList_name a string to indicates label of the matrix.
  *  \return A pointer to a structure of type frgList.
 */
struct frgList* allocate_frgList( int sizeL, int sizeM, int sizeC,char* frgList_name);

/*! \fn free_frgList(struct frgList* frg, int sizeL)
 *  \brief Release memory allocated for a struct frgList.
 *	\param conf pointer to struct frgList.
 *  \param sizeL an integer to indicate number of rows in the matrix.
 */
void free_frgList(struct frgList* frg, int sizeL);

/*! \fn allocate_matrix_model(struct MyModel *molecule,int size)
 *  \brief Allocate memory for a struct MyModel (global struct that contains all variables and matrices needed in the program).
 *	\param molecule A pointer to a struct MyModel.
 *  \param size an integer to indicate size of the matrices of the struct.
 */
void allocate_matrix_model(struct MyModel *molecule,int size);

/*! \fn free_matrix_model(struct MyModel *molecule,int size)
 *  \brief Allocate memory for a structure MyModel (global structure that contains all variables and matrices needed in the program).
 *	\param molecule A pointer to a structure MyModel.
 *  \param size an integer to indicate size of the matrices of the structure.
 */
void free_matrix_model(struct MyModel *molecule,int size);

/*! \fn allocate_tabHash()
 *  \brief Allocate memory for a vector of type nodHash. A hash table that stores conformation using hash function.
 *  \return A pointer to a vector of type nodHash.
 */
struct nodHash*  allocate_tabHash();

/*! \fn free_mymodel(struct MyModel *molecule)
 *  \brief Release memory allocated to a structure of type MyModel.
 *	\param molecule A pointer to a structure of type MyModel.
 */
void free_mymodel(struct MyModel *molecule);

/*! \fn free_graph(struct MyGraph *graph)
 *  \brief Release memory allocated to a structure of type MyGraph.
 *	\param graph A pointer to a structure of type MyGraph.
 */
void free_graph(struct MyGraph *graph);

/*! \fn free_graphList(struct MyTrajGraph *graphTraj)
 *  \brief Release memory allocated to a structure of type MyTrajGraph.
 *	\param graphTraj A pointer to a structure of type MyTrajGraph.
 */
void free_graphList(struct MyTrajGraph *graphTraj);

/*! \fn free_LIFO(struct path *F)
 *  \brief Get the first element of the list F (LIFO : last in first out).
 *	\param F A pointer to a LIFO.
 */
int get_elt_LIFO(struct path* F);

/*! \fn delet_elt_LIFO(struct path* F)
 *  \brief Delete the first element of the LIFO.
 *	\param F A pointer to a LIFO.
 */
struct path* delet_elt_LIFO(struct path* F);

/*! \fn free_LIFO(struct path *F)
 *  \brief Release memory allocated to a LIFO.
 *	\param F A pointer to a LIFO.
 */
void free_LIFO(struct path *F);

/*! \fn get_elt_FIFO(struct path* F)
 *  \brief Get the first element of the queue (FIFO).
 *	\param F A pointer to a FIFO.
*/
int get_elt_FIFO(struct path* F);

/*! \fn free_FIFO(struct path *F)
 *  \brief Release memory allocated to a FIFO.
 *	\param F A pointer to a FIFO.
 */
void free_FIFO(struct path *F);

/*! \fn free_tabHash(struct nodHash* tabHash, int size)
 *  \brief Release memory allocated to a hash table.
 *  \param tabHash A pointer to a vector of type nodHash.
 *	\param size The size of the vector tabHash.
 */
void free_tabHash(struct nodHash* tabHash, int size);

/*! \fn free_tabConf(struct nodConf* tabConf, int size)
 *  \brief Release the memory allocated to the vector tabConf.
 *  \param tabConf A pointer to a vector of type nodConf. 
 *	\param size The size of the vector tabConf.
 */
void free_tabConf(struct nodConf* tabConf, int size);

/*! \fn free_tabCC(struct nodCC* tabCC, int size)
 *  \brief Release the memory allocated to the vector tabCC.
 *	\param tabCC A pointer to a vector of type nodCC.
 *	\param size The size of the vector tabCC.
 */
void free_tabCC(struct nodCC* tabCC, int size);

/*! \fn free_listEng(struct nodEng *listEng)
 *  \brief Release the memory allocated to a structure of type listEng.
 *	\param listEng A pointer to a structure of type nodEng.
 */
void free_listEng(struct nodEng *listEng );

/*! \fn free_graphPoss(struct GraphModel* graphPoss)
 *  \brief Release the memory allocated to structure of type GraphModel.
 *  \param graphPoss A pointer to a structure of type GraphModel.
 */
void free_graphPoss(struct GraphModel* graphPoss);

/*! \fn free_level(struct conf_elt* level)
 *  \brief Release the memory allocated to a structure of type conf_elt.
 *  \param level A pointer to a structure of type conf_elt.
 */
void free_level(struct conf_elt* level);

/*! \fn free_path(struct DijkModel path)
 *  \brief Release the memory allocated to a structure of type DijkModel. A list of the shortest path using Dijkstra algorithm. 
 *  \param path An element of type DijkModel.
 */
void free_path(struct DijkModel path);

/*! \fn free_nodDijk(struct nodDijk *list)
 *  \brief Release the memory allocated to a structure of type nodDijk. 
 *  \param list A pointer to a list of elements of type nodDijk.
 */
void free_nodDijk(struct nodDijk *list);

/*! \fn free_tabHBAtoms(struct nodHB* tabHBAtoms,int size)
 *  \brief Release the memory allocated for a vector of type nodHB.
 *  \param tabHBAtoms A pointer to a vector of type nodHB.
 *	\param size The size of the vector tabHBAtoms.
 */
void free_tabHBAtoms(struct nodHB* tabHBAtoms,int size);

#endif //MEM_H
