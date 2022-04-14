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

/*! \file check.h
    \brief This file contains all headers and constants to display and verify content of variables.
*/
#ifndef CHECK_H
#define CHECK_H 

/*=======================Verification & tests===============*/
// #define DISP_FILE_LIST /*...Uncomment this label if you want to check the list of files in molecule->inputFNList...*/
//#define  DIS_BE /*...Uncomment this label if you want to check ...*/
//#define VERIF_VAL /*...Uncomment this label if you want to check ...*/
//#define VERIF_CONF /*...Uncomment this label if you want to check ...*/
#define VERIF_MEM_ACCESS /*...Uncomment this label if you want to check ...*/
#define VERIF_FILE_ACCESS /*...Uncomment this label if you want to check ...*/
// #define VERIF_HBONDS /*...Uncomment this label if you want to check ...*/
//#define DISP_MAT /*...Uncomment this label if you want to check ...*/
//#define DISP_CONF /*...Uncomment this label if you want to check ...*/
//#define DISP_COV /*...Uncomment this label if you want to check ...*/
// #define DISP_IMG /*...Uncomment this label if you want to check ...*/
// #define DISP_TAB_HASH /*...Uncomment this label if you want to check ...*/
    
/*! \fn  display_fileList(struct MyModel *molecule)
 *  \brief Display the content of the variable molecule->inputFNList : a list of files names.
 *  \param molecule A pointer to a structure of type MyModel.
 */
void display_fileList(struct MyModel *molecule);

/*! \fn display_cov_bonds(struct MyModel *molecule)
 *  \brief Display the covalent bonds based on  molecule->covBond.
 *  \param molecule A pointer to a structure of type MyModel.
 */
void display_cov_bonds(struct MyModel *molecule);

/*! \fn display_matrix(int **matrix, int sizeL,int sizeC, char matrix_name[256] )
 *  \brief Display the content of a two dimensional matrix of integers.
 *  \param matrix A pointer to a two dimensional matrix of integers. 
 *  \param sizeL An integer to indicate number of rows in the matrix.
 *  \param sizeC An integer to indicate number of columns in the matrix.
 *	\param matrix_name A string which indicates label of the matrix.
 */
void display_matrix(int **matrix, int sizeL,int sizeC, char matrix_name[256] );

/*! \fn display_img(struct atom* img, int size, char img_name[256])
 *  \brief Display the content of the vector img. It contains a list of atoms with their chemical types and their (x,y,z) positions.
 *  \param img  A pointer to a structure of type img. 
 *  \param size The size of the vector img.
 *	\param img_name A string which indicates the label of the vector.
 */
void display_img(struct atom* img, int size, char img_name[256]);

/*! \fn display_conf(char *conf, int size)
 *  \brief Display the content of the set of characters conf and count the number of '1' in it.  
 *  \param conf A pointer to a vector of characters. 
 *	\param size The size of the vector.
 */
void display_conf(char *conf, int size);

/*! \fn display_tabHash(struct nodHash* tabHash)
 *  \brief Display the content of the hash table which contains list of conformation hashed
 *  \param tabHash A pointer to a structure of type nodHash. A hash table. 
 */
void display_tabHash(struct nodHash* tabHash);

/*! \fn check_hbond_conf(struct confList* conf) 
 *  \brief Check and display if there is an atom which has been more than 2 times donor or more than 2 times acceptor. 
 *  \param conf A pointer to a structure of type confList. Pointer to a conformation.
 */
void check_hbond_conf(struct confList* conf);


#endif //CHECK_H

