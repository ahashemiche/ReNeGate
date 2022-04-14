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

/*! \file covalent.h
    \brief This file contains all headers related to treatment of covalent bonds
*/

#ifndef COVALENT_H
#define COVALENT_H

/*!	\fn get_cov(struct MyModel* molecule)
 * 	\brief Compute covalent bonds for the molecule at instant molecule->num_img. The function update the "molecule->covBond" field. 
 * 	\param molecule A pointer to a structure of MyModel.  
 * 	\return True if there was a change in covalent bonds.
 */
bool get_cov(struct MyModel* molecule);

/*! \fn check_cov_bond(struct MyModel* molecule)
 * 	\brief Check the covalent bonds found in the get_cov function at instant "molecule->num_img".
 *	This consists in adding what is missed and deleting what is over.
 * 	\param molecule A pointer to a structure of MyModel. The procedure update the "molecule->covBond" field.  
 */
void check_cov_bond(struct MyModel* molecule);

/*! void get_inter_bonds(struct MyModel* molecule)
 * \brief Provide the internal covalent bonds NI-NI of molecule.
 *	It uses the covalent bonds identified with get_cov function (adjacency matrix "molecule->covBond").
 *	A covalent bond is considered as internal, if it is related with at least two other covalent bonds. 
 * 	The procedure updates the "molecule->interBond" field.
 * \param molecule A pointer to a structure of MyModel.   
 */
void get_inter_bonds(struct MyModel* molecule);

/*! \fn check_inter_node(struct MyModel* molecule, int a, int b)
 *  \brief Check if a-b is an internal covalent bond. It use the covalent bonds adjacency matrix.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param a An atom index 
 *	\param b An atom index
 *  \return True if a-b is an internal covalent bond. 
 */
bool check_inter_node(struct MyModel* molecule, int a, int b);

/*! \fn verif_change_cov(struct MyModel *molecule,int **covBond,int **covCopy,int size)
 *  \brief Verify if there have been a change in the covalent bonds. 
 *	We suppose that covbond and covCopy have the same order and number of atoms. 
 * 	The comparison is done by comparing the adjacency matrices.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param covBond A pointer to a 2 dimensional matrix of integers.
 *	\param covCopy A pointer to a 2 dimensional matrix of integers. A copy of covBond matrix.
 *	\param size The size of the matrix
 *  \return True if there is a difference between the two matrices covbond and covCopy.
 */
bool verif_change_cov(struct MyModel *molecule,int **covBond,int **covCopy,int size);

/*! \fn verif_cov_bond(struct MyModel* molecule, int i, int j)
 *  \brief Check if there is a bond between atoms i and j. It verify the type of the atoms i and j and compute 
 *	the distance using the covalent radii.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param i An atom index.
 *	\param j An atom index.
 *  \return True if i and j can form a covalent bond. 
 */
bool verif_cov_bond(struct MyModel* molecule, int i, int j) ;

/*! \fn get_covr_value(struct MyModel* molecule, char atomType)
 *  \brief Get the covalent radius of the atom with type atomType. 
 *  \param molecule A pointer to a structure of MyModel.
 *  \param atomType The type of atom ('C' for Carbon, 'N' for Nitrogen, etc.)
 *  \return The covalent radius value (in Angstroms)
 */
double get_covr_value(struct MyModel* molecule, char* atomName);

/*! \fn sens_Hbond(struct MyModel* molecule, int index)
 *  \brief Verify the orientation of H-bonds (proton transfer or not)
 *  \param molecule A pointer to a structure of MyModel.
 *  \param index The index of the hydrogen atom to be analysed. 
 *  \return 1 if there is a proton transfer, 0 else.
 */
int sens_Hbond(struct MyModel* molecule, int index);

/*! \fn add_bond(struct MyModel *molecule, int i, int j)
 *  \brief Add a covalent bond between two atoms with indexes i and j. It uses the molecule->covBond adjacency matrix,
 *	and updates the number of bonds formed for each atom.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param i An atom index.
 *	\param j An atom index.
 */
void add_bond(struct MyModel *molecule, int i, int j);

/*! \fn del_bond(struct MyModel *molecule, int i, int j)
 *  \brief Delete a covalent bond between two atoms with indexes i and j
 *  \param molecule A pointer to a structure of MyModel. It uses the molecule->covBond adjacency matrix,
 *	and updates the number of bonds formed for each atom.
 *  \param i An atom index.
 *	\param j An atom index.
 */
void del_bond(struct MyModel *molecule, int i, int j);

/*! \fn check_nearb(struct MyModel *molecule, int i,int type)
 *  \brief Check if i has a covalent bond with an atom with respect to some conditions.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param i An atom index.
 *	\param type  The type of condition: '+' for an atom with more covalent bonds that he has to form, 
 *	and '-' for atom with lack of covalent bonds.
 *  \return The index of the nearest atom to atom with index i.
 */
int check_nearb(struct MyModel *molecule, int i,int type);


/*! \fn del_bond(struct MyModel *molecule, int i, int j)
 *  \brief Check if there is any triangle structure formed by the covalent bonds and delete it
 *  \param molecule A pointer to a structure of MyModel. 
 */
void check_triangle_bond(struct MyModel* molecule);    
#endif //COVALENT_H