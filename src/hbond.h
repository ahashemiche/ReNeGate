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

/*! \file hbond.h
    \brief This file contains all headers related to treatment of hydrogen bonds.
*/

#ifndef HBOND_H
#define HBOND_H 

/*! \fn get_orbits(struct MyModel *molecule,double alpha)
 *  \brief Compute the orbits around the hydrogen atoms. This updates the molecule->bondHList field.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param alpha A constant used to define the orbits radius. 
 */
void get_orbits(struct MyModel *molecule,double alpha);

/*! \fn update_donHacc_orbit(struct MyModel *molecule,int indexH, int *oldOrbit, int size)
 *  \brief Update donHacc if orbit has changed
 *  \param molecule A pointer to a structure of MyModel.
 *	\param indexH  The index of the hydrogen atom that has a change in its orbit.
 *	\param oldOrbit A pointer to a vector of integers. The list of atoms in an old orbit of the hydrogen atom.
 *  \param size The size of the donHacc matrix.
 */
void update_donHacc_orbit(struct MyModel *molecule,int indexH, int *oldOrbit, int size);

/*! \fn exit_atom_inOrbit(int *tabOrbit, int  size, int index)
 *  \brief Check if the atom with position index belongs to the orbit.
 *  \param tabOrbit A pointer to a vector of integers. A list of atoms in an orbit.
 *  \param size The size of the vector (orbit).
 *  \return True if index belongs to tabOrbit.
 */
bool exit_atom_inOrbit(int *tabOrbit, int  size, int index);

/*! \fn init_Hbonds(struct MyModel *molecule)
 *  \brief Initialize the H-bonds list. This initialize the molecule->bondHList field.
 *  \param molecule A pointer to a structure of MyModel.
 */
void init_Hbonds(struct MyModel *molecule);

/*! \fn get_don(int H,struct MyModel *molecule)
 *  \brief Get the index of donor related to the hydrogen atom with index H
 *  \param H An index of a hydrogen atom.
 *  \param molecule A pointer to a structure of MyModel.
 *  \return The index of the heavy atom (donor) bonded to atom H.
 */
int get_don(int H,struct MyModel *molecule);

/*! \fn verify_Hbond(struct MyModel *molecule,int D, int H, int A,int *state, char *sens)
 *  \brief Check if an H bond can be created between atoms of indexes D (heavy atom), H (hydrogen atom), and A (acceptor).
 *	In this function, the proton transfer is also checked.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param D The index of the heavy atom (donor).
 *	\param H The index of the hydrogen atom.
 *	\param A The index of the acceptor atom. 
 *	\param state[in,out] The state of the H-bond: 1 if it is formed, 0 if atoms are in position of H-bond, -1 for no H-bond.
 *	\param sens[in,out] The sense of the H-bond: 1 if a proton transfer, 0 else.
 *  \return True if a hydrogen bond can be created between the three atoms D, H, and A.
 */
bool verify_Hbond(struct MyModel *molecule,int D, int H, int A,int *state, char *sens);

/*! \fn check_new_Hbond(struct MyModel *molecule, int index)
 *  \brief Check if a new H-bond is created with atom of position index
 *  \param molecule A pointer to a structure of MyModel.
 *  \param index The index of an hydrogen atom.
 *  \return The index of the atom that verify the H-bond conditions, it returns -1 if no atom respects these conditions.
 */
int check_new_Hbond(struct MyModel *molecule, int index);

/*! \fn change_cov(struct MyModel *molecule, int index)
 *  \brief Change on the covalent bonds related to the hydrogen atom index. This is called if proton transfer has been occurred. 
 *  \param molecule A pointer to a structure of MyModel.
 *  \param index The index of hydrogen atom concerned by the proton transfer.
 */
void change_cov(struct MyModel *molecule, int index);

/*! \fn update_AH(struct MyModel *molecule, int index)
 *  \brief Update the H-bonds states (Donor,Acceptor) of hydrogen atom of position index 
 *  \param molecule A pointer to a structure of MyModel.
 *  \param index The index of hydrogen atom.
 */
void update_AH(struct MyModel *molecule, int index);

/*! \fn get_indexP(struct MyModel *molecule, int index, int index2)
 *  \brief Get index of atom with index2 in the list of orbit of hydrogen atom index.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param index The index of an hydrogen atom.
 *	\param index2 An index of an atom.
 *  \return The index of index2 in the orbit of the hydrogen atom index.
 */
int get_indexP(struct MyModel *molecule, int index, int index2);

/*! \fn imp_hbond(struct MyModel *molecule, int i, int j, char sens)
 *  \brief Mark the atoms involved in a H-bond.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param i An index of an atom.
 *  \param j An index of an atom.
 *  \param sense The orientation of the H-bond formed by the atoms i (heavy atom) and j (acceptor).
 */
void imp_hbond(struct MyModel *molecule, int i, int j, char sens);

/*! \fn get_Hbonds_dynamics(struct MyModel *molecule)
 *  \brief Dynamic analysis of H-bonds (proton transfer included).
 *	This function identify the H-bonds formed at molecule->num_img. It uses the molecule->bondHList field.
 *  \param molecule A pointer to a structure of MyModel.
 *  \return True if the H-bonds list has been changed according to the previous list ( previous snapshot).
 */
bool get_Hbonds_dynamics(struct MyModel *molecule);

/*! \fn add_Hbond(struct MyModel* molecule,int bondH[3], int state)
 *  \brief Add dynamics of H-bond bondH to molecule->bondHdyn.
 *  \param molecule A pointer to a structure of MyModel.
 *	\param bondH A vector containing the indexes of atoms involved in H-bond (Donor, Hydrogen,Acceptor).
 *  \param state The state of the H-bond.
 */
void add_Hbond(struct MyModel* molecule,int bondH[3], int state);

/*! \fn add_new_Hbond(struct MyModel* molecule,int bondH[3], int state)
 *  \brief Add a new H-bond to the "bondHdyn" list.
 *  \param molecule A pointer to a structure of MyModel.
 *	\param bondH A vector containing the indexes of atoms involved in H-bond (Donor, Hydrogen,Acceptor).
 *  \param state The state of the H-bond.
 */
void add_new_Hbond(struct MyModel* molecule,int bondH[3], int state);

/*! \fn exist_Hbond(struct MyModel* molecule, int bondH[3])
 *  \brief Check if the H-bond bondH already exists in bondHdyn list.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param bondH A vector containing the indexes of atoms involved in H-bond (Donor, Hydrogen,Acceptor).
 *  \return True if the H-bond bondH already exists.
 */
int exist_Hbond(struct MyModel* molecule, int bondH[3]);

/*! \fn update_donHacc(int **donHacc,int D,int H, int A, int state,char sens, int size)
 *  \brief Update the matrix of donor-hydrogen-acceptor.
 *  \param donHacc A pointer to a 2 dimensional matrix of integers. 
 *  \param D The index of the heavy atom (donor).
 *	\param H The index of the hydrogen atom.
 *	\param A The index of the acceptor atom. 
 *  \param state The state of the H-bond.
 *  \param sense The orientation of the H-bond formed by the atoms D (heavy atom), H (hydrogen), and A (acceptor).
 *  \param size The size of the donHacc matrix.
 */
void update_donHacc(int **donHacc,int D,int H, int A, int state,char sens, int size);

/*! \fn verif_change_Hbonds(struct MyModel *molecule,int **donHacc,int **copy_donHacc,int size)
 *  \brief Verify if there has been a change in the hydrogen bonds.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param donHacc A pointer to a 2 dimensional matrix of integers. 
 *	\param copy_donHacc A pointer to a 2 dimensional matrix of integers. A copy of donHacc matrix.
 *	\param size The size of the donHacc matrix. 
 *  \return Check if the H-bonds list has been changed, by comparing the donHacc and copy_donHacc matrices.
 */
bool verif_change_Hbonds(struct MyModel *molecule,int **donHacc,int **copy_donHacc,int size);

/*! \fn update_donHacc_ptr(int **donHacc,int *numOrb, int indexH,int indexOrb,int size)
 *  \brief If there has been a proton transfer, set all atoms in the orbit (different of indexOrb) of indexH to -1.
 *  \param donHacc A pointer to a 2 dimensional matrix of integers.
 *	\param numOrb A pointer to a vector of integers. The orbit of the atom indexH.
 *	\param indexH The index of the hydrogen atom that has a change in its orbit.
 *	\param indexOrb The atom index, concerned by the hydrogen bond.
 *  \param size The size of the donHacc matrix
 */
void update_donHacc_ptr(int **donHacc,int *numOrb, int indexH,int indexOrb,int size);

/*! \fn verif_pot_Hbond(struct MyModel *molecule, int D,int H, int A )
 *  \brief Verify if three atoms D,H and A can potentially form a hydrogen bond.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param D The index of the heavy atom (donor).
 *	\param H The index of the hydrogen atom.
 *	\param A The index of the acceptor atom. 
 *  \return True if the three atoms D,H and A can potentially form a hydrogen bond. This by comparing their distance.

 */
bool verif_pot_Hbond(struct MyModel *molecule, int D,int H, int A );

/*! \fn clean_Hbond_dynamics(struct MyModel* molecule)
 *  \brief Clean Hbonds according percentage of appearance and disappearance.
 *	\bug : this function is not used yet in our program.
 *  \param molecule A pointer to a structure of MyModel.
 */
void clean_Hbond_dynamics(struct MyModel* molecule);


/*! \fn check_nb_acc_per_hbond(struct MyModel *molecule)
 *  \brief Check that each hydrogen atom has at most 2 hydrogen bond 
 *  \param molecule A pointer to a structure of MyModel.
 */
void check_nb_acc_per_hbond(struct MyModel *molecule);

#endif //HBOND_H