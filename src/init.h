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

/*! \file init.h
    \brief This file contains all headers related to initialization 
*/

#ifndef INIT_H
#define INIT_H 

/*! \fn init_model(struct MyModel* molecule)
 *  \brief initialize the model, create result directories and initialize the adjacency matrices.
 *  \param molecule A pointer to a structure of type MyModel.  
 */
void init_model(struct MyModel* molecule);

/*! \fn get_nbAtom(char inputFN[])
 *  \brief Compute the size of one snapshot in the trajectory (#atoms).
 *  \param inputFN[] The file name of the trajectory to be analysed.
 *  \return The size of one snapshot (number of atoms)
 */
int get_nbAtom(char inputFN[]); 

/*! \fn get_nbMolecule(char inputFN[])
 *  \brief Compute the size of trajectory (#snapshots).
 *  \param inputFN[] The file name of the trajectory to be analysed.
 *  \return The number of snapshots in the trajectory.
 */
int get_nbMolecule(char inputFN[]); 

/*! \fn get_ion_atoms(struct MyModel *molecule)
 *  \brief Compute the number of ions (cations, anions) in the molecular system (one snapshot).
 *	In this version we take into account : Li, Cl, Ar. One can add more ions.
 *  \param molecule A pointer to a structure of type MyModel.
 *  \return The number of ions (cations, anions) in the molecular system.
 */
int get_ion_atoms(struct MyModel *molecule);

/*! \fn get_metal_atoms(struct MyModel *molecule)
 *  \brief Compute the number of metals in the molecular system (one snapshot).
 *	In this version we take into account : Mn. One can add more metals.
 *  \param molecule A pointer to a structure of type MyModel.
 *  \return The number of metals in the molecular system.
 */
int get_metal_atoms(struct MyModel *molecule);

/*! \fn get_molecule(FILE *traj,struct MyModel *molecule)
 *  \brief Read the current snapshot. This function uses the molecule->img field.
 *  \param traj The input file of a trajectory to be read.
 *  \param molecule A pointer to a structure of type MyModel.
 */
void get_molecule(FILE *traj,struct MyModel *molecule);

/*! \fn treat_snapshot(FILE *traj,struct MyModel* molecule)
 *  \brief Treat the current snapshot by computing bonds (depending on the level set by the user). 
 *	Check if there has been a change in bonds according to the previous snapshot.
 *  \param traj The input file of a trajectory to be read.
 *  \param molecule A pointer to a structure of type MyModel.
 *  \return True if there has been a change in bonds between the current snapshot and the previous one.
 */
bool treat_snapshot(FILE *traj,struct MyModel* molecule);

/*! \fn treat_trajectory(struct MyModel* molecule)
 *  \brief Analyse the dynamics of the trajectory in terms of conformations. 
 * 	The function indeed, analyse all snapshots in the trajectory, identify the conformation (construct the mixed graphs) and its transitions.
 *  \param molecule A pointer to a structure of type MyModel.
 */
void treat_trajectory(struct MyModel* molecule);

/*! \fn get_molecule_name(struct MyModel* molecule)
 *  \brief Get the Empirical formula according the atoms of the molecular system. This functions uses molecule->name field.
 *  \param molecule A pointer to a structure of type MyModel.
 */
void get_molecule_name(struct MyModel* molecule);

/*! \fn get_new_ref(struct atom* img, int size)
 *  \brief Change the reference snapshot. 
 *  \param img A pointer to a structure of type atom.
 *  \param size The size of the snapshot (number of atoms)
 *  \return A pointer to a structure of type atom. This pointer will be the new reference snapshot.
 */
struct atom* get_new_ref(struct atom* img, int size);

/*! \fn add_atom(char *atomName, int nbAtom)
 *  \brief Create chain with name of atom and its number of occurrences.
 *  \param atomName A name of an atom.
 *  \param nbAtom The number of occurrence of atomName. 
 *  \return A chain with the name of atom concatenated to its number of occurrences.
 */
char* add_atom(char *atomName, int nbAtom);


#endif //INIT_H