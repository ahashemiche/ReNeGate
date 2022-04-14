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

/*! \file intermolecular.h
    \brief This file contains all headers related to treatment of intermolecular electrostatic interactions.
*/

#ifndef INTERMOLECULAR_H
#define INTERMOLECULAR_H 

/*! \fn get_ion(struct MyModel* molecule)
 *  \brief Identify the intermolecular interactions and check if there has been a change with the previous snapshot.
 *  \param molecule A pointer to a structure of MyModel.
 *  \return True if there has been a change in the intermolecular interactions. 
 */
bool get_ion(struct MyModel* molecule);

/*! \fn verif_ion_bond(struct MyModel* molecule, int index1, int index2)
 *  \brief Verify if atoms of indexes index1 and index2 form an intermolecular interaction.
 *  \param molecule A pointer to a structure of MyModel.
 *  \param index1 An index of an atom.
 *  \param index2 An index of an atom.
 *	\return True if atoms of indexes index1 and index2 form an intermolecular interaction.
 */
bool verif_ion_bond(struct MyModel* molecule, int index1, int index2); //check the parameter 	


#endif //INTERMOLECULAR_H