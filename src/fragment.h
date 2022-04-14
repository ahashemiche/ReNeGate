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

/*! \file fragment.h
    \brief This file contains all headers related to treatment of fragments.
*/

#ifndef FRAGMENT_H
#define FRAGMENT_H 

/*! \fn void get_frag_dynamics(struct MyModel* molecule)
 *  \brief Get the list of fragments that have been observed along the trajectory
 *  \param molecule A pointer to a structure of MyModel.
 */
void get_frag_dynamics(struct MyModel* molecule);

/*! \fn  
 *  \brief  
 *  \param  
 *  \return 
 */
struct frgList* create_fragment(struct confList* conf,int **listfrag,int index, int nbConf);

/*! \fn  frgList* exist_frag(struct frgList* currfrg, struct frgList* fragments)
 *  \brief  
 *  \param  
 *  \return 
 */
struct frgList* exist_frag(struct frgList* currfrg, struct frgList* fragments);

/*! \fn  
 *  \brief  
 *  \param  
 *  \return 
 */
void construct_dreadnaut_file_frg(struct frgList* frgA,struct frgList* frgB, char inputFN[]);

/*! \fn  
 *  \brief  
 *  \param  
 *  \return 
 */
void save_frg(struct MyModel* molecule, struct frgList* frg);

/*! \fn  
 *  \brief  
 *  \param  
 *  \return 
 */
void draw_atom(FILE *output, struct atom *img, int **CB, int nbAtom, int index);

/*! \fn  
 *  \brief  
 *  \param  
 *  \return 
 */
void save_frg_list(struct frgList* fragments, int nbFrg, int nbConf, char *outputFN);

/*! \fn  
 *  \brief  
 *  \param  
 *  \return 
 */
void get_motifs(struct MyModel *molecule);

/*! \fn  
 *  \brief  
 *  \param  
 *  \return 
 */

/*! \fn  
 *  \brief  
 *  \param  
 *  \return 
 */

#endif //FRAGMENT_H