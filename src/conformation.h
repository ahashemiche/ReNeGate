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

/*! \file conformation.h
    \brief This file contains all headers related to the treatments of the conformational dynamics analysis.
*/

#ifndef CONFORMATION_H
#define CONFORMATION_H

/*! \fn add_transit(struct MyModel* molecule,struct confList** pred, int numImg)
 *  \brief Analyse the current conformation (numImg) according to the conformations list (add it to the list, if it is not identified yet).
 *  \param molecule A pointer to a structure of type MyModel.
 *  \param pred A pointer to a structure of type confList. This contains the previous conformation.
 *	\param numImg The number of the current step (snapshot).
 *  \return A pointer to a structure of type confList. It is a pointer to the head of the conformations list.
 */
struct confList* add_transit(struct MyModel* molecule,struct confList** pred, int numImg);

/*! \fn update_conf_states(struct MyModel* molecule)
 *	\brief Update the list of conformations. Check if each conformation is stable or an intermediate state. 
 *  \param molecule A pointer to a structure of type MyModel. It uses molecule->conformations field.
 */
void update_conf_states(struct MyModel* molecule);

/*! \fn verif_state(struct confList* head, struct confList* courState, int nbAtom, int level, int *transf);
 *  \brief Verify if the conformation courState already exists in the conformations list (head).
 *  \param head A pointer to a structure of type confList. It is a pointer to the head of the conformations list.
 *  \param courState A pointer to a structure of type confList (current conformation). 
 *	\param nbAtom The number of atoms in the molecular system.
 *	\param level The level of comparison, a 4 bits number that indicates which bonds are taken into account in the comparison.  
 *	\param [out]transf A 7 bits number that contains changes in bonds in the current conformation courState.  
 *  \return A pointer to a structure of type confList of the conformation corresponding to the courState, NULL else. 
 */
struct confList* verif_state(struct confList* head, struct confList* courState, int nbAtom, int level, int *transf, int **bondDyn, int nbAtomDyn);

/*! \fn verif_isom(struct confList* confA,struct confList* confB,int level, int nbAtom,int *transf)
 *  \brief Verify if confA and confB are isomorphic or not. Comparison of adjacency matrices and isomorphism tests are performed.
 *  \param confA A pointer to a structure of type confList.
 *  \param confB A pointer to a structure of type confList.
 *	\param level The level of comparison, a 4 bits number that indicates which bonds are taken into account in the comparison.  
 *	\param nbAtom The number of atoms in the molecular system.
 *	\param [out]transf A 7 bits number that contains changes in bonds in the current conformation courState.  
 *  \return True if conformations confA and confB are isomorphic.
 */
bool verif_isom(struct confList* confA,struct confList* confB,int level, int nbAtom,int *transf, int **bondDyn, int nbAtomDyn);


/*! \fn add_bondDyn(int **bondDyn, int index1, int index2, int type)
 *  \brief Add the bond between atoms index1 and index2 as a dynamic bond
 *  \param bondDyn A pointer to a 2 dimensional matrix.
 *  \param index1 An index of an atom.
 *	\param index2 An index of an atom. 
 *	\param type the type of bond (1: CB, 2: HB, 3: IB, 4: MB).
 */
void add_bondDyn(int **bondDyn, int index1, int index2, int type);

/*! \fn construct_dreadnaut_file(struct confList* confA,struct confList* confB, int level, char inputFN[])
 *  \brief construct a dreadnaut file for the two conformation graphs, in order to apply the nauty program for the isomorphism.
 *  \param confA A pointer to a structure of type confList.
 *  \param confB A pointer to a structure of type confList.
 *	\param level The level of comparison, a 4 bits number that indicates which bonds are taken into account in the comparison.  
 *	\param [in,out]inputFN The name of the dreadnaut file. This file will contains the canonical graphs of confA and confB.
 */
void construct_dreadnaut_file(struct confList* confA,struct confList* confB, int level, char inputFN[]);

/*! \fn create_state(struct MyModel* molecule)
 *  \brief Create a state for the current conformation. Based upon molecule->covBond , molecule->bondHList, molecule->ionBond.
 *  \param molecule A pointer to a structure of type MyModel.
 *  \return A pointer to a structure of type confList. 
 */
struct confList* create_state(struct MyModel* molecule); 

/*! \fn add_state(struct MyModel* molecule,struct confList *newConf,int numImg)
 *  \brief Add the new conformation newConf at the "end" of the conformations list (molecule->conformations).
 *  \param molecule A pointer to a structure of type MyModel.
 *  \param newConf A pointer to a structure of type confList. The new conformation to add.
 *	\param numImg The number of the current step (snapshot).
 *  \return A pointer to a structure of type confList. It is a pointer to the head of the conformations list.
 */
struct confList* add_state(struct MyModel* molecule,struct confList *newConf,int numImg);

/*! \fn verif_succ(struct confList* head, int confName)
 *  \brief Verify if a node head has confName as successor.
 *  \param head A pointer to a structure of type confList.
 *  \param confName A name of a conformation.
 *  \return A pointer to structure of type successor where the conformation's name is confName.
 */
struct successor* verif_succ(struct confList* head, int confName);

/*! \fn new_succ(struct confList* head,struct confList* conf,int transf)
 *  \brief Create a new successor at the beginning of head->succ list of head.
 *  \param head A pointer to a structure of type confList.
 *  \param conf A pointer to a structure of type confList. The new successor of head.
 *	\param transf The difference between conf and head. 
 *  \return A pointer to structure of type successor. The beginning of head->succ list of head.
 */
struct successor* new_succ(struct confList* head,struct confList* conf,int transf);

/*! \fn get_bondcovdyn(struct MyModel* molecule,struct confList** conf)
 *  \brief Get rotational axes of the molecule in conformation conf (cov-bonds dynamics except leaves).
 *  \param molecule A pointer to a structure of type MyModel. This function changes the molecule->interBond matrix.
 *  \param conf A pointer to a structure of type confList.
 */
void get_bondcovdyn(struct MyModel* molecule,struct confList** conf);

/*! \fn add_isthm(struct MyModel *molecule, int index1, int index2)
 *  \brief Add an isthmus (index1-index2) to the isthmus list.
 *  \param molecule A pointer to a structure of type MyModel. This function uses molecule->isthList.
 *  \param index1 An index of an atom.
 *	\param index2 An index of an atom. 
 */
void add_isthm(struct MyModel *molecule, int index1, int index2);

/*! \fn verif_exist_isthm(int isthList[NB_MaxIsthm], int index1, int index2)
 *  \brief Verify if the isthmus (index1-index2) already exists in isthList list.
 *  \param isthList A vector of integers. The list of isthmus. 
 *  \param index1 An index of an atom.
 *	\param index2 An index of an atom. 
 *  \return The index of the isthmus in isthList list or -1 if not found.
 */
int verif_exist_isthm(int isthList[NB_MaxIsthm], int index1, int index2);

/*! \fn update_HbondList(struct MyModel *molecule, int index)
 *  \brief Delete the H-bond that has been never formed (always in position of HB).
 *  \param molecule A pointer to a structure of type MyModel.
 *  \param index The index of the H-Bond concerned. 
 */
void update_HbondList(struct MyModel *molecule, int index);

/*! \fn update_isthList(struct MyModel *molecule)
 *  \brief Update isthmus list to if there has been conformational rotation or a simple rotation (without change of conformations).
 *  \param molecule A pointer to a structure of type MyModel. This uses molecule->isthList and molecule->conformations.
 *  \return True if there has been a conformational rotation. An isthmus change from one conformation to another.
 */
bool update_isthList(struct MyModel *molecule);

/*! \fn update_ConfList(struct MyModel *molecule)
 *  \brief Update the list of conformations if there are similar one, once the H-bond were updated (some H-bonds may be deleted).
 *  \param molecule A pointer to a structure of type MyModel. This uses molecule->conformations.
 */
void update_ConfList(struct MyModel *molecule);

/*! \fn merge_period(struct confList* confA,struct confList* confB)
 *  \brief Merge all periods of confA in confB. It changes confA->imgList and deletes confB->imgList.
 *  \param confA A pointer to a structure of type confList.
 *  \param confB A pointer to a structure of type confList.
 */
void merge_period(struct confList* confA,struct confList* confB);

/*! \fn shift_succ( struct confList* head, struct successor* succ)
 *  \brief Delete succ from the list of successors of Head
 *  \param head A pointer to a structure of type confList.
 *  \param succ A pointer to a structure of type successor.
 */
void shift_succ( struct confList* head, struct successor* succ);

 /*! \fn get_suiv(struct MyModel *molecule,int numImg,int *index)
 *  \brief Find the conformation that appears at snapshot b+1.
 *  \param molecule A pointer to a structure of type MyModel.
 *  \param numImg The number of the current step (snapshot).
 *	\param [out]index The index in imgList of appearance/ disappearance of conformation.
 *  \return The pointer to the conformation which appears at step (snapshot) b+1.
 */
struct confList* get_suiv(struct MyModel *molecule,int numImg,int *index);

/*! \fn get_snap_index(struct confList *conf,int numImg)
 *  \brief Find in imgList of conformation conf the index of snapshot b+1
 *  \param conf A pointer to a structure of type confList.
 *  \param numImg The number of the current step (snapshot).
 *  \return The index of snapshot b+1 in imgList of conformation conf.
 */
int get_snap_index(struct confList *conf,int numImg);

/*! \fn shift_periods(struct confList *conf,int index)
 *  \brief Delete the period of index in imgList of conformation conf.
 *  \param conf A pointer to a structure of type confList.
 *  \param index An index in imgList of conformation conf.
 *  \return A pointer to structure of type confList. It returns the conformation conf with a new list of snapshots (periods).
 */
struct confList*  shift_periods(struct confList *conf,int index);

/*! \fn add_period(struct confList *conf,int index)
 *  \brief add a period in imgList of conformation conf, at index.
 *  \param conf A pointer to a structure of type confList.
 *  \param index An index in imgList of conformation conf.
 *  \return A pointer to structure of type confList. It returns the conformation conf with a new list of snapshots (periods).
 */
struct confList*  add_period(struct confList *conf,int index);

/*! \fn add_conf(struct MyModel* molecule,struct confList** pred, int numImg)
 *  \brief Add a conformation to the global conformation list. Based upon molecule->covBond , molecule->bondHList, molecule->ionBond.
 *  \param molecule A pointer to a structure of type MyModel.
 *  \param pred A pointer to a structure of type confList. The previous conformation.
 *	\param numImg The number of the current step (snapshot).
 *  \return A pointer to a structure of type confList. Update the conformations list according to the current conformation.
 */
struct confList* add_conf(struct MyModel* molecule,struct confList** pred, int numImg);

#endif //CONFORMATION_H