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

/*! \file io.h
    \brief This file contains all headers related to...
*/
#ifndef IO_H
#define IO_H 

/*! \fn save_event(struct MyModel *molecule, int index, int index2, char type ,char mode) 
 *  \brief Save a change that has been occurred in the trajectory. 
 *	Change includes : bonds appearance or disappearance, proton transfer, change in reference snapshot. 
 *  \param molecule A pointer to a structure of MyModel.
 *	\param index An index of an atom.
 *	\param index2 An index of an atom.
 *	\param type The type of change : '+' if bond created , '-' if broken.
 *	\param mode The mod used to open the output file. 
 */
void save_event(struct MyModel *molecule, int index, int index2, char type ,char mode);

/*! \fn save_snapshot(struct MyModel *molecule, struct confList* conf)
 *  \brief Save one snapshot for conformation conf into "file.xyz" file.
 *  \param molecule A pointer to a structure of MyModel.
 *	\param conf A pointer to a structure of type confList.
 */
void save_snapshot(struct MyModel *molecule, struct confList* conf);

/*! \fn save_cov_bonds(struct MyModel *molecule)
 *  \brief Save the covalent bonds molecule->covBond into "file.cov" file.
 *  \param molecule A pointer to a structure of MyModel. This function uses molecule->covBond field.
 */
void save_cov_bonds(struct MyModel *molecule);

/*! \fn save_Hbond_dynamics(struct MyModel *molecule)
 *  \brief Save H-bond dynamics into "file.Hbonds" and draw results using Gnuplot
 *  \param molecule A pointer to a structure of MyModel. This function uses molecule->bondHdyn field.
 */
void save_Hbond_dynamics(struct MyModel *molecule);

/*! \fn save_conformations(struct MyModel *molecule)
 *  \brief Save the identified conformations into "file.csv". This file will contain all details about 
 *	the global trajectory and for each conformation identified.
 *  \param molecule A pointer to a structure of MyModel. This function uses molecule->conformations.
 */
void save_conformations(struct MyModel *molecule);

/*! \fn save_dyn_desc(struct MyModel* molecule)
 *  \brief Save list of conformations found in the whole trajectories. This is used when many trajectories are analysed. 
 *  \param molecule A pointer to a structure of MyModel. This function uses molecule->conformations.
 */
void save_dyn_desc(struct MyModel* molecule);

/*! \fn save_conf_bonds(struct MyModel* molecule,struct confList* conf)
 *  \brief Save the bonds of a conformation conf into "conformation/bonds_conf/conf.txt" file. 
 *  \param molecule A pointer to a structure of MyModel.
 *	\param conf A pointer to a structure of type confList.
 *	\param type The type of the analysis performed ; single, multiple, isomers...
 */
void save_conf_bonds(struct MyModel* molecule,struct confList* conf, char type);

/*! \fn save_isom_bonds(struct ModelIso *isomList)
 *  \brief Save isomers list into "conf2iso.txt" file.
 *  \param isomList A pointer to a structure of ModelIso.
 */
void save_isom_bonds(struct ModelIso *isomList);

/*! \fn save_isom_xyz(char inputDir[500], int isom_num, int conf_num)
 *  \brief Save Cartesian coordinates for isomer isom_num that corresponds to conformation conf_num. 
 *  \param inputDir The directory name of trajectories. 
 *	\param isom_num The name of the concerned isomer.
 *	\param conf_num The name of the concerned conformation.
 */
void save_isom_xyz(char inputDir[500], int isom_num, int conf_num);

/*! \fn save_isom_traj(struct ModelIso *isomList)
 *  \brief Save the dynamics of isomers for each trajectory. This uses the intermediate files generated in the analyses of trajectories.
 *  \param isomList A pointer to a structure of ModelIso.
 */
void save_isom_traj(struct ModelIso *isomList);

/*! \fn draw_globalGraph(struct MyTrajGraph graphTraj)
 *  \brief Save the graphs of transitions (different paths explored) for the trajectories analysed.
 *  \param graphTraj An instance of a structure of type MyTrajGraph. This function uses graphTraj.graphList
 */
void draw_globalGraph(struct MyTrajGraph graphTraj);

/*! \fn draw_trajGraph(struct MyGraph *graph, char graphTrajFN[256], int num_traj)
 *  \brief Draw the graph of transition file (GraphViz file) related to trajectory num_traj
 *  \param graph A pointer to a structure of type MyGraph.
 *	\param graphTrajFN The file name of the output graphviz file. 
 *	\param num_traj A number of trajectory.
 */
void draw_trajGraph(struct MyGraph *graph, char graphTrajFN[256], int num_traj);

/*! \fn plot_Hbonds(struct MyModel *molecule)
 *  \brief Plot the time evolution of H-bonds identified along the analysed trajectory, using Gnuplot.
 *  \param molecule A pointer to a structure of MyModel.
 */
void plot_Hbonds(struct MyModel *molecule);


/*! \fn save_evol_conf(struct MyModel *molecule)
 *  \brief Save the order of appearance of conformations.
 *  \param molecule A pointer to a structure of MyModel.
 */
void save_evol_conf(struct MyModel *molecule);

/*! \fn plot_conf_periods(struct MyModel *molecule)
 *  \brief Plot the time evolution of the identified conformations, using Gnuplot.
 *  \param molecule A pointer to a structure of MyModel.
 */
void plot_conf_periods(struct MyModel *molecule);

/*! \fn draw_mixed_graph(struct MyModel *molecule, char type)
 *  \brief Draw the mixed graphs for each conformation identified along the trajectory (ies), using GraphViz.
 *  \param molecule A pointer to a structure of MyModel.
 *	\param type The type of the analysis performed ; single, multiple, isomers...
 */
void draw_mixed_graph(struct MyModel *molecule, char type);

/*! \fn draw_trans_graph(struct MyModel *molecule)
 *  \brief Draw the graph of transition using GraphViz
 *  \param molecule A pointer to a structure of MyModel.
 */
void draw_trans_graph(struct MyModel *molecule);

/*! \fn draw_distribution(char inputFN[], char outputFN[])
 *  \brief Draw the distribution of conformations according to number of H-bonds formed, using gnuplot.
 *  \param inputFN The input file name (data).
 *	\param outputFN The output file name.
 */
void draw_distribution(char inputFN[], char outputFN[]) ;

/*! \fn save_cc(char *outputFN, struct nodCC* tabCC,int size, int *maxCC)
 *  \brief Save the connected components found and their size.
 *  \param outputFN The output file name.
 *	\param tabCC A pointer to structure if type nodCC. It contains the list of connected components.
 *	\param size The size of vector tabCC. 
 *	\param [in,out]maxCC The size of the biggest connected component.
 */
void save_cc(char *outputFN, struct nodCC* tabCC,int size, int *maxCC);

/*! \fn save_global_values(char *outputFN,struct GraphModel* graphPoss)
 *  \brief Save global values for graph of possible conformations : 
 *	- Energy min,
 *	- Energy max,
 *	- Number of conformations, 
 *	- Number of connected components,etc.
 *  \param outputFN The output file name.
 *	\param graphPoss A pointer to structure of type GraphModel.
 */
void save_global_values(char *outputFN,struct GraphModel* graphPoss);

/*! \fn save_real_energy(struct MyModel* molecule)
 *	\brief Save the real energy of conformations found in simulation, from the xyz file.	
 *  \param molecule A pointer to a structure of MyModel.
 */
void save_real_energy(struct MyModel* molecule);

/*! \fn save_diff_with_ref(struct MyModel* molecule)
 *	\brief Save all the differences in terms of bonds/interactions between the reference an other conformers	
 *  \param molecule A pointer to a structure of MyModel.
 */
void save_diff_with_ref(struct MyModel* molecule);


/*! \fn get_label_conf(int confName, int size, char *fileName)
 *	\brief Return string which represents the label of the conformer confName in the graph of comparaison with the reference
 *  \param confName An integer which represents the number of a conformer.
 *  \param size An integer which represents the size of the trajectory (number of steps).
 *  \param fileName A string which represent the file name of the list of conformers with there energy's values 
 *  \return  A pointer to char which contains the string that represents the label of the conformer  
 */
char *get_label_conf(int confName, int size, char *fileName);

/*! \fn draw_pathConf(char *outputDir,struct MyModel *molecule,char *conf, int indexF,int indexdiff)
 *  \brief Draw the mixed graph of conformation beloging to a path between two conformations.
 *  \param outputDir The output directory name.
 *	\param molecule A pointer to a structure of MyModel.
 *	\param conf The name of conformation.
 *	\param indexF The index used to generate the output file name.
 *	\param indexdiff The index of the H-bonds that make conformation conf different.
 */
void draw_pathConf(char *outputDir,struct MyModel *molecule,char *conf, int indexF,int indexdiff);


/*! \fn check_atom_involved(struct MyModel *molecule, struct confList* conf, int index);
 *	\brief Check if atom with index index is involved in dynamic bond or not.	
 *	\param molecule A pointer to a structure of MyModel.
 *	\param conf A pointer to a structure of type confList.
 *	\param index An index of an atom.
 *  \return A positif integer (1,2,3 or 4) if the atom with index is involved  a dynamic bonds, 0 else.
 */
int check_atom_involved(struct MyModel *molecule, struct confList* conf, int index);


/*! \fn
 *	\brief Add the list of conformers where the bond index1-index2 , of type bondType is present (output file _conf.txt)
 *	\param outputF A pointer to file.
 *	\param index1 An index of an atom.
 *	\param index2 An index of an atom.
 *	\param bondType Type of the bond index1-index2
 *	\param conf A pointer to a structure of type confList.
 */
void add_conf_list(FILE *outputF,int index1,int index2, int bondType, struct confList  *conf);


/*! \fn
 *	\brief  
 *	\param conf A pointer to a structure of type confList
 *	\param index An index of an atom.
 */
int get_coord_num(struct confList *conf, int index);

#endif //IO_H