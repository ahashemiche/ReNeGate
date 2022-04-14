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

/*! \file usage.h
    \brief This file contains all headers related to the usage of the program.
*/

#ifndef USAGE_H

/**
 * \fn  usage(int argc, char *argv[])
 * \brief Read arguments introduced by user and set the corresponding options.
 * \param argc  number of parameters introduced by the user.
 * \param argv  table containing all parameters introduced by the user.
 * \return  \e the list of options to be performed by the Getaway analyses.
 */
struct MyOpt usage(int argc, char *argv[]);

/*  \fn struct MyOpt default_val() 
 *	\brief The function set the default values to analyse the trajectories?
 *  \return A pointer to a structure of type MyOpt
 */
struct MyOpt default_val();

/*!
 * \fn display_usage(char *argv[])
 * \brief  Display the usages of the different parts of the Getaway program 
 * \details Getaway program is divided into 3 modules:
 * + singanalysis : this module consists in analysing a single trajectory
 * + multanalysis : this module consists in analysing many trajectories simultaneously 
 * + interface : this module consists in analysing trajectory in condensed matter with more treatments ()
 * \param argv  table containing all parameters introduced by the user.
 */
void display_usage(char *argv[]);

/*!
 * \fn get_anal_opt(struct MyModel* molecule,struct MyOpt opt )
 * \brief Set the different options and threshold values in the global variable \a molecule according to \a opt
 * \param molecule global variable used to contain all parameters and results of the analyses
 * \param opt  list of options set by the user
 */
void get_anal_opt(struct MyModel* molecule,struct MyOpt opt );

/*!
 * \fn get_parameter(struct MyModel *molecule)
 * \brief Get the threshold values used to calculate the different types bonds (covalent, hydrogen bonds and electrostatic interactions)
 * \param molecule global variable used to contain all parameters and results of the analyses
 */
void get_parameter(struct MyModel *molecule);
   
#define USAGE_H 
#endif //USAGE_H