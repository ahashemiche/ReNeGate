#*****************************************************************************
#    GATEwAY - GrAph ThEory for conformAtional dYnamics 
#
#    Copyright (c) 2014-2022 Sana Bougueroua
#                  2020-2022 Ali Hashemi
#    Please cite:  J. Chem. Phys. 2018, 149 (18), 184102.         (DOI 10.1063/1.5045818 )
#		   	
#    This file written by Sana Bougueroua.
#
#    ---------------------------------------------------------------------------
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************


#References and Variables

#Compiler
CC ?=gcc
#Flags for complication 
CFLAGS ?= -Wall 
#Flags for linker 
LDFLAGS ?= -lm
#Directory of code files (.c,.h)
SRC_DIRS ?= ./src
#Directory of object files (.o)
BUILD_DIR ?= ./obj
#Directory of executable files 
BIN_DIR ?= ./bin
TARGET_EXEC ?= singanalysis

SRCS := $(shell find $(SRC_DIRS) -name *.c )
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)


#======
#Linker
#======
$(BIN_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS) 


#===========
#Compilation
#===========
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

#========
#Example
#========
.PHONY: example

example:
	./bin/singanalysis -s ./example/sample_MnNN_Alkoxide_KBr.xyz -c -d 10 

#========
#Cleaning
#========
.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)
	$(RM) -r $(BIN_DIR)/$(TARGET_EXEC)


-include $(DEPS)

MKDIR_P ?= mkdir -p


#=================================
#Installation of required softwares
#=================================
# 3 external free software are used in our code : dreadnaut (nauty for isomorphism), gnuplot (plots), dot and neato (graphviz for graphs)
# For MacOS 
# sudo port install nauty 
# sudo port install gnuplot
# sudo port install graphviz

# For Linux 
# sudo apt install nauty 
# sudo apt install gnuplot
# sudo apt install graphviz
