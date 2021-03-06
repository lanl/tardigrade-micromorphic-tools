# Micromorphic tools library
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   February 21, 2020
#
# This is the makefile for test_micromorphic_tools.cpp it tests the 
# micromorphic continuum tools.

include ../../../../config.mk

# List of executables
EXECUTABLE=test_micromorphic_tools
TESTOBJECT=test_micromorphic_tools.o
TESTSOURCE=test_micromorphic_tools.cpp

OBJECT=micromorphic_tools.o
SOURCE=../../micromorphic_tools.cpp
HEADER=../../micromorphic_tools.h

ERROROBJECT=error_tools.o

OBJECTS = $(TESTOBJECT)
OBJECTS += $(OBJECT) $(ERROROBJECT)

INC+=-I../..
LOCAL_CFLAGS=-ggdb

#Set the mode due to the method variable
DBG =
ifeq ($(METHOD),dbg)
  DBG = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g
endif
ifeq ($(METHOD),oprof)
  DBG = -pg
endif

# Makefile rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(TESTOBJECT) $(OBJECT) $(ERROROBJECT)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=c++11

$(TESTOBJECT): $(TESTSOURCE) $(HEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $< -c -o $@ $(CFLAGS) -lstdc++ -std=c++11

$(OBJECT): $(SOURCE) $(HEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(CONSTITUTIVEOBJECT): $(CONSTITUTIVESOURCE) $(CONSTITUTIVEHEADER) $(VECTORHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

clean:
	rm -f $(EXECUTABLE)
	rm $(OBJECTS)

.PHONY: all clean
