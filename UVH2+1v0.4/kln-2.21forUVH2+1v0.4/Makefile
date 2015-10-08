# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

#include Makefile.arch
VERSION       = 2.21
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"
CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC -g
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -shared

#------------------------------------------------------------------------------

HWORLD        = kln.exe

INC  = Bases.h Integral.h Regge96.h OverLap.h\
       CentralityCut.h Random.h PyRand.h\
       Saturation.h KLNModel.h\
       FTDipole.h\
       UnintegPartonDist.h KLNfunc.h KLNGamfunc.h\
       MCSaturation.h\
       Particle.h\

SRC  = ../paramreader.cpp\
       Bases.cxx Integral.cxx Regge96.cxx OverLap.cxx\
       CentralityCut.cxx Random.cxx PyRand.cxx\
       Saturation.cxx KLNModel.cxx\
       FTDipole.cxx\
       MCSaturation.cxx\
       main.cxx

OBJ          = $(SRC:.cxx=.o)

PROGRAMS      = kln.exe

OBJS         += $(GUITESTO) $(GUIVIEWERO) $(TETRISO)
PROGRAMS     += $(GUITEST) $(GUIVIEWER) $(TETRISSO)

INCLUDES = -I. -Igsl/include/
LIBS = -L. -lm -Lgsl/lib/ -lgsl -lgslcblas
SAVE = Makefile $(INC) $(SRC)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    Aclock Hello Tetris

all:            $(PROGRAMS)

$(HWORLD):      $(OBJ)
		$(LD) $(LDFLAGS) $(INCLUDES) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"


clean:
		@rm -f $(OBJ) core

distclean:      clean
		@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations
		@rm -rf cxx_repository
		-@cd RootShower && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)

###

Event.$(ObjSuf): Event.h
MainEvent.$(ObjSuf): Event.h
Bases.o: Bases.h
Integral.o: Integral.h
Regge96.o: Regge96.h
OverLap.o: OverLap.h
CentralityCut.o: CentralityCut.h
Saturation.o: Saturation.h
KLNModel.o: KLNModel.h
FTDipole.o: FTDipole.h


.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

dist: $(SAVE)
	tar cfvz kln-$(VERSION).tar.gz $(SAVE)

