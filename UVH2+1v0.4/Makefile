PWD = `pwd`

DEBUG =
OPTIMIZATION = -O3
FLOWTRACE =
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE)
COMPILER = g++
LIBS = -Lgsl/lib -lm -lgsl -lgslcblas
INCLUDES = -I. -Igsl/include/ 

.SUFFIXES: .o .cpp .h

.cpp.o :
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c $*.cpp
S1 =\
paramreader.cpp\
UVH2+1.cpp
S2 =\
paramreader.cpp\
convert.cpp\
convertnew.cpp
S3 =\
extract.cpp
S3b =\
paramreader.cpp\
extractfull.cpp
S3c =\
paramreader.cpp\
extractfullangle.cpp
S4 =\
paramreader.cpp\
initE.cpp
S5 =\
paramreader.cpp\
givedist.cpp
S6 =\
paramreader.cpp\
prereso.cpp
S7 =\
paramreader.cpp\
preresofull.cpp
S7b =\
paramreader.cpp\
preresofull.cpp
OBJ1 =\
paramreader.o\
UVH2+1.o
OBJ2 =\
paramreader.o\
convert.o
OBJ2b=\
paramreader.o\
convertnew.o
OBJ2c=\
paramreader.o\
convertfull.o
OBJ2d=\
paramreader.o\
convertlinear.o
OBJ3 =\
extract.o
OBJ3b =\
paramreader.o\
extractfull.o
OBJ3c =\
paramreader.o\
extractfullangle.o
OBJ4=\
paramreader.o\
initE.o
OBJ5=\
paramreader.o\
givedist.o
OBJ6 =\
paramreader.o\
prereso.o
OBJ7 =\
paramreader.o\
preresofull.o
OBJ7b =\
paramreader.o\
preresofullangle.o
EXE1 =\
vh2
EXE2 =\
convert
EXE2b =\
convertnew
EXE2c =\
convertfull
EXE2d =\
convertlinear
EXE3 =\
extract
EXE3b =\
extractfull
EXE3c =\
extractfullangle
EXE4 =\
initE
EXE5 =\
givedist
EXE6 =\
prereso
EXE7 =\
preresofull
EXE7b =\
preresofullangle

$(EXE1) : $(OBJ1) 
$(EXE2) : $(OBJ2) 
$(EXE2b) : $(OBJ2b) 
$(EXE2c) : $(OBJ2c) 
$(EXE2d) : $(OBJ2d)
$(EXE3) : $(OBJ3) 
$(EXE3b) : $(OBJ3b)
$(EXE3c) : $(OBJ3c)
$(EXE4) : $(OBJ4)
$(EXE5) : $(OBJ5)
$(EXE6) : $(OBJ6)
$(EXE7) : $(OBJ7)
$(EXE7b) : $(OBJ7b)

vh2:	
	$(COMPILER) $(OBJ1) -o $(EXE1)  $(LIBS)

convert:
	$(COMPILER) $(OBJ2) -o $(EXE2)  $(LIBS)

convertnew:
	$(COMPILER) $(OBJ2b) -o $(EXE2b)  $(LIBS)

convertfull:
	$(COMPILER) $(OBJ2c) -o $(EXE2c)  $(LIBS)

convertlinear:
	$(COMPILER) $(OBJ2d) -o $(EXE2d)  $(LIBS)

extract:
	$(COMPILER) $(OBJ3) -o $(EXE3)  $(LIBS)

extractfull:
	$(COMPILER) $(OBJ3b) -o $(EXE3b)  $(LIBS)

extractfullangle:
	$(COMPILER) $(OBJ3c) -o $(EXE3c)  $(LIBS)

initE:
	$(COMPILER) $(OBJ4) -o $(EXE4)  $(LIBS)

givedist:
	$(COMPILER) $(OBJ5) -o $(EXE5) $(LIBS)
	
prereso:
	$(COMPILER) $(OBJ6) -o $(EXE6) $(LIBS)

preresofull:
	$(COMPILER) $(OBJ7) -o $(EXE7) $(LIBS)

preresofullangle:
	$(COMPILER) $(OBJ7b) -o $(EXE7b) $(LIBS)

# clean up misc files
clean :
	rm -f $(EXE1) $(EXE2) $(EXE2b) $(EXE2c) $(EXE2d) $(EXE3) $(EXE4) $(EXE5) $(EXE6) *\.o *~ #* core 
