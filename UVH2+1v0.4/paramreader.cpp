#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <stdlib.h>


// external vars defined in solve.cpp which are loaded here
extern int    NUMT,UPDATE,SNAPUPDATE;
extern long int STEPS;
extern double AT,EPS,TINIT,ETAOS,TF,TSTART,COEFF;
extern double B,L1COEF,L2COEF;
extern double IC;
extern int PTASIZE,PHIPASIZE;
extern int FREEZE;
extern double PTMAX, TRIEPS, TRIANGLE, QUADEPS, QUADANGLE;
extern double QUINTEPS, QUINTANGLE, SEXEPS, SEXANGLE, SEPTEPS, SEPTANGLE;
extern double BIEPS, BIANGLE;

using namespace std;

// this workhorse examines a key to see if corresponds to a var we are setting
// and then attempts to set the var corresponding to key by converting value to the
// appropriate type.  lots of hardcoding here
void setParameter(char *key, char *value) {
	// integer params
	if (strcmp(key,"NUMT")==0) NUMT=atoi(value);
	if (strcmp(key,"B")==0) B=atof(value);
	if (strcmp(key,"L1COEF")==0) L1COEF=atof(value);
	if (strcmp(key,"L2COEF")==0) L2COEF=atof(value);
	if (strcmp(key,"TINIT")==0) TINIT=atof(value);
	if (strcmp(key,"STEPS")==0) STEPS=atoi(value);
	if (strcmp(key,"UPDATE")==0) UPDATE=atoi(value);
	if (strcmp(key,"SNAPUPDATE")==0) SNAPUPDATE=atoi(value);
	if (strcmp(key,"AT")==0) AT=atof(value);
	if (strcmp(key,"EPS")==0) EPS=atof(value);
	if (strcmp(key,"ETAOS")==0) ETAOS=atof(value);
	if (strcmp(key,"COEFF")==0) COEFF=atof(value);
	if (strcmp(key,"TF")==0) TF=atof(value);
	if (strcmp(key,"TSTART")==0) TSTART=atof(value);
	if (strcmp(key,"IC")==0) IC=atof(value);
	if (strcmp(key,"PTASIZE")==0) PTASIZE=atoi(value);
	if (strcmp(key,"PHIPASIZE")==0) PHIPASIZE=atoi(value);
	if (strcmp(key,"FREEZE")==0) FREEZE=atoi(value);
	if (strcmp(key,"PTMAX")==0) PTMAX=atof(value);
	if (strcmp(key,"BIEPS")==0) BIEPS=atof(value);
	if (strcmp(key,"BIANGLE")==0) BIANGLE=atof(value);
	if (strcmp(key,"TRIEPS")==0) TRIEPS=atof(value);
	if (strcmp(key,"TRIANGLE")==0) TRIANGLE=atof(value);
	if (strcmp(key,"QUADEPS")==0) QUADEPS=atof(value);
	if (strcmp(key,"QUADANGLE")==0) QUADANGLE=atof(value);
	if (strcmp(key,"QUINTEPS")==0) QUINTEPS=atof(value);
	if (strcmp(key,"QUINTANGLE")==0) QUINTANGLE=atof(value);
	if (strcmp(key,"SEXEPS")==0) SEXEPS=atof(value);
	if (strcmp(key,"SEXANGLE")==0) SEXANGLE=atof(value);
	if (strcmp(key,"SEPTEPS")==0) SEPTEPS=atof(value);
	if (strcmp(key,"SEPTANGLE")==0) SEPTANGLE=atof(value);
	return;
}

//
// This routine assumes that paramters are in text file with
// each parameter on a new line in the format 
//
// PARAMKEY	PARAMVALUE
//
// The PARAMKEY must begin the line and only tabs and spaces
// can appear between the PARAMKEY and PARAMVALUE.
// 
// Lines which begin with 'commentmarker' defined below are ignored
//
void readParameters(const char *filename) {
		
	string commentmarker = "//"; 
	char space = ' '; 
	char tab = '\t';

	int maxline = 128; // maximum line length used in the buffer for reading
	char buffer[maxline];
	ifstream paramFile(filename);
	
	while(!paramFile.eof()) {
		paramFile.getline(buffer,maxline,'\n');
		string line = buffer; int length = strlen(buffer);
		if (line.substr(0,commentmarker.length())!=commentmarker && line.length()>0) {
			char key[32]="",value[32]=""; int founddelim=0;
			for (int i=0;i<length;i++) {
				if (buffer[i]==space || buffer[i]==tab) founddelim=1;
				else {
					if (founddelim==0) key[strlen(key)] = buffer[i];
					else value[strlen(value)] = buffer[i];
				}
			}
			if (strlen(key)>0 && strlen(value)>0) {
				setParameter(key,value);
				cout << key << " = " << value << endl;
			}
		}
	}
	
	return;	
}

