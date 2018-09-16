/*-----------------------------------------------------------------------------
Compile With: "clang++ NM4P/rk4.cpp neuron.cpp 2neuron.cpp -o 2neuron"
-----------------------------------------------------------------------------*/

// we gotta add some threading for neuron 1 and neuron 2 eventually


#include <iostream>
#include "neuron.h"
#include "time.h"
#include <unistd.h> // timed delays only!
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <string.h>


using namespace std;

void timeStep(Time & T, Neuron & preSyn, Neuron & postSyn);
IType selectInput(char * arg);

const int N_DATA_POINTS = 10000;
const double TIME_STEP = .01;


double preI[N_DATA_POINTS];
double preVm[N_DATA_POINTS];

double preN[N_DATA_POINTS];
double preM[N_DATA_POINTS];
double preH[N_DATA_POINTS];

double preI_Na[N_DATA_POINTS];
double preI_K[N_DATA_POINTS];
double preI_Leak[N_DATA_POINTS];

double postI[N_DATA_POINTS];
double postVm[N_DATA_POINTS];

double timeAxis[N_DATA_POINTS];

int main(int argc, char *argv[])
{
    cout << "2NEURON | David Roccapriore | Spring 2017" << endl;

    Neuron n1("PRE");
    Neuron n2("POST");


    IType input = selectInput(argv[2]);
    cout << "Input Type: " << input << endl;
    n1.setInputCurrent(input);
    n1.setOutputCurrent(EPSP);

    n1.setPostSynapticNeuron(&n2); // couple neurons together


    Time T(TIME_STEP); // .01 ms per step
    for (int t = 0; t < N_DATA_POINTS; t++)
    {
        timeStep(T, n1, n2);
    }

    // print formatted csv data
    string filename = argv[1];
    ofstream cout(filename.c_str());
    for (int t = 0; t < N_DATA_POINTS; t++)
    {
        cout << timeAxis[t] << ",";
        cout << preI[t] << ",";
        cout << preVm[t] << ",";
        cout << preN[t] << ",";
        cout << preM[t] << ",";
        cout << preH[t] << ",";
        cout << preI_Na[t] << ",";
        cout << preI_K[t] << ",";
        cout << preI_Leak[t] << ",";        
        cout << postI[t] << ",";
        cout << postVm[t] << ",";
        cout << endl;
    }
    string cmd = "octave plotter.m " + filename;
    system(cmd.c_str());
}

void timeStep(Time & T, Neuron & preSyn, Neuron & postSyn)
{
    double t = T.getNeuronTime();
    double tau = T.getNeuronTimeStep();

    // T.logTime();


    double preSynPotential = preSyn.getPotential(t, tau);
    double postSynPotential = postSyn.getPotential(t, tau);


    int timeIndex = T.getTicks();
    timeAxis[timeIndex]     = t;
    preI[timeIndex]         = preSyn.getCurrent(t);
    preVm[timeIndex]        = preSynPotential;

    postI[timeIndex]         = postSyn.getCurrent(t);
    postVm[timeIndex]        = postSynPotential;

    preN[timeIndex]         = preSyn.getN();
    preM[timeIndex]         = preSyn.getM();
    preH[timeIndex]         = preSyn.getH();

    preI_Na[timeIndex]      = preSyn.getI_Na();
    preI_K[timeIndex]       = preSyn.getI_K();
    preI_Leak[timeIndex]    = preSyn.getI_Leak();

    // preSyn.logCurrentState(t);
    // postSyn.logCurrentState(t);

    // usleep(1000000);

    T.increment();
}


IType selectInput(char * arg)
{
    string x = arg;
    if (x == "ZERO"){return ZERO;}
    if (x == "CONST"){return CONST;}
    if (x == "PULSE"){return PULSE;}
    if (x == "STEP"){return STEP;}
    throw exception();
}