#ifndef NEURON_H
#define NEURON_H

#include <vector>
#include <string>
#include <math.h>

using namespace std;

enum IType {ZERO, CONST, PULSE, STEP, EPSP, IPSP}; // EPSP and IPSP are time dependent functions based on on the presyn neuron

const double V_Na =     50;    // Current going through Sodium Channels (mV)
const double V_K =      -77;    // Current going through Potassium Channels (mV)
const double V_Leak =   -54.4;   // Current going through Leak Channels (mV)

const double g_Na =     120;//.5e-6;    // Sodium Channel Conductance (mS/cm^2)
const double g_K =      36;/// 10e-6;     // Potassium Channel Conductance (mS/cm^2)
const double g_Leak =   .3;///.10e-6;    // Leak Channel Conductance (mS/cm^2)

const double Cm = 1;            // Capacitance of cell (microF/cm2)

class Neuron
{
    private:
        string ID;      // an indentifier string

        double Vm;      // Membrane Potential
        double I;      // Membrane Current

        double n;
        double m;
        double h;
                
        IType inputCurrent;     // Type of current going into the Neuron
        IType outputCurrent;    // Type of current going into the Post Synaptic Neuron as a result of Action Potential from this neuron

        Neuron * postSynapticNeuron; // Pointer to next neuron

        bool pulseOn;

        double postSynapticPulseTime;   // the duration of a post synaptic 'current'
        double postSynapticPulseBegin;  // the start time of a post synaptic 'current'

        double getCurrentPulse(double t);
        double getCurrentStep(double t);
        double getCurrentEPSP(double t);
        double getCurrentIPSP(double t);

    public:
        Neuron(string name);
        void setPostSynapticNeuron(Neuron * nextNeuron);
        void setInputCurrent(IType input);
        void setOutputCurrent(IType output);
        void logCurrentState(double t);


        double getCurrent(double t); // given the time, returns the current in microAmp/cm^2
        double getPotential(double t, double tau);
        double getN(){return n;}
        double getM(){return m;}
        double getH(){return h;}

        double getI_Na(){return g_Na*(Vm-V_Na)*h*pow(m,3);}
        double getI_K(){return g_K*(Vm-V_K)*pow(n,4);}
        double getI_Leak(){return g_Leak*(Vm-V_Leak);}
};

#endif