
#include "neuron.h"
#include <iostream>
#include "NM4P/NumMeth.h"

using namespace std;

enum XType {N, M, H};
double alphaX(double u, XType x);
double betaX(double u, XType x);
double tauX(double u, XType x);
double initX(double u, XType x);
double getX(double u, XType x, double t);
void dVm_dt_RK(double x[], double t, double param[], double deriv[]);
void dnmh_dt_RK(double x[], double t, double param[], double deriv[]);
double constrain(double k, double min, double max)
{
    if (k < min){return min;}
    if (k > max){return max;}
    return k;
}


double alphaM(double V);
double betaM(double V);
double alphaH(double V);
double betaH(double V);
double alphaN(double V);
double betaN(double V);
double phi(double T);


Neuron::Neuron(string name)
{
    this->ID = name;
    this->Vm = -65; //mV
    m=.052;h=.596;n=.317;
    postSynapticPulseTime = 1;
}

void Neuron::logCurrentState(double t)
{
    cout << "Neuron: "  << this->ID << endl;
    cout << "Vm: "      << this->Vm << endl;
    cout << "I:  "      << this->I << endl;
    cout << "N:  "      << this->n << endl;
    cout << "M:  "      << this->m << endl;
    cout << "H:  "      << this->h << endl;

    cout << "I Na:    " << getI_Na() << endl;
    cout << "I K:     " << getI_K() << endl;
    cout << "I Leak:  " << getI_Leak() << endl;
}

void Neuron::setPostSynapticNeuron(Neuron * nextNeuron)
{
    this->postSynapticNeuron = nextNeuron;
    nextNeuron->setInputCurrent(this->outputCurrent);
}

void Neuron::setInputCurrent(IType input)
{
    this->inputCurrent = input;

}

void Neuron::setOutputCurrent(IType output)
{
    this->outputCurrent = output;
}

// These all take in a time t in milliseconds, the actual simulation time of neuron system
double Neuron::getCurrent(double t)
{
    double I = 0;

    switch(this->inputCurrent)
    {
        case ZERO:
            I = 0;
            break;
        case CONST:
            I = 10; // 20 microAmp / cm^2
            break;
        case PULSE:
            I = this->getCurrentPulse(t);
            break;
        case STEP:
            I = this->getCurrentStep(t);
            break;
        case EPSP:
            I = this->getCurrentEPSP(t);
            break;
        case IPSP:
            I = this->getCurrentIPSP(t);
            break;
    }
    return I;
}

double Neuron::getCurrentPulse(double t)
{
    double Iconst = 10;

    double off_interval = 15;
    double on_interval = 5;
    double modt = fmod(t, on_interval+off_interval);

    if (modt < on_interval)
    {
        return Iconst;
    }
    else {return -20;}
}
double Neuron::getCurrentStep(double t)
{
    const int dI = 2;
    const int dT = 7;
    int modt = fmod(t, 5*dT);

    if (modt < 1*dT){return 0*dI;}
    if (modt < 2*dT){return 1*dI;}
    if (modt < 3*dT){return 2*dI;}
    if (modt < 4*dT){return 3*dI;}
    return 4*dI;
}

double Neuron::getCurrentEPSP(double t)
{
    double Iconst = 10;

    double beginT = this->postSynapticPulseBegin;
    double endT = beginT + this->postSynapticPulseTime;

    if (t < endT)   {return Iconst;}
    else            {return -20;}
}
double Neuron::getCurrentIPSP(double t)
{
    double Iconst = -22;

    double beginT = this->postSynapticPulseBegin;
    double endT = beginT + this->postSynapticPulseTime;

    if (t < endT)   {return Iconst;}
    else            {return 10;}
}

double Neuron::getPotential(double t, double tau)
{
    // calculate N, M, H:
    double y [4] = {0, n, m, h};
    double param [4]; 
    param[0] = Vm;
    param[1] = n;
    param[2] = m;
    param[3] = h;

    rk4(y, 3, t, tau, &dnmh_dt_RK, param);
    n = constrain(y[1], 0, 1);
    m = constrain(y[2], 0, 1);
    h = constrain(y[3], 0, 1);

    // calculate Vm:
    double x [2] = {Vm, Vm};
    param[0] = I = this->getCurrent(t); // I
    param[1] = n;
    param[2] = m;
    param[3] = h;


    rk4(x, 1, t, tau, &dVm_dt_RK, param);

    Vm = constrain(x[1], -100, 100);


    // if Vm > 35, send EPSP to postsynaptic
    if (Vm > 35 && postSynapticNeuron != NULL)
    {
        postSynapticNeuron->postSynapticPulseBegin = t;
    }


    // cout << "New Vm: " << x[1] << endl;

    // experiment:
    // Vm = (g_Na*V_Na + g_K*V_K + g_Leak*V_Leak)/(g_Na + g_K + g_Leak);


    return Vm;
}

double alphaX(double u, XType x)
{
    double alpha;
    switch (x)
    {
        case N:
            alpha = (0.1-0.01*u)/(exp(1-(0.1*u))-1);
            break;
        case M:
            alpha = (2.5-0.1*u)/(exp(2.5-(0.1*u))-1);
            break;
        case H:
            alpha =  0.07*exp(-u/20);
            break;
    }
    return alpha;
}

double betaX(double u, XType x)
{
    double beta;
    switch (x)
    {
        case N:
            beta = 0.125*exp(-u/80);
            break;
        case M:
            beta = 4*exp(-u/18);
            break;
        case H:
            beta = 1/(exp(3-(0.1*u))+1);
            break;
    }
    return beta;
}

double tauX(double u, XType x)
{
    return pow((alphaX(u, x) + betaX(u, x)), 2);
}

double initX(double u, XType x)
{
    return alphaX(u, x)/(alphaX(u, x) + betaX(u, x));
}

double getX(double u, XType x, double t)
{
    // double x_0 = initX(u, x);

    // double dx_dt = (-1/tauX(u, M))*(x-x_0);
    // return dx_dt*t + x_0;
    switch(x)
    {
        case N : return N;
        case M : return M;
        case H : return H;
    }

}

void dVm_dt_RK(double x[], double t, double param[], double deriv[])
{
    double adjust = 32;//30;

    double Vm = x[0];
    double I = param[0];
    double n = param[1];
    double m = param[2];
    double h = param[3];

    deriv[0] = deriv[1] = ((I+adjust) - g_Na*h*(Vm-V_Na)*pow(m,3)-g_K*(Vm-V_K)*pow(n,4)-g_Leak*(Vm-V_Leak))/Cm;
}

void dnmh_dt_RK(double x[], double t, double param[], double deriv[])
{
    double V = param[0];
    double n = param[1];
    double m = param[2];
    double h = param[3];

    deriv[1] = alphaN(V)*(1-n)-betaN(V)*n; //n'
    deriv[2] = alphaM(V)*(1-m)-betaM(V)*m; //m'
    deriv[3] = alphaH(V)*(1-h)-betaH(V)*h; //h'
}

double phi(double T){return 1;}
double alphaM(double V){return phi(0)*.1*(V+40)/(1-exp(-(V+40)/10));}
double betaM(double V){return phi(0)*4*exp(-(V+65)/18);}
double alphaH(double V){return phi(0)*.07*exp(-(V+65)/20);}
double betaH(double V){return phi(0)*1/(1+exp(-(V+35)/10));}
double alphaN(double V){return phi(0)*.01*(V+55)/(1-exp(-(V+55)/10));}
double betaN(double V){return phi(0)*.125*exp(-(V+65)/80);}