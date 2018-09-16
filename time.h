#ifndef TIME_H
#define TIME_H

#include <iostream>
#include <time.h>


using namespace std;

class Time
{
    private:
        int tickTime;           // Number of 'ticks' incremented during each for loop 
        double neuronTimeStep;  // The amount of time between ticks, in milliseconds (?)
        double neuronTime;      // Simulation Time, measured in milliseconds (?)
        double computerTime;    // Real Time, measured in seconds
        time_t startTime;



    public: 
        Time(double timeStep)
        {
            neuronTimeStep = timeStep;
            resetTime();
        }
        void increment()
        {   tickTime++;
            neuronTime = tickTime*neuronTimeStep;
            computerTime = difftime(time(NULL), startTime); //real time at BEGINNING of increment, not end!
        }
        int getTicks(){return tickTime;}
        double getNeuronTime(){return neuronTime;}
        double getComputerTime(){return computerTime;}
        void resetTime()
        {   tickTime = 0;
            neuronTime = 0;
            time(&startTime); // update the computer time for t = 0
        }
        void logTime()
        {
            cout << endl << "Tick Time:     " << tickTime <<  endl;
            cout << "Neuron Time:   " << neuronTime << " ms" << endl;
            cout << "Computer Time: " << computerTime << " s" << endl;
        }
        double getNeuronTimeStep(){return neuronTimeStep;}
};

#endif