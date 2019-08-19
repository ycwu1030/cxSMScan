#ifndef PHASES_H
#define PHASES_H

#include <vector>
#include <stdio.h>


typedef std::vector<std::pair<double, double> > Phase;

class Phases
{
public:
    Phases(){};
    ~Phases(){};
    
    void Clear();
    void PrintLast();
    void Set(double T, double LocalMinima[10][2], double NMinima);
    std::vector<double> Temperatures;
    std::vector<Phase> Minima;
};



#endif
