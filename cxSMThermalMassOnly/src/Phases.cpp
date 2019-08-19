#include "Phases.h"

using namespace std;

void Phases::Clear()
{
    Temperatures.clear();
    Minima.clear();
}

void Phases::Set(double T, double LocalMinima[10][2], double NMinima)
{
    Phase current;
    for (int i = 0; i < NMinima; ++i)
    {
        current.push_back(make_pair(LocalMinima[i][0],LocalMinima[i][1]));
    }
    Temperatures.push_back(T);
    Minima.push_back(current);
}
void Phases::PrintLast()
{
    int len = Temperatures.size();
    printf("@Temperature:\t%f", Temperatures[len-1]);
    for (int i = 0; i < Minima[len-1].size(); ++i)
    {
        printf("\t (%f,%f)",Minima[len-1][i].first,Minima[len-1][i].second);
    }
    printf("\n");
}