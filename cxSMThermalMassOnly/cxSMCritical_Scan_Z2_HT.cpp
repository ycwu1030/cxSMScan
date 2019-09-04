#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "SM_cxSM.h"
#include "Phases.h"
#include "Utilities.h"

using namespace std;
double RandomReal(double min, double max)
{
    return (rand()/(double)RAND_MAX)*(max-min)+min;
}
int main(int argc, char const *argv[])
{
    CXSM model = CXSM(500,450,450,0.1,1000.0);
    // model.FindLocalMinimumT0();
    // model.PrintLocalMinimumT0();
    double VS;
    double MHH,MHA;
    double d2;
    double del2;
    double theta;
    double a1;
    double phiSC,phiHB,phiSB,Tc;
    double phiSCd,phiHBd,phiSBd,Tcd;
    double results[4];
    // double resultsd[3];
    int status;
    // int statusd;
    double InitialGuess[3];
    bool GOOD;
    int id = atoi(argv[2]);
    char filename[500];
    sprintf(filename,"%s",argv[1]);
    ofstream output(filename);
    string logs;
    srand (time(NULL)+113*id);
    for (int i = 0; i < 1000 ; )
    {
        // VS = RandomReal(0,150);
        MHH = RandomReal(65,150);
        MHA = RandomReal(65,2000);
        d2 = RandomReal(0,20);
        del2 = RandomReal(-20,20);
        // theta = pow(10,RandomReal(-4.0,-1.0));
        // a1 = pow(RandomReal(-100,100),3);
        // model.SetInput(VS,MHH,MHA,theta,a1);
        model.SetInputZ2(MHH,MHA,d2,del2);
        if (!model.CheckUnitarity()||!model.CheckStability()||!model.CheckGlobalMinimum())
        {
            continue;
        }
        model.FindLocalMinimumT0();
        // GOOD = false;
        // for (int j = 0; j < model.NLocalMinimaT0ZERO; ++j)
        // {
        //     GOOD=GOOD||(model.CheckHessianMatrix(0.0,model.LocalMinimumT0ZERO[j])&&model.VT0(0.0,model.LocalMinimumT0ZERO[j])<0);
        // }
        // if (!GOOD)
        // {
        //     continue;
        // }
        status = FindCriticalPointsZ2(&model,results,logs,1000,true); 
        if (status != GSL_SUCCESS)
        {
            continue;
        }
        output<<i<<"  "<<MHH<<"  "<<MHA<<"  "<<d2<<"  "<<del2<<"  "<<logs;
        output<<endl;
        i++;
    }
    output.close();
    return 0;
}
