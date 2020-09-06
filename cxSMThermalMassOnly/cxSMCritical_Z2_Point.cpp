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
    string logs;

        // VS = RandomReal(0,150);
        MHH = 129.208;
        MHA = 949.971;
        d2 = 9.46504;
        del2 = 2.7516599999999998;
        model.SetInputZ2(MHH,MHA,d2,del2);
        if (!model.CheckUnitarity()||!model.CheckStability()||!model.CheckGlobalMinimum())
        {
            cout<<"Not Good"<<endl;
            return 0;
        }
        model.FindLocalMinimumT0();
        model.PrintLocalMinimumT0();
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
            return 0;
        }
        cout<<"  "<<MHH<<"  "<<MHA<<"  "<<d2<<"  "<<del2<<"  "<<logs<<endl;
    return 0;
}
