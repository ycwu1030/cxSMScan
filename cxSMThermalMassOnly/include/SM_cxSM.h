#ifndef CXSM_H
#define CXSM_H

#include "CubicSolver.h"
#include <stdio.h>
#include <vector>

#ifndef alpha1
#define alpha1 127.9
#endif

#ifndef MZ
#define MZ 91.1876
#endif

#ifndef MZ2
#define MZ2 (MZ*MZ)
#endif

#ifndef GF
#define GF 0.0000116637
#endif

#ifndef MT
#define MT 173.5
#endif

#ifndef MT2
#define MT2 (MT*MT)
#endif

#ifndef MHL
#define MHL 125.0
#endif

#ifndef MHL2
#define MHL2 (MHL*MHL)
#endif

class SM
{
public:
    SM();
    ~SM(){};

    double MW;
    double thetaW;  
    double alpha;
    double vev;     // (sqrt(2)GF)^(-0.5)
    double ee;      // ee = sqrt(4*Pi*alpha)
    double g_weak;  // g = ee/sw
    double gp_hyper; // g' = ee/cw
    double yt; // sqrt(2)mt/vev;

    double MW2;
    
};


class CXSM:public SM
{
public:
    CXSM();
    CXSM(double VSin, double MHHin, double MHAin, double thetain, double a1in);
    ~CXSM(){};

    void SetInput(double VSin=500, double MHHin=400, double MHAin=400, double thetain=0.1, double a1in=0);
    void SetInputZ2(double MHHin=400,double MHAin=450,double d2in=0.5, double del2in = 0.3);
// T=0 Potential and derivative
    double VT0(double phiH, double phiS);
    double dVT0dh(double phiH, double phiS);
    double dVT0ds(double phiH, double phiS);
    double d2VT0dh2(double phiH, double phiS);
    double d2VT0ds2(double phiH, double phiS);
    double d2VT0dhds(double phiH, double phiS);

//Hessian Matrix
    void GetHessian(double phiH, double phiS, double T);
    bool CheckHessianMatrix(double phiH, double phiS, double T);
    // Without T, they are used for LO potential
    void GetHessian(double phiH, double phiS);
    bool CheckHessianMatrix(double phiH, double phiS);
    
//Leading Order Stability and Unitarity
    bool CheckStability();
    bool CheckUnitarity(double MAX=0.5);
    bool CheckGlobalMinimum();

//Find LO local minimum position, used as starting points for critical point finding.
    void SolveCubicEquation(double A[4],double *results, int &NSolution);
    void FindLocalMinimumT0();
    void FindLocalMinimumT(double T);
    void PrintLocalMinimumT0();
    void PrintLocalMinimumT(double T);
    void PrintPotentialParameter();

// Potential Parameters:
    double mu2;
    double lam;
    double del2;
    double a1;
    double b1;
    double b2;
    double d2;

// Physical Parameters:
    double MHH;
    double MHH2;
    double MHA;
    double MHA2;
    double VS;
    double theta;

// Hessian Matrix:
    double HM11;
    double HM12;
    double HM22;

// Local minimum @ T=0
    CubicSolver Solver;
    int NLocalMinimaT0ZERO;
    int NLocalMinimaT0NONZERO;
    double LocalMinimumT0ZERO[3];
    double LocalMinimumT0NONZERO[3][2];
    int NLocalMinimaTZERO;
    int NLocalMinimaTNONZERO;
    double LocalMinimumTZERO[3];
    double LocalMinimumTNONZERO[3][2];

// Simple Thermal Potential with only thermal mass corrections.
    double Pih(double T);
    double dPihdT(double T);
    double Pis(double T);
    double dPisdT(double T);
    double VT(double phiH, double phiS, double T);
    double dVTdh(double phiH, double phiS, double T);
    double dVTds(double phiH, double phiS, double T);
    double dVTdT(double phiH, double phiS, double T);
    double d2VTdh2(double phiH, double phiS, double T);
    double d2VTds2(double phiH, double phiS, double T);
    double d2VTdhds(double phiH, double phiS, double T);
    double d2VTdhdT(double phiH, double phiS, double T);
    double d2VTdsdT(double phiH, double phiS, double T);
};


#endif
