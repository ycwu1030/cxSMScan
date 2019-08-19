#include <iostream>
#include <cmath>
#include "SM_cxSM.h"

using namespace std;

SM::SM()
{
    alpha = 1.0/alpha1;
    ee = sqrt(4*Pi*alpha);
    vev = pow(sqrt(2)*GF,-0.5);
    double A = sqrt(Pi*alpha)*vev;
    thetaW = asin(2*A/MZ)/2.0;
    MW = A/sin(thetaW);
    MW2 = MW*MW;
    g_weak = ee/sin(thetaW);
    gp_hyper = ee/cos(thetaW);
    yt = sqrt(2)*MT/vev;
#ifdef DEBUG
    cout<<"A:  "<<A<<endl;
    cout<<"vev: "<<vev<<endl;
    cout<<"thetaW: "<<thetaW<<endl;
    cout<<"MW: "<<MW<<endl;
    cout<<"yt: "<<yt<<endl;
    cout<<"g_weak: "<<g_weak<<endl;
    cout<<"gp_hyper: "<<gp_hyper<<endl;
#endif
}

CXSM::CXSM()
{
    SetInput();
}
CXSM::CXSM(double VSin, double MHHin, double MHAin, double thetain, double a1in)
{
    SetInput(VSin, MHHin, MHAin, thetain, a1in);
}
void CXSM::SetInput(double VSin, double MHHin, double MHAin, double thetain, double a1in)
{
    VS = VSin;
    MHH = MHHin;
    MHH2 = MHH*MHH;
    MHA = MHAin;
    MHA2 = MHA*MHA;
    theta = thetain;
    a1 = a1in;

    mu2 = (vev*cos(2*theta)*(MHH2-MHL2)-vev*(MHH2+MHL2)+VS*sin(2*theta)*(MHH2-MHL2))/(4*vev);
    b2 = -(4*sqrt(2)*a1-2*MHA2*VS-vev*sin(2*theta)*(MHH2-MHL2)+VS*cos(2*theta)*(MHH2-MHL2)+VS*(MHH2+MHL2))/(2*VS);
    lam = (cos(2*theta)*(MHL2-MHH2)+MHH2+MHL2)/(4*vev*vev);
    del2 = sin(2*theta)*(MHL2-MHH2)/vev/VS;
    b1 = -(sqrt(2)*a1+MHA2*VS)/VS;
    d2 = (2*sqrt(2)*a1+VS*cos(2*theta)*(MHH2-MHL2)+VS*(MHH2+MHL2))/(pow(VS,3));
}

void CXSM::SetInputZ2(double MHHin,double MHAin,double d2in, double del2in)
{
    VS = 0.0;
    MHH = MHHin;
    MHH2 = MHH*MHH;
    MHA = MHAin;
    MHA2 = MHA*MHA;
    theta = 0.0;
    a1 = 0.0;

    d2 = d2in;
    del2 = del2in;
    lam = MHL2/2.0/vev/vev;
    b1 = MHH2 - MHA2;
    b2 = MHH2 + MHA2 - del2/2.0*vev*vev;
    mu2 = -lam*vev*vev;
}

double CXSM::VT0(double phiH, double phiS)
{
    return sqrt(2)*a1*phiS + (b1+b2)/4.0*phiS*phiS + d2/16.0*pow(phiS,4) + mu2/2.0*phiH*phiH + lam/4.0*pow(phiH,4) + del2/8.0*pow(phiH*phiS,2);
}

double CXSM::dVT0dh(double phiH, double phiS)
{
    return phiH*(mu2+lam*phiH*phiH+del2*phiS*phiS/4.0);
}

double CXSM::dVT0ds(double phiH, double phiS)
{
    return (4.0*sqrt(2)*a1+phiS*(2.0*b1+2.0*b2+d2*phiS*phiS+del2*phiH*phiH))/4.0;
}

double CXSM::d2VT0dh2(double phiH, double phiS)
{
    return mu2 + 3.0*lam*phiH*phiH + del2/4.0*phiS*phiS;
}

double CXSM::d2VT0ds2(double phiH, double phiS)
{
    return (2.0*b1+2.0*b2+3.0*d2*phiS*phiS+del2*phiH*phiH)/4.0;
}

double CXSM::d2VT0dhds(double phiH, double phiS)
{
    return del2/2.0*phiH*phiS;
}

double CXSM::Pih(double T)
{
    return ((3*pow(g_weak,2)+pow(gp_hyper,2))/16 + yt*yt/4 + del2/24 + lam/2)*T*T;
}
double CXSM::dPihdT(double T)
{
    return 2*((3*pow(g_weak,2)+pow(gp_hyper,2))/16 + yt*yt/4 + del2/24 + lam/2)*T;
}
double CXSM::Pis(double T)
{
    return (d2+del2)/12*T*T;
}
double CXSM::dPisdT(double T)
{
    return (d2+del2)/6*T;
}
double CXSM::VT(double phiH, double phiS, double T)
{
    return VT0(phiH,phiS) + (Pih(T)*phiH*phiH+Pis(T)*phiS*phiS)/2.;
}

double CXSM::dVTdh(double phiH, double phiS, double T)
{
    return dVT0dh(phiH,phiS) + Pih(T)*phiH;
}

double CXSM::dVTds(double phiH, double phiS, double T)
{
    return dVT0ds(phiH,phiS) + Pis(T)*phiS;
}
double CXSM::dVTdT(double phiH, double phiS, double T)
{
    return (dPihdT(T)*phiH*phiH + dPisdT(T)*phiS*phiS)/2.0;
}
double CXSM::d2VTdh2(double phiH, double phiS, double T)
{
    return d2VT0dh2(phiH, phiS) + Pih(T);
}

double CXSM::d2VTds2(double phiH, double phiS, double T)
{
    return d2VT0ds2(phiH, phiS) + Pis(T);
}

double CXSM::d2VTdhds(double phiH, double phiS, double T)
{
    return d2VT0dhds(phiH, phiS);
}
double CXSM::d2VTdhdT(double phiH, double phiS, double T)
{
    return dPihdT(T)*phiH;
}
double CXSM::d2VTdsdT(double phiH, double phiS, double T)
{
    return dPisdT(T)*phiS;
}

bool CXSM::CheckStability()
{
    return (lam>0)&&(d2>0)&&((lam*d2>del2*del2)||del2>0);
}

bool CXSM::CheckUnitarity(double MAX)
{
    double EigenA0[4];
    EigenA0[0] = d2/2;
    EigenA0[1] = 2*lam;
    EigenA0[2] = (d2+6*lam-sqrt(d2*d2-12*d2*lam+2*del2*del2+36*lam*lam))/2;
    EigenA0[3] = (d2+6*lam+sqrt(d2*d2-12*d2*lam+2*del2*del2+36*lam*lam))/2;
    bool good=true;
    for (int i = 0; i < 4; ++i)
    {
        good*=(abs(EigenA0[i])<MAX*16.0*Pi);
        if (!good)
        {
            return good;
        }
    }

    return good;
}

void CXSM::PrintPotentialParameter()
{
    cout<<"Potential Parameter: "<<endl;
    cout<<"mu2: "<<mu2<<endl;
    cout<<"b2: "<<b2<<endl;
    cout<<"lam: "<<lam<<endl;
    cout<<"del2: "<<del2<<endl;
    cout<<"b1: "<<b1<<endl;
    cout<<"d2: "<<d2<<endl;
    cout<<"a1: "<<a1<<endl;
}

void CXSM::SolveCubicEquation(double A[4], double *results, int &NSolution)
{
    NSolution = 0;
    Solver.Solve(A[3],A[2],A[1],A[0]);
    double vSTemp;
    if (Solver.STATE==ONEREAL||Solver.STATE==THREEEQUALREAL)
    {
        vSTemp = Solver.SOLUTIONS[0];
        {
            results[NSolution]=vSTemp;
            NSolution++;
        }
    }
    else if (Solver.STATE==ONETWO)
    {
        for (int i = 0; i < 2; ++i)
        {
            vSTemp = Solver.SOLUTIONS[i];
            // if (CheckHessianMatrix(0.0,vSTemp,0.0))
            {
                results[NSolution]=vSTemp;
                NSolution++;
            }
        }
    }
    else if (Solver.STATE==THREEREAL)
    {
        for (int i = 0; i < 3; ++i)
        {
            vSTemp = Solver.SOLUTIONS[i];
            // if (CheckHessianMatrix(0.0,vSTemp,0.0))
            {
                results[NSolution]=vSTemp;
                NSolution++;
            }
        }
    }
    else
    {
#ifndef DEBUG
        std::cout<<"Warning in Cubic Solver 1"<<endl;
#endif
    }
}

void CXSM::FindLocalMinimumT0()
{
    NLocalMinimaT0ZERO=0; // Track we are getting which solution
    NLocalMinimaT0NONZERO=0;
    double results[3];
    // First solve the case with vH = 0
    double AA[4] = {4.0*sqrt(2)*a1,2.0*(b1+b2),0.0,d2};
    SolveCubicEquation(AA,results,NLocalMinimaT0ZERO);
    for (int i = 0; i < NLocalMinimaT0ZERO; ++i)
    {
        LocalMinimumT0ZERO[i] = results[i];
    }
    // Now For vH != 0 
    AA[3] = d2 - del2*del2/4.0/lam;
    AA[2] = 0.0;
    AA[1] = 2.0*(b1+b2)-del2*mu2/lam;
    AA[0] = 4.0*sqrt(2)*a1;
    int Ntemp=0;
    double vH2Temp;
    SolveCubicEquation(AA,results,Ntemp);
    for (int i = 0; i < Ntemp; ++i)
    {
        vH2Temp = -(mu2+del2/4.0*results[i]*results[i])/lam;
        if (vH2Temp>0)
        {
            LocalMinimumT0NONZERO[NLocalMinimaT0NONZERO][0] = sqrt(vH2Temp);
            LocalMinimumT0NONZERO[NLocalMinimaT0NONZERO][1] = results[i];
            NLocalMinimaT0NONZERO++;
        }
    }
}

void CXSM::FindLocalMinimumT(double T)
{
    NLocalMinimaTZERO=0; // Track we are getting which solution
    NLocalMinimaTNONZERO=0;
    double results[3];
    // First solve the case with vH = 0
    double AA[4] = {4.0*sqrt(2)*a1,2.0*(b1+b2)+4.0*Pis(T),0.0,d2};
    SolveCubicEquation(AA,results,NLocalMinimaTZERO);
    for (int i = 0; i < NLocalMinimaTZERO; ++i)
    {
        LocalMinimumTZERO[i] = results[i];
    }
    // Now For vH != 0 
    AA[3] = d2 - del2*del2/4.0/lam;
    AA[2] = 0.0;
    AA[1] = 2.0*(b1+b2) + 4.0*Pis(T) -del2*(mu2+Pih(T))/lam;
    AA[0] = 4.0*sqrt(2)*a1;
    int Ntemp=0;
    double vH2Temp;
    SolveCubicEquation(AA,results,Ntemp);
    for (int i = 0; i < Ntemp; ++i)
    {
        vH2Temp = -(mu2+Pih(T)+del2/4.0*results[i]*results[i])/lam;
        if (vH2Temp>0)
        {
            LocalMinimumTNONZERO[NLocalMinimaTNONZERO][0] = sqrt(vH2Temp);
            LocalMinimumTNONZERO[NLocalMinimaTNONZERO][1] = results[i];
            NLocalMinimaTNONZERO++;
        }
    }
}

void CXSM::PrintLocalMinimumT0()
{
    cout<<"The Local Minima for T0 are: "<<endl;
    bool Hessian;
    for (int i = 0; i < NLocalMinimaT0ZERO; ++i)
    {
        Hessian=CheckHessianMatrix(0,LocalMinimumT0ZERO[i]);
        if (!Hessian)
        {
            continue;
        }
        cout<<i+1<<"\tvH:\t"<<0.0<<"\tvS:\t"<<LocalMinimumT0ZERO[i]<<"\t\t"<<Hessian<<"\tV:"<<VT0(0,LocalMinimumT0ZERO[i])<<endl;
    }
    for (int i = 0; i < NLocalMinimaT0NONZERO; ++i)
    {
        Hessian=CheckHessianMatrix(LocalMinimumT0NONZERO[i][0],LocalMinimumT0NONZERO[i][1]);
        if (!Hessian)
        {
            continue;
        }
        cout<<i+1+NLocalMinimaT0ZERO<<"\tvH:\t"<<LocalMinimumT0NONZERO[i][0]<<"\tvS:\t"<<LocalMinimumT0NONZERO[i][1]<<"\t\t"<<Hessian<<"\tV:"<<VT0(LocalMinimumT0NONZERO[i][0],LocalMinimumT0NONZERO[i][1])<<endl;
    }

}
void CXSM::PrintLocalMinimumT(double T)
{
    cout<<"The Local Minima for T: "<<T<<" are: "<<endl;
    bool Hessian;
    for (int i = 0; i < NLocalMinimaTZERO; ++i)
    {
        Hessian=CheckHessianMatrix(0,LocalMinimumTZERO[i],T);
        if (!Hessian)
        {
            continue;
        }
        cout<<i+1<<"\tvH:\t"<<0.0<<"\tvS:\t"<<LocalMinimumTZERO[i]<<"\t\t"<<Hessian<<"\tV:"<<VT0(0,LocalMinimumTZERO[i])<<endl;
    }
    for (int i = 0; i < NLocalMinimaTNONZERO; ++i)
    {
        Hessian=CheckHessianMatrix(LocalMinimumTNONZERO[i][0],LocalMinimumTNONZERO[i][1],T);
        if (!Hessian)
        {
            continue;
        }
        cout<<i+1+NLocalMinimaTZERO<<"\tvH:\t"<<LocalMinimumTNONZERO[i][0]<<"\tvS:\t"<<LocalMinimumTNONZERO[i][1]<<"\t\t"<<Hessian<<"\tV:"<<VT0(LocalMinimumTNONZERO[i][0],LocalMinimumTNONZERO[i][1])<<endl;
    }

}

void CXSM::GetHessian(double phiH, double phiS, double T)
{
    HM11 = d2VTdh2(phiH,phiS,T);
    HM12 = d2VTdhds(phiH,phiS,T);
    HM22 = d2VTds2(phiH,phiS,T);
}

bool CXSM::CheckHessianMatrix(double phiH, double phiS, double T)
{
    GetHessian(phiH, phiS, T);
    return ((HM11+HM22)>0)&&(HM11*HM22>HM12*HM12);
}

void CXSM::GetHessian(double phiH, double phiS)
{
    HM11 = d2VT0dh2(phiH, phiS);
    HM12 = d2VT0dhds(phiH, phiS);
    HM22 = d2VT0ds2(phiH, phiS);
}

bool CXSM::CheckHessianMatrix(double phiH, double phiS)
{
    GetHessian(phiH, phiS);
    return ((HM11+HM22)>0)&&(HM11*HM22>HM12*HM12);
}

bool CXSM::CheckGlobalMinimum()
{
    FindLocalMinimumT0();
    double vinput = VT0(vev,VS);
    double vtemp;
    for (int i = 0; i < NLocalMinimaT0ZERO; ++i)
    {
        vtemp = VT0(0.0,LocalMinimumT0ZERO[i]);
        if (vtemp<vinput)
        {
            return false;
        }
    }
    for (int i = 0; i < NLocalMinimaT0NONZERO; ++i)
    {
        if (abs(vev - LocalMinimumT0NONZERO[i][0])<0.1&&abs(VS - LocalMinimumT0NONZERO[i][1]))
        {
            continue;
        }
        vtemp = VT0(LocalMinimumT0NONZERO[i][0],LocalMinimumT0NONZERO[i][1]);
        if (vtemp < vinput)
        {
            return false;
        }
    }
    return true;
}
