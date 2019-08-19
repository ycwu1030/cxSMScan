#ifndef Utilities_H
#define Utilities_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <string>
#include "SM_cxSM.h"
#include "Phases.h"

double MAX(double x[3],int Number);
// For finding the critical points
int CriticalF(const gsl_vector * x, void *model, gsl_vector * f);
int CriticalDF(const gsl_vector * x, void *model, gsl_matrix *J);
int CriticalFDF(const gsl_vector * x, void *model, gsl_vector * f, gsl_matrix * J);
// int FindCriticalPointsSingle(CXSM *model, double *res, int MAXITER=1000, bool usingdf=true);
int FindCriticalPoints(CXSM *model, double *res, std::string &logs, int MAXITER=1000, bool usingdf=true);
int FindCriticalPointsZ2(CXSM *model, double *res, int MAXITER=1000, bool usingdf=false);
// For finding the minimum
struct MinParam
{
    CXSM *model;
    double T;
};
bool CloseToEachOther(double a, double b);
double MinF(const gsl_vector * x, void *param);
void MinDF(const gsl_vector * x, void *param, gsl_vector *df);
void MinFDF(const gsl_vector * x, void *param, double *f, gsl_vector * df);
int FindPhasesSingle(CXSM *model, double T, double StartingPoints[2], double Results[2], int MAXITER, bool usingdf);
int FindPhases(CXSM *model, double T, Phases *modelPhase, int MAXITER=500, bool usingdf=true);
int FindPhasesZ2(CXSM *model, double T, Phases *modelPhase, int MAXITER, bool usingdf);
#endif