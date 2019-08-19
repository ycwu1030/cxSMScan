#include "Utilities.h"
double MAX(double x[3],int Number)
{
    if (Number == 1)
    {
        return x[0];
    }
    else if (Number == 2)
    {
        return x[0]>x[1]?x[0]:x[1];
    }
    else
    {
        double temp = x[0];
        if (temp < x[1]) temp = x[1];
        if (temp < x[2]) temp = x[2];
        return temp;
    }
}

int CriticalF(const gsl_vector *x, void *model, gsl_vector * f)
{
    CXSM *mymodel = (CXSM *) model;
    const double phiSC = gsl_vector_get(x,0);
    const double phiHB = gsl_vector_get(x,1);
    const double phiSB = gsl_vector_get(x,2);
    const double Tc = gsl_vector_get(x,3);

    const double y0 = mymodel->dVTds(0.0,phiSC,Tc);
    const double y1 = mymodel->dVTdh(phiHB,phiSB,Tc);
    const double y2 = mymodel->dVTds(phiHB,phiSB,Tc);
    const double y3 = mymodel->VT(0.0,phiSC,Tc)-mymodel->VT(phiHB,phiSB,Tc);
    gsl_vector_set(f,0,y0);
    gsl_vector_set(f,1,y1);
    gsl_vector_set(f,2,y2);
    gsl_vector_set(f,3,y3);

    return GSL_SUCCESS;
}
int CriticalDF(const gsl_vector *x, void *model, gsl_matrix *J)
{
    CXSM *mymodel = (CXSM *) model;
    const double phiSC = gsl_vector_get(x,0);
    const double phiHB = gsl_vector_get(x,1);
    const double phiSB = gsl_vector_get(x,2);
    const double Tc = gsl_vector_get(x,3);

    const double df00 = mymodel->d2VTds2(0.0,phiSC,Tc);
    const double df01 = 0.0;
    const double df02 = 0.0;
    const double df03 = mymodel->d2VTdsdT(0.0,phiSC,Tc);

    const double df10 = 0.0;
    const double df11 = mymodel->d2VTdh2(phiHB,phiSB,Tc);
    const double df12 = mymodel->d2VTdhds(phiHB,phiSB,Tc);
    const double df13 = mymodel->d2VTdhdT(phiHB,phiSB,Tc);

    const double df20 = 0.0;
    const double df21 = mymodel->d2VTdhds(phiHB,phiSB,Tc);
    const double df22 = mymodel->d2VTds2(phiHB,phiSB,Tc);
    const double df23 = mymodel->d2VTdsdT(phiHB,phiSB,Tc);

    const double df30 = mymodel->dVTds(0.0,phiSC,Tc);
    const double df31 = mymodel->dVTdh(phiHB,phiSB,Tc);
    const double df32 = mymodel->dVTds(phiHB,phiSB,Tc);
    const double df33 = mymodel->dVTdT(0.0,phiSC,Tc)-mymodel->dVTdT(phiHB,phiSB,Tc);

    gsl_matrix_set (J, 0, 0, df00);
    gsl_matrix_set (J, 0, 1, df01);
    gsl_matrix_set (J, 0, 2, df02);
    gsl_matrix_set (J, 0, 3, df03);

    gsl_matrix_set (J, 1, 0, df10);
    gsl_matrix_set (J, 1, 1, df11);
    gsl_matrix_set (J, 1, 2, df12);
    gsl_matrix_set (J, 1, 3, df13);

    gsl_matrix_set (J, 2, 0, df20);
    gsl_matrix_set (J, 2, 1, df21);
    gsl_matrix_set (J, 2, 2, df22);
    gsl_matrix_set (J, 2, 3, df23);

    gsl_matrix_set (J, 3, 0, df30);
    gsl_matrix_set (J, 3, 1, df31);
    gsl_matrix_set (J, 3, 2, df32);
    gsl_matrix_set (J, 3, 3, df33);

    return GSL_SUCCESS;
}
int CriticalFDF(const gsl_vector *x, void *model, gsl_vector * f, gsl_matrix * J)
{
    CriticalF(x,model,f);
    CriticalDF(x,model,J);
    return GSL_SUCCESS;
}
void print_state(size_t iter, gsl_multiroot_fdfsolver * s)
{
  printf ("iter = %3u x = %.3f %.3f %.3f %.3f "
          "f(x) = % .3e % .3e %.3e %.3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2),
          gsl_vector_get (s->f, 3)
          );
}
void print_state(size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = %.3f %.3f %.3f %.3f "
          "f(x) = % .3e % .3e %.3e %.3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2),
          gsl_vector_get (s->f, 3)
          );
}
int MultiRootIterFDF(CXSM *model, double *InitialGuess, double *res, int MAXITER)
{
    const size_t n = 4;
    int status;
    size_t iter = 0;
    gsl_vector * x = gsl_vector_alloc(n);
    for (int i = 0; i < n; ++i)
    {
        gsl_vector_set(x, i, InitialGuess[i]);
    }
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;
    gsl_multiroot_function_fdf FF = {&CriticalF,&CriticalDF,&CriticalFDF,n,model};
    // T = gsl_multiroot_fdfsolver_gnewton;
    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc(T,n);
    gsl_multiroot_fdfsolver_set(s, &FF, x);
#ifdef DEBUG
    print_state(iter,s);
#endif
    do{
        iter++;
        status = gsl_multiroot_fdfsolver_iterate(s);
#ifdef DEBUG
        print_state(iter,s);
#endif
        if (status) break;
        status = gsl_multiroot_test_residual(s->f,1e-6);
    } while (status == GSL_CONTINUE && iter < MAXITER);
#ifdef DEBUG
    print_state(iter,s);
#endif
    for (int i = 0; i < n; ++i)
    {
        res[i] = gsl_vector_get(s->x,i);
    }
    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x);  
    return status;
}
int MultiRootIterF(CXSM *model, double *InitialGuess, double *res, int MAXITER)
{
    const size_t n = 4;
    int status;
    size_t iter = 0;
    gsl_vector * x = gsl_vector_alloc(n);
    for (int i = 0; i < n; ++i)
    {
        gsl_vector_set(x, i, InitialGuess[i]);
    }
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    gsl_multiroot_function FF = {&CriticalF,n,model};
    // T = gsl_multiroot_fdfsolver_gnewton;
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,n);
    gsl_multiroot_fsolver_set(s, &FF, x);
#ifdef DEBUG
    print_state(iter,s);
#endif
    do{
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        // status = gsl_multiroot_fsolver_iterate(s);
#ifdef DEBUG
        print_state(iter,s);
#endif
        if (status) break;
        status = gsl_multiroot_test_residual(s->f,1e-6);
    } while (status == GSL_CONTINUE && iter < MAXITER);
#ifdef DEBUG
    print_state(iter,s);
#endif
    for (int i = 0; i < n; ++i)
    {
        res[i] = gsl_vector_get(s->x,i);
    }
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);  
    return status;
}
int FindCriticalPointsSingle(CXSM *model, double *InitialGuess, double *res, int MAXITER, bool usingdf)
{
    if (usingdf)
    {
        return MultiRootIterFDF(model, InitialGuess, res, MAXITER);
    }
    else
    {
        return MultiRootIterF(model, InitialGuess, res, MAXITER);
    }
}
int FindCriticalPoints(CXSM *model, double *res, std::string &logs,int MAXITER, bool usingdf)
{
    double InitialGuess[4];
    int TotalStatus = GSL_CONTINUE;
    int status;
    double phiSC,phiHB,phiSB,Tc;
    logs = "";
    char temp[200];
    for (int izero = 0; izero < model->NLocalMinimaT0ZERO; ++izero)
    {
        InitialGuess[0] = model->LocalMinimumT0ZERO[izero];
        InitialGuess[1] = model->vev;
        InitialGuess[2] = model->VS;
        InitialGuess[3] = 100.0;
        status = FindCriticalPointsSingle(model, InitialGuess, res, MAXITER, usingdf);
        if (status != GSL_SUCCESS)
        {
            continue;
        }
        phiSC = res[0];
        phiHB = res[1];
        phiSB = res[2];
        Tc = res[3];
        // if (CloseToEachOther(0.0,phiHB)&&CloseToEachOther(phiSC,phiSB))
        if (CloseToEachOther(0.0,phiHB))
        {
            continue;
        }
        if (Tc<0 || !model->CheckHessianMatrix(0.0,phiSC,Tc) || !model->CheckHessianMatrix(phiHB,phiSB,Tc))
        {
            continue;
        }
        sprintf(temp,"STATUS: %d; Tc: %f, From (0.0,%f)->(%f,%f)",status,Tc,phiSC,phiHB,phiSB);
        logs += temp;
        TotalStatus = GSL_SUCCESS;
    }
    return TotalStatus;
}
// For finding the Z2 case minimum
double MinZ2F(const gsl_vector * x, void *param)
{
    CXSM *mymodel = ((MinParam *) param)->model;
    double Temp = ((MinParam *) param)->T;
    const double phiHB = gsl_vector_get(x,0);
    const double phiSB = gsl_vector_get(x,1);

    return mymodel->VT(phiHB,phiSB,Temp);
}
void MinZ2DF(const gsl_vector * x, void *param, gsl_vector *df)
{
    CXSM *mymodel = ((MinParam *) param)->model;
    double Temp = ((MinParam *) param)->T;
    const double phiHB = gsl_vector_get(x,0);
    const double phiSB = gsl_vector_get(x,1);
    // mymodel->SetUpMassDerivatives(phiHB,phiSB);
    gsl_vector_set(df, 0, mymodel->dVTdh(phiHB,phiSB,Temp));
    gsl_vector_set(df, 1, mymodel->dVTds(phiHB,phiSB,Temp));
}
void MinZ2FDF(const gsl_vector * x, void *param, double *f, gsl_vector * df)
{
    *f = MinZ2F(x, param);
    MinZ2DF(x,param,df);
}

int FindPhasesSingleZ2(CXSM *model, double T, double StartingPoints[2], double Results[2], int MAXITER, bool usingdf)
{
    MinParam p = {model, T};
    size_t iter = 0;
    int status;
    double size;
    gsl_vector *x;
    x = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, StartingPoints[0]);
    gsl_vector_set (x, 1, StartingPoints[1]);
    if (usingdf)
    {
        const gsl_multimin_fdfminimizer_type *T;
        gsl_multimin_fdfminimizer *s;
        gsl_multimin_function_fdf my_func;
        my_func.n = 2;
        my_func.f = MinZ2F;
        my_func.df = MinZ2DF;
        my_func.fdf = MinZ2FDF;
        my_func.params = &p;
        T = gsl_multimin_fdfminimizer_conjugate_fr;
        s = gsl_multimin_fdfminimizer_alloc (T, 2);
        gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 0.1);
        do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate (s);
            if (status)
                break;
            status = gsl_multimin_test_gradient (s->gradient, 1e-3);
        }while (status == GSL_CONTINUE && iter < MAXITER);
        if (status == GSL_SUCCESS)
        {
            Results[0] = gsl_vector_get(s->x, 0);
            Results[1] = gsl_vector_get(s->x, 1);
        }
        gsl_multimin_fdfminimizer_free (s);
    }
    else
    {
        const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
        gsl_multimin_fminimizer *s = NULL;
        gsl_vector *ss;
        gsl_multimin_function my_func;
        my_func.n = 2;
        my_func.f = MinZ2F;
        my_func.params = &p;
        ss = gsl_vector_alloc(2);
        gsl_vector_set_all (ss, 1.0);
        s = gsl_multimin_fminimizer_alloc (T, 2);
        gsl_multimin_fminimizer_set (s, &my_func, x, ss);
        do
        {
            iter++;
            status = gsl_multimin_fminimizer_iterate (s);
            if (status)
                break;
            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 0.1);
        }while (status == GSL_CONTINUE && iter < MAXITER);
        if (status == GSL_SUCCESS)
        {
            Results[0] = gsl_vector_get(s->x, 0);
            Results[1] = gsl_vector_get(s->x, 1);
        }
        gsl_vector_free(ss);
        gsl_multimin_fminimizer_free (s);
    }

    gsl_vector_free (x);
    return status;

}
bool CloseToEachOther(double a, double b)
{
    if (abs(a-b)<1.0)
    {
        return true;
    }
    return false;
}

bool DuplicatedMinimum(double results[2],double TMinimum[6][2],int NTMinima)
{
    for (int i = 0; i < NTMinima; ++i)
    {
        if (CloseToEachOther(results[0],TMinimum[i][0])&&CloseToEachOther(results[1],TMinimum[i][1]))
        {
            return true;
        }
    }
    return false;
}
int FindPhasesZ2(CXSM *model, double T, Phases *modelPhase, int MAXITER, bool usingdf)
{
    model->FindLocalMinimumT0();
    double T0Minimum[6][2];
    int NT0Minima;
    double TMinimum[6][2];
    int NTMinima=0;
    int status;
    double results[2];
    NT0Minima = model->NLocalMinimaT0ZERO + model->NLocalMinimaT0NONZERO;
    for (int i = 0; i < model->NLocalMinimaT0ZERO; ++i)
    {
        T0Minimum[i][0] = 0.0;
        T0Minimum[i][1] = model->LocalMinimumT0ZERO[i];
    }
    for (int i = 0; i < model->NLocalMinimaT0NONZERO; ++i)
    {
        T0Minimum[i+model->NLocalMinimaT0ZERO][0] = model->LocalMinimumT0NONZERO[i][0];
        T0Minimum[i+model->NLocalMinimaT0ZERO][1] = model->LocalMinimumT0NONZERO[i][1];
    }
    for (int i = 0; i < NT0Minima; ++i)
    {
        if (T0Minimum[i][0]<-20.0||T0Minimum[i][1]<-20.0)
        {
            continue;
        }
        status = FindPhasesSingleZ2(model,T,T0Minimum[i],results,MAXITER,usingdf);
        if (status == GSL_SUCCESS)
        {
            if(!DuplicatedMinimum(results,TMinimum,NTMinima)&&results[0]>-20&&results[1]>-20)
            {
                TMinimum[NTMinima][0] = results[0];
                TMinimum[NTMinima][1] = results[1];
                NTMinima++;
            }
        }
    }
    modelPhase->Set(T,TMinimum,NTMinima);
    return NTMinima;
}

int CriticalZ2F(const gsl_vector * x, void *model, gsl_vector * f)
{
    CXSM *mymodel = (CXSM *) model;
    const double phiSC = gsl_vector_get(x,0);
    const double phiHB = gsl_vector_get(x,1);
    const double T = gsl_vector_get(x,2);

    const double y0 = mymodel->dVTds(0,phiSC,T);
    const double y1 = mymodel->dVTdh(phiHB,0,T);
    const double y2 = mymodel->VT(0,phiSC,T)-mymodel->VT(phiHB,0,T);
    gsl_vector_set(f,0,y0);
    gsl_vector_set(f,1,y1);
    gsl_vector_set(f,2,y2);

    return GSL_SUCCESS;
}
int CriticalZ2DF(const gsl_vector * x, void *model, gsl_matrix *J)
{
    CXSM *mymodel = (CXSM *) model;
    const double phiSC = gsl_vector_get(x,0);
    const double phiHB = gsl_vector_get(x,1);
    const double T = gsl_vector_get(x,2);

    const double df00 = mymodel->d2VTds2(1e-12,phiSC,T);
    const double df01 = 1e-12;
    const double df02 = mymodel->dPisdT(T)*phiSC;
    const double df10 = 1e-12;
    const double df11 = mymodel->d2VTdh2(phiHB,1e-12,T);
    const double df12 = mymodel->dPihdT(T)*phiHB;
    const double df20 = mymodel->dVTds(1e-12,phiSC,T);
    const double df21 = mymodel->dVTdh(phiHB,0,T);
    const double df22 = (mymodel->dPisdT(T)*phiSC*phiSC/2.0)-(mymodel->dPihdT(T)*phiHB*phiHB)/2.0;

    gsl_matrix_set (J, 0, 0, df00);
    gsl_matrix_set (J, 0, 1, df01);
    gsl_matrix_set (J, 0, 2, df02);
    gsl_matrix_set (J, 1, 0, df10);
    gsl_matrix_set (J, 1, 1, df11);
    gsl_matrix_set (J, 1, 2, df12);
    gsl_matrix_set (J, 2, 0, df20);
    gsl_matrix_set (J, 2, 1, df21);
    gsl_matrix_set (J, 2, 2, df22);
    return GSL_SUCCESS;
}
int CriticalZ2FDF(const gsl_vector * x, void *model, gsl_vector * f, gsl_matrix * J)
{
    CriticalZ2F(x, model, f);
    CriticalZ2DF(x, model, J);

    return GSL_SUCCESS;
}
void print_state_z2 (size_t iter, gsl_multiroot_fdfsolver * s)
{
  printf ("iter = %3u x = %.3f %.3f %.3f "
          "f(x) = % .3e % .3e %.3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          // gsl_vector_get (s->x, 3),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2));
}
int FindCriticalPointsZ2(CXSM *model, double *res, int MAXITER, bool usingdf)
{
    const size_t n = 3;
    int status;
    size_t i, iter = 0;
    double phiHC,phiSC,phiHB,phiSB,Tc;
    double InitialGuess[n];
    gsl_vector * x = gsl_vector_alloc(n);
    model->FindLocalMinimumT0();
    InitialGuess[1]=model->vev;
    InitialGuess[0]=MAX(model->LocalMinimumT0ZERO,model->NLocalMinimaT0ZERO);
    InitialGuess[2]=100.0;
    for (int j = 0; j < n; ++j)
    {
        gsl_vector_set(x,j,InitialGuess[j]);
    }
    if (usingdf)
    {
        const gsl_multiroot_fdfsolver_type *T;
        gsl_multiroot_fdfsolver *s;
        gsl_multiroot_function_fdf FF = {&CriticalZ2F,&CriticalZ2DF,&CriticalZ2FDF,n,model};
        iter = 0;
        T = gsl_multiroot_fdfsolver_gnewton;
        s = gsl_multiroot_fdfsolver_alloc(T,n);
        gsl_multiroot_fdfsolver_set(s, &FF, x);
#ifdef DEBUG
        print_state_df(iter,s);
#endif
        do{
            iter++;
            status = gsl_multiroot_fdfsolver_iterate(s);
            // status = gsl_multiroot_fsolver_iterate(s);
#ifdef DEBUG
            print_state_df(iter,s);
#endif
            if (status) break;
            status = gsl_multiroot_test_residual(s->f,n*1.0);
        } while (status == GSL_CONTINUE && iter < MAXITER);
#ifdef DEBUG
        print_state_df(iter,s);
#endif
        for (int i = 0; i < n; ++i)
        {
            res[i] = gsl_vector_get(s->x,i);
        }
        gsl_multiroot_fdfsolver_free(s);
    }
    else
    {
        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;
        gsl_multiroot_function FF = {&CriticalZ2F,n,model};
        T = gsl_multiroot_fsolver_dnewton;
        s = gsl_multiroot_fsolver_alloc(T,n);
        gsl_multiroot_fsolver_set(s, &FF, x);
        iter = 0;
#ifdef DEBUG
        print_state(iter,s);
#endif
        do{
            iter++;
            status = gsl_multiroot_fsolver_iterate(s);
            // status = gsl_multiroot_fsolver_iterate(s);
#ifdef DEBUG
            print_state(iter,s);
#endif
            if (status) break;
            status = gsl_multiroot_test_residual(s->f,n*1.0);
        } while (status == GSL_CONTINUE && iter < MAXITER);
#ifdef DEBUG
        print_state(iter,s);
#endif
        for (int i = 0; i < n; ++i)
        {
            res[i] = gsl_vector_get(s->x,i);
        }
        gsl_multiroot_fsolver_free(s);
    }
    gsl_vector_free(x);  
    return status;
}
