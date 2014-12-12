%{
#include "ghmm/bayesian_hmm.h"
%}

typedef struct ghmm_hyperparameters{
    ghmm_density_t type; //emission whose parameters we are estimating
    int num; //number of emissions for example for normal there are 2, normal and gamma

    union{
        ghmm_c_emission *continuous;//array of emissions 
        double *discrete;
    }emission;

}ghmm_hyperparameters;
%extend ghmm_hyperparameters{
    ghmm_hyperparameters(){ return calloc(1, sizeof(ghmm_hyperparameters)); }

    int alloc(ghmm_density_t type, int dim);
    ghmm_hyperparameters(ghmm_density_t type, int dim){
        ghmm_hyperparameters *hyper = calloc(1, sizeof(ghmm_hyperparameters));
        ghmm_hyperparameters_alloc(hyper, type, dim);
        return hyper;
    }

    void set_normal(double normal_mean,
        double normal_variance, double gamma_a, double gamma_b);

    void set_normal_truncated(
        double normal_lower_boundry, double normal_upper_boundry, 
        double normal_mean, double normal_variance,
        double gamma_lower_boundry, double gamma_upper_boundry, double gamma_a, double gamma_b);

    void sample_emission_discrete(double *emission, int dim);

    void sample_emission(ghmm_c_emission *emission);

    void set_discrete(double *counts, int dim);

    ~ghmm_hyperparameters(){ ghmm_hyperparameters_free(self); }
}


// All the hyperparameters for a model
typedef struct ghmm_bayes_hmm{
    //hyperparameters
    int continuous;//XXX 1 if continuous change later
    int dim;
    int N;
    int *M;
    ghmm_hyperparameters** params; //[state][mixture]
    double **A;
    double *pi;
}ghmm_bayes_hmm;
%extend ghmm_bayes_hmm{

ghmm_bayes_hmm(int dim, int n, int *m){
    ghmm_bayes_hmm *bay = calloc(1, sizeof(ghmm_bayes_hmm));
    bay->dim = dim;
    ghmm_bayes_hmm_alloc(bay, n, m);
    return bay;
}

void set_hyperparameter(int i, int j, ghmm_hyperparameters *p){
    self->params[i][j] = *p; 
}
    
int alloc(int N, int *M);

double* sample_initial();

double* sample_transistions(int i);

ghmm_dmodel* sample_model_discrete();

ghmm_cmodel* sample_model();

}
