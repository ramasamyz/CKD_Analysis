/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/model.c
*       Authors:  Benhard Knab, Bernd Wichern, Benjamin Georgi, Alexander Schliep
*
*       Copyright (C) 1998-2004 Alexander Schliep
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
*                               Berlin
*
*       Contact: schliep@ghmm.org
*
*       This library is free software; you can redistribute it and/or
*       modify it under the terms of the GNU Library General Public
*       License as published by the Free Software Foundation; either
*       version 2 of the License, or (at your option) any later version.
*
*       This library is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*       Library General Public License for more details.
*
*       You should have received a copy of the GNU Library General Public
*       License along with this library; if not, write to the Free
*       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*
*       This file is version $Revision: 2304 $
*                       from $Date: 2013-05-31 13:48:13 -0400 (Fri, 31 May 2013) $
*             last change by $Author: ejb177 $.
*
*******************************************************************************/

#include "bayesian_hmm.h"
#include "mes.h"
#include "randvar.h"
#include "math.h"

int ghmm_hyperparameters_alloc(ghmm_hyperparameters *params, ghmm_density_t type, int dim){
#define CUR_PROC "alloc_hyperparameters"
    switch(type){
        //postpone boundries until set
        case(normal):
            //normal-gamma (mean, var, a, b) b = scale
            params->type = type;
            params->num = 2;
            ARRAY_CALLOC(params->emission.continuous,2);
            params->emission.continuous[0].type = normal;
            params->emission.continuous[1].type = gamma_density;
            return 0;
        case(discrete):
            ARRAY_CALLOC(params->emission.discrete, 1);
            params->type = type;
            params->num = dim;
            return 0;
        default:
            goto STOP;//not supported
    }
STOP:
   return -1;
#undef CUR_PROC
}

//use INFINITY for no boundries
//XXX need to support sampling from truncated normals in randvar.c
void ghmm_hyperparameters_set_normal_truncated(ghmm_hyperparameters *params,
        double normal_lower_boundry, double normal_upper_boundry, 
        double normal_mean, double normal_variance,
        double gamma_lower_boundry, double gamma_upper_boundry, double gamma_a, double gamma_b){

    ghmm_hyperparameters_set_normal(params, normal_mean, normal_variance, 
        gamma_a, gamma_b);

    //mean boundries 
    if(normal_lower_boundry == -INFINITY && normal_upper_boundry != INFINITY){

        params->emission.continuous[0].max = normal_upper_boundry;
        params->emission.continuous[0].type = normal_left;
    }
    else if(normal_lower_boundry != -INFINITY && normal_upper_boundry == INFINITY){

        params->emission.continuous[0].min = normal_lower_boundry;
        params->emission.continuous[0].type = normal_right;
    }
    else if (normal_lower_boundry == -INFINITY && normal_upper_boundry == INFINITY) {
        params->emission.continuous[0].type = normal;
    }
    else{
        params->emission.continuous[0].type = truncated_normal;
        params->emission.continuous[0].min = normal_lower_boundry;
        params->emission.continuous[0].max = normal_upper_boundry;
    }

    //variance boundries
    if(gamma_lower_boundry == INFINITY && gamma_upper_boundry == INFINITY){
        params->emission.continuous[1].type = gamma_density;
    }
    else{
        params->emission.continuous[1].min = gamma_lower_boundry;
        params->emission.continuous[1].max = gamma_upper_boundry;
        params->emission.continuous[1].type = gamma_truncated;
    }

}

void ghmm_hyperparameters_set_normal(ghmm_hyperparameters *params, double normal_mean,
        double normal_variance, double gamma_a, double gamma_b){
    params->type = normal;
    params->emission.continuous[0].mean.val = normal_mean;
    params->emission.continuous[0].variance.val = normal_variance;
    params->emission.continuous[1].alpha = gamma_a;
    params->emission.continuous[1].beta = gamma_b;
}

void ghmm_hyperparameters_set_discrete(ghmm_hyperparameters *params, double *counts, int dim){
#define CUR_PROC "ghmm_set_discrete_hyperparamters"
    int i;
    ARRAY_MALLOC(params->emission.discrete, dim);
    for(i=0;i<dim;i++){
        params->emission.discrete[i] = counts[i];
    }
STOP: return;
#undef CUR_PROC
}

ghmm_dmodel* ghmm_bayes_hmm_sample_model_discrete(ghmm_bayes_hmm *mo){
#define CUR_PROC "ghmm_sample_model"
    //XXX Mixtures, 1 transitions, multivariate
    int i,j;
    ghmm_dmodel* hmm;
    int *inDeg;
    int *outDeg;
    ARRAY_MALLOC(inDeg, mo->N);
    ARRAY_MALLOC(outDeg, mo->N);
    for(i = 0; i < mo->N; i++){
        inDeg[i] = outDeg[i] = mo->N;
    }
    hmm = ghmm_dmodel_calloc(mo->dim, mo->N, GHMM_kDiscreteHMM, inDeg, outDeg);
   
    ghmm_dstate *s =  hmm->s; 
    double **tmp;
    ARRAY_MALLOC(tmp, mo->N);

    //get sample mats for transition
    for(i = 0; i < mo->N; i++){
        tmp[i] = ghmm_bayes_hmm_sample_transistions(mo, i);
    }

    //get in/out state transitions
    int out_state;
    int in_state;
    for(i=0; i<mo->N; i++){
        for(j=0; j<mo->N; j++){
            if(tmp[i][j] != 0){
                out_state = s[i].out_states++;
                in_state  = s[j].in_states++;
                s[i].out_id[out_state] = j;
                s[i].out_a[out_state] = tmp[i][j];
                s[j].in_id[in_state] = i;
                s[j].in_a[in_state]= tmp[i][j];
            }
        }
    }

   hmm->s = s;
   ighmm_cmatrix_free(&tmp, mo->N);


    //params of states
    double *pi = ghmm_bayes_hmm_sample_initial(mo);
    double *emission;

    for(i = 0; i < mo->N; i++){
        ARRAY_MALLOC(emission, mo->dim );
        ghmm_hyperparameters_sample_emission_discrete(&(mo->params[i][0]), emission, mo->dim);
        hmm->s[i].b = emission;
        hmm->s[i].pi = pi[i];
    }

    return hmm;
STOP:
    return NULL;
#undef CUR_PROC
}

ghmm_cmodel* ghmm_bayes_hmm_sample_model(ghmm_bayes_hmm *mo){
#define CUR_PROC "ghmm_sample_model"
    //XXX Mixtures, 1 transitions, multivariate
    int i,j;
    ghmm_cmodel* hmm;
    hmm = ghmm_cmodel_calloc(mo->N, GHMM_kContinuousHMM, 1);
   
    ghmm_cstate *s =  hmm->s; 
    double **tmp;
    ARRAY_MALLOC(tmp, mo->N);
    hmm->M = 1;
    hmm->cos = 1;

    //get sample mats for transition
    for(i = 0; i < mo->N; i++){
        tmp[i] = ghmm_bayes_hmm_sample_transistions(mo, i);
    }

    //alloc 
    for(i=0; i<mo->N; i++){
        ARRAY_MALLOC(s[i].in_id, mo->N);
        ARRAY_MALLOC(s[i].out_id, mo->N);
        s[i].in_a = ighmm_cmatrix_alloc(mo->N, mo->N);
        s[i].out_a = ighmm_cmatrix_alloc(mo->N, mo->N);
    }

    //get in/out state transitions
    int out_state;
    int in_state;
    for(i=0; i<mo->N; i++){
        for(j=0; j<mo->N; j++){
            if(tmp[i][j] != 0){
                out_state = s[i].out_states++;
                in_state  = s[j].in_states++;
                s[i].out_id[out_state] = j;
                s[i].out_a[0][out_state] = tmp[i][j];
                s[j].in_id[in_state] = i;
                s[j].in_a[0][in_state]= tmp[i][j];
            }
        }
    }

    //realloc in out states
    for(i=0; i<mo->N; i++){
            ARRAY_REALLOC(s[i].in_id, s[i].in_states);
            ARRAY_REALLOC(s[i].out_id, s[i].out_states);
    }
    ighmm_cmatrix_free(&tmp, mo->N);
    hmm->s = s;


    //params of states
    double *pi = ghmm_bayes_hmm_sample_initial(mo);
    ghmm_c_emission *emission;
    double *weights;
    for(i = 0; i < mo->N; i++){
        ARRAY_MALLOC(weights, 1);
        ARRAY_MALLOC(emission,1);
        if(mo->dim > 1){
            ghmm_c_emission_alloc(emission, mo->dim);
        }
        weights[0] = 1.0;
        ghmm_hyperparameters_sample_emission(&(mo->params[i][0]), emission);//M
        emission->dimension = mo->dim;
        hmm->s[i].e = emission;
        hmm->s[i].pi = pi[i];
        hmm->s[i].c = weights;
        hmm->s[i].M = 1;
    }

    return hmm;
STOP:
    return NULL;
#undef CUR_PROC
}

void ghmm_hyperparameters_free(ghmm_hyperparameters *params){
#define CUR_PROC "free_hyperparameters"
    switch(params->type){
        case(normal):
            //m_free(params->emission.continuous);
        default:
            return;
    }
#undef CUR_PROC
}

void ghmm_hyperparameters_sample_emission_discrete(ghmm_hyperparameters *params, double *emission, int dim){
    ighmm_rand_dirichlet(0, dim, params->emission.discrete, emission);
}

void ghmm_hyperparameters_sample_emission(ghmm_hyperparameters *params, ghmm_c_emission *emission){
    switch(params->type){
        case(normal):
            {
                //printf("%x sample_emission: %f\n", params,  params->emission.continuous[1].beta);
                emission->type = normal;
                double precision = ighmm_rand_gamma(params->emission.continuous[1].alpha,
                        1/params->emission.continuous[1].beta, 0);
                emission->mean.val = ighmm_rand_normal(params->emission.continuous[0].mean.val, 
                        1/(precision * params->emission.continuous[0].variance.val), 0);
                emission->variance.val = 1/precision;
            }
            break;
        default:
            return;
    }
}


double* ghmm_bayes_hmm_sample_initial(ghmm_bayes_hmm *bayes){
#define CUR_PROC "ghmm_sample_initial"
    double *sample;
    ARRAY_CALLOC(sample, bayes->N);
    ighmm_rand_dirichlet(0, bayes->N, bayes->pi, sample);
    return sample;
STOP:
    return NULL;
#undef CUR_PROC
}

double* ghmm_bayes_hmm_sample_transistions(ghmm_bayes_hmm *bayes, int i){
#define CUR_PROC "ghmm_sample_transitions"
    double *prior = bayes->A[i];
    double *sample;
    ARRAY_MALLOC(sample, bayes->N);
    ighmm_rand_dirichlet(0, bayes->N, prior, sample);
    return sample;
STOP:
    return NULL;
#undef CUR_PROC
}

//M is mixture for continuous, alphabet size for discrete
int ghmm_bayes_hmm_alloc(ghmm_bayes_hmm *mo, int N, int *M){
#define CUR_PROC "ghmm_alloc_bayes_hmm"
    int i;
    mo->N = N;
    ARRAY_CALLOC(mo->M, N);
    for(i=0;i<N;i++)
        mo->M[i] = M[i];
    ARRAY_CALLOC(mo->params, N);
    for(i=0;i<N;i++){
        ARRAY_CALLOC(mo->params[i], M[i]);
    }
    mo->A = ighmm_cmatrix_alloc(N,N);
    ARRAY_CALLOC(mo->pi, N);

    return 0;
STOP:
    printf("error\n");
    return -1;
#undef CUR_PROC
}
