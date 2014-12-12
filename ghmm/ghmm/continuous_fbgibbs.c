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

#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <math.h>
#include <float.h>
#include "ghmm.h"
#include "mprintf.h"
#include "mes.h"
#include "matrix.h"
#include "randvar.h"
#include "ghmm_internals.h"
#include "rng.h"
#include "sfoba.h"
#include "continuous_fbgibbs.h"
#include "smodel.h"
#include "block_compression.h"

//data for each state for posterior
typedef struct sample_emission_data{
    int emitted;
    union{
        double val;
        double *vec;
    }mean;
    union{
        double val;
        double *vec;
        double **mat;
    }variance;
    double a;
    double b;
} sample_emission_data;

//data for model posterior
typedef struct ghmm_sample_data{
    double **transition;
    sample_emission_data **state_data;  //[state][mixture] rename to emission_data
}ghmm_sample_data;




int ghmm_alloc_sample_data(ghmm_bayes_hmm *mo, ghmm_sample_data *data){
#define CUR_PROC "ghmm_alloc_sample_data"
    int i;
    data->transition = ighmm_cmatrix_alloc(mo->N, mo->N);
    ARRAY_MALLOC(data->state_data, mo->N);
    for(i = 0; i < mo->N; i++){
        ARRAY_MALLOC(data->state_data[i], mo->M[i]);
        /*  XXX dim >1
        if(mo->dim > 1){  //init matrices for data for states
            for(i = 0; i < mo->M[i]; i++){ //each distribution in the mixture needs to alloc
                ghmm_alloc_emission_data(data->state_data[i][j], ghmm_bayes_hmm->params[i][j])
            }
        }
        */
    }
    return 0;
STOP:
    return -1;
#undef CUR_PROC
}
 
void ghmm_clear_emission_data(sample_emission_data *data){
    //XXX multivariate clear matrices instead
    data->emitted = 0;
    data->mean.val = 0;
    data->variance.val = 0;
    data->a = 0;
    data->b = 0;
}          

void ghmm_clear_sample_data(ghmm_sample_data * data, ghmm_bayes_hmm *bayes){
    int i, j;
    for(i = 0; i < bayes->N; i++){
        for(j = 0; j < bayes->M[i]; j++){
            ghmm_clear_emission_data(&data->state_data[i][j]);
        }
        for(j=0; j<bayes->N; j++){
            data->transition[i][j] = 0;
        }
    }
}
       
void ghmm_get_emission_data_first_pass(sample_emission_data *data, ghmm_density_t type,
        double *observation){
    //XXX add dist
    data->emitted++;
    switch(type){
        case(normal):
            data->mean.val += *observation;
            return;
        default://not supported
            return;
    }
}

void ghmm_get_emission_data_second_pass(sample_emission_data *data, ghmm_density_t type,
        double *observation){
    double tmp;
    //XXX add dist
    switch(type){
        case(normal):
            tmp = *observation - data->mean.val;
            data->variance.val += tmp*tmp;
            return;
        default://not supported
            return;
    }
}

void ghmm_get_sample_data(ghmm_sample_data *data, ghmm_bayes_hmm *bayes,int *Q, double *O, int T){
    int i;
    for(i=0; i<T-1; i++){
        data->transition[Q[i]][Q[i+1]]++;
    }
    for(i=0; i<T-1; i++){
        ghmm_get_emission_data_first_pass(&(data->state_data[Q[i]][0]),
                bayes->params[Q[i]][0].type, O+i);
    }

    for(i=0; i<bayes->N; i++){
        if(data->state_data[i][0].emitted>0)
            data->state_data[i][0].mean.val /= data->state_data[i][0].emitted;
    }

    for(i=0;i<T;i++){
        // a second pass through the observations since sample mean is needed
        ghmm_get_emission_data_second_pass(&data->state_data[Q[i]][0],
                bayes->params[Q[i]][0].type, &O[i]);
    }
}
   
/* using data colected in sample_emission_data sample from posterior distribution*/
void ghmm_update_emission(sample_emission_data *data, ghmm_hyperparameters *params,
        ghmm_c_emission *emission){
    //XXX add dist
    switch(params->type){
        case(normal):
            {
                double mean, var, a, b;
                double tmp;

                //var
                var = params->emission.continuous[0].variance.val + data->emitted;

                //mean
                mean = params->emission.continuous[0].variance.val * params->emission.continuous[0].mean.val;
                mean += data->emitted*data->mean.val;
                mean /= (params->emission.continuous[0].variance.val + data->emitted );

                //alpha
                a = params->emission.continuous[1].alpha + data->emitted/2;
                
                //beta
                tmp = data->mean.val - params->emission.continuous[0].mean.val;
                b = params->emission.continuous[1].beta;
                b += .5*data->variance.val;
                b += (data->emitted*params->emission.continuous[0].variance.val/2*tmp*tmp)/
                    (data->emitted+params->emission.continuous[0].variance.val);

                // sample from posterior hyperparameters
                if(params->emission.continuous[1].type == gamma_truncated){//XXX truncated_gamma 
                    tmp = ighmm_rand_truncated_gamma(params->emission.continuous[1].min, 
                            params->emission.continuous[1].max, a, 1/b, 0);
                }
                else{
                    tmp = ighmm_rand_gamma(a, 1/b, 0);
                }

                emission->variance.val = 1/tmp;
                if(params->emission.continuous[0].type == normal_right){

                    emission->mean.val = ighmm_rand_normal_right(params->emission.continuous[0].min,
                            mean, 1/(var*tmp),0);

                }
                else if (params->emission.continuous[0].type == normal_left) {
                    emission->mean.val = -ighmm_rand_normal_right(-params->emission.continuous[0].max, -mean,
                            1/(var*tmp),0);


                }
                else if (params->emission.continuous[0].type == truncated_normal) {
                     emission->mean.val = ighmm_rand_truncated_normal(params->emission.continuous[0].min,
                             params->emission.continuous[0].max, mean,
                            1/(var*tmp),0);

                }
                else{
                    emission->mean.val = ighmm_rand_normal(mean, 1/(var*tmp),0);
                }
                return;
            }
        default:
            return;
    }
}

void ghmm_update_model(ghmm_cmodel *mo, ghmm_bayes_hmm *bayes, ghmm_sample_data *data){
    int i, k;
    double tmp_n[mo->N];
    double tmp2_n[mo->N];

    //emission
    for(i=0; i<bayes->N; i++){
        //XXX mixture
        //for(k=0; k<bayes->M[k]; k++){
        k=0;
        ghmm_update_emission(&data->state_data[i][k], &bayes->params[i][k],&mo->s[i].e[k]);
        //}
    }

    //add prior for A, Pi
    for(i = 0; i < bayes->N; i++){
        tmp_n[i] = 0;
        for(k=0;k<bayes->M[i];k++){
            tmp_n[i] += data->state_data[i][k].emitted; 
        }
        for(k = 0; k < bayes->N; k++){
            data->transition[i][k] += bayes->A[i][k];
        }
    }

    //Pi
    ighmm_rand_dirichlet(0, mo->N, tmp_n, tmp2_n);
    for(k=0;k<mo->N;k++){
        mo->s[k].pi = tmp2_n[k];
    }

    //A
    for(i=0;i<mo->N;i++){
        ighmm_rand_dirichlet(0, mo->N, data->transition[i], tmp_n);
        for(k = 0; k < mo->N; k++){
            ghmm_cmodel_set_transition(mo, i, k, 0, tmp_n[k]);
        }
    }
}

//b is precomputed emissions used only for block compression
void ghmm_cmodel_fbgibbstep (ghmm_cmodel * mo, double *O, int len,int *Q, double** alpha, 
        double***pmats, double***b){
    int i,j,k;
    for(i = 0; i < len; i++){
        for(j = 0; j < mo->N; j++){
            alpha[i][j] = 0;
            for(k = 0; k < mo->N; k++){
               pmats[i][j][k] = 0;
            }
        }
    }

    double scale[len];
    double logP;
    ghmm_cmodel_forwardgibbs(mo, O, len, b, alpha, scale, &logP, pmats);
    sampleStatePath(mo->N, alpha[len-1], pmats, len, Q);
}



int** ghmm_cmodel_fbgibbs(ghmm_cmodel *mo, ghmm_bayes_hmm *bayes, ghmm_cseq* seq,
         int burnIn, int seed){
#define CUR_PROC "ghmm_cmodel_fbgibbs"
    GHMM_RNG_SET (RNG, seed);
    int max_seq = ghmm_cseq_max_len(seq);
    double **alpha = ighmm_cmatrix_alloc(max_seq,mo->N);
    double ***pmats = ighmm_cmatrix_3d_alloc(max_seq, mo->N, mo->N);
    int **Q; 
    ARRAY_CALLOC(Q, seq->seq_number);
    int seq_iter;
    for(seq_iter = 0; seq_iter < seq->seq_number; seq_iter++){
        ARRAY_CALLOC(Q[seq_iter], seq->seq_len[seq_iter]);
    }

    ghmm_sample_data data;
    ghmm_alloc_sample_data(bayes, &data);
    ghmm_clear_sample_data(&data, bayes);//XXX swap parameter
    for(; burnIn > 0; burnIn--){
        for(seq_iter = 0; seq_iter < seq->seq_number; seq_iter++){
            ghmm_cmodel_fbgibbstep(mo,seq->seq[seq_iter],seq->seq_len[seq_iter], Q[seq_iter],
                    alpha, pmats, NULL);
            ghmm_get_sample_data(&data, bayes, Q[seq_iter], seq->seq[seq_iter], 
                    seq->seq_len[seq_iter]); 
            ghmm_update_model(mo, bayes, &data);
            ghmm_clear_sample_data(&data, bayes);
        }
    }
    ighmm_cmatrix_free(&alpha, max_seq);
    ighmm_cmatrix_3d_free(&pmats, max_seq,mo->N);
    return Q;
STOP:
    return NULL; //XXX error handle
#undef CUR_PROC
}
//================== block compression ====================================

/*
based off of 
fast mcmc sampling for hidden markov models to determine copy number variations
md pavel mahmud, alexander schliep
*/

/* XXX mixture */
//[state][power]
//assumes normal dist
void precompute_blocks(ghmm_cmodel *mo, double ** mean, double **std,
        double **transition, int max_block_len){
    int i, j;
    for(i=0; i < mo->N; i++){
        mean[i][0] = 1;
        mean[i][1] = (mo->s+i)->e->mean.val * (mo->s+i)->e->mean.val;
        std[i][0] = sqrt( (mo->s+i)->e->variance.val);
        std[i][1] = std[i][0] * sqrt( 2 * M_PI  );
        //printf("std %e\n",  std[i][0]);
        transition[i][0] = 1;
        transition[i][1] = (mo->s+i)->out_a[0][i];

        //printf("transition %d %d %e\n", i, 0, transition[i][0]);
        //printf("transition %d %d %e\n", i, 1, transition[i][1]);
        for(j = 2; j < max_block_len; j++){
            std[i][j] = std[i][j-1] * std[i][1];
            transition[i][j] = transition[i][j-1]*transition[i][1];
            mean[i][j] = mean[i][j-1]+mean[i][1];
            //printf("transition %d %d %e\n", i, j, std[i][j]);
        }
    }
}

//assumes normal dist
//precomputes b for forward algo
void precompute_block_emission(ghmm_cmodel *mo, block_stats *stats, 
        int max_block_len, double ***b){
#define CUR_PROC "precalculate_block_emission"

    //precompute intermediate values
    double **mean2, **std, **transition;
    mean2 = ighmm_cmatrix_alloc(mo->N, max_block_len+1);
    std = ighmm_cmatrix_alloc(mo->N, max_block_len+1);
    transition = ighmm_cmatrix_alloc(mo->N, max_block_len+1);
    precompute_blocks(mo, mean2, std, transition, max_block_len+1);

    int t, i;
    double exponent;
    for(t = 0; t < stats->total; t++){
        for(i = 0; i < mo->N; i++){
            exponent = -1 * ( stats->moment2[t] - 2*stats->moment1[t] *
                    (mo->s+i)->e->mean.val + mean2[i][stats->length[t]] ) /
                    (2 * (mo->s+i)->e->variance.val);
                    
            b[t][i][1] = transition[i][stats->length[t]] * exp( exponent ) / 
                std[i][stats->length[t]];
        }
    }

    ighmm_cmatrix_free(&mean2, mo->N);
    ighmm_cmatrix_free(&std, mo->N);
    ighmm_cmatrix_free(&transition, mo->N);
STOP:
    //XXX ERROR
    return;
#undef CUR_PROC
}

void ghmm_get_sample_data_compressed(ghmm_sample_data *data, ghmm_bayes_hmm *bayes,
        int *Q, double *O, int T, block_stats *stats){
    int i,j,index;

    //T is the number of blocks
    for(i=0; i<T-1; i++){
        data->transition[Q[i]][Q[i+1]]++;
        data->transition[Q[i]][Q[i]] += stats->length[i];
    }
    data->transition[Q[T-1]][Q[T-1]] += stats->length[T-1];

    index = 0;
    //index is used to iterate through each observation in each block   
    for(i=0; i<T; i++){
        for(j = 0; j < stats->length[i]; j++){
            ghmm_get_emission_data_first_pass(&(data->state_data[Q[i]][0]),
                    bayes->params[Q[i]][0].type, O+index);
            index++;
        }
    }

    for(i=0; i<bayes->N; i++){
        if(data->state_data[i][0].emitted>0)
            data->state_data[i][0].mean.val /= data->state_data[i][0].emitted;
    }

    index = 0;
    for(i=0;i<T;i++){
        for(j = 0; j < stats->length[i]; j++){
            ghmm_get_emission_data_second_pass(&data->state_data[Q[i]][0],
                    bayes->params[Q[i]][0].type, O+index);
            index++;
        }
    }
}

int** ghmm_cmodel_cfbgibbs(ghmm_cmodel *mo, ghmm_bayes_hmm *bayes, ghmm_cseq* seq,
         int burnIn, int seed, double width, double delta, int max_len_permitted){
#define CUR_PROC "ghmm_cmodel_fbgibbs_compressed"
    GHMM_RNG_SET (RNG, seed);

    int seq_iter;
    int max_seq = ghmm_cseq_max_len(seq);

    //create blocks
    block_stats *stats[seq->seq_number];
    for(seq_iter = 0; seq_iter < seq->seq_number; seq_iter++){
        stats[seq_iter] = compress_observations(seq->seq[seq_iter], seq->seq_len[seq_iter],
                width*delta, delta);
        stats[seq_iter] = merge_observations(seq->seq[seq_iter], seq->seq_len[seq_iter],
                width, max_len_permitted, stats[seq_iter]);
    }

    //get max_block_len
    int max_block_len[seq_iter];
    int i;
    for(seq_iter = 0; seq_iter < seq->seq_number; seq_iter++){
        max_block_len[seq_iter] = stats[seq_iter]->length[0];
        for(i = 1; i < stats[seq_iter]->total; i++){
            if(max_block_len[seq_iter] < stats[seq_iter]->length[i])
                max_block_len[seq_iter] = stats[seq_iter]->length[i];
        }
    }

    double **alpha = ighmm_cmatrix_alloc(max_seq,mo->N);
    double ***pmats = ighmm_cmatrix_3d_alloc(max_seq, mo->N, mo->N);


    //used to precalculate the emision probabilities using compressed stats
    double ***b[seq->seq_number];
    for(seq_iter = 0; seq_iter < seq->seq_number; seq_iter++){
        b[seq_iter] = ighmm_cmatrix_3d_alloc(stats[seq_iter]->total, mo->N, mo->M+1);
    }

    int **Q; 
    ARRAY_CALLOC(Q, seq->seq_number);
    for(seq_iter = 0; seq_iter < seq->seq_number; seq_iter++){
        ARRAY_CALLOC(Q[seq_iter], seq->seq_len[seq_iter]);
    }

    ghmm_sample_data data;
    ghmm_alloc_sample_data(bayes, &data);
    ghmm_clear_sample_data(&data, bayes);//XXX swap parameter 
    for(; burnIn > 0; burnIn--){
        for(seq_iter = 0; seq_iter < seq->seq_number; seq_iter++){
            precompute_block_emission(mo, stats[seq_iter], max_block_len[seq_iter], b[seq_iter]);
            ghmm_cmodel_fbgibbstep(mo,seq->seq[seq_iter], stats[seq_iter]->total, Q[seq_iter], 
                    alpha, pmats, b[seq_iter]);
            ghmm_get_sample_data_compressed(&data, bayes, Q[seq_iter], seq->seq[seq_iter], 
                    stats[seq_iter]->total, stats[seq_iter]); 
            ghmm_update_model(mo, bayes, &data);
            ghmm_clear_sample_data(&data, bayes);
        }
    }

    ighmm_cmatrix_free(&alpha, max_seq);
    ighmm_cmatrix_3d_free(&pmats, max_seq,mo->N);
    for(seq_iter = 0; seq_iter < seq->seq_number; seq_iter++){
        ighmm_cmatrix_3d_free(b+seq_iter, stats[seq_iter]->total, mo->N);
        free_block_stats(stats+seq_iter);
    }

    return Q;
STOP:
    return NULL; //XXX error handle
#undef CUR_PROC
}

//======================== end block compression =======================================

/*XXX to reduce code duplication  maybe add pmats as param to original and ignore if NULL */ 

/*----------------------------------------------------------------------------*/
/*                     modified forwards to get the cdf pmats                 */
/*----------------------------------------------------------------------------*/

static double sfoba_stepforward_gibbs (ghmm_cstate * s, double *alpha_t, int osc,
                                 double b_omega, double* pmats, int N)
{
  int i, id, prv;
  double value = 0.0;
  prv = s->in_id[0];
  for (i = 0; i < s->in_states; i++) {
    id = s->in_id[i];
    pmats[id] = s->in_a[osc][i] * alpha_t[id];
    value += pmats[id];
    
    //fill in values of pmats previous id to current id
    for(; prv < id; prv++){
        pmats[prv+1] += pmats[prv];
    }
    prv = id;
  }
  for(prv+=1;prv<N;prv++){
      pmats[prv] += pmats[prv-1];
  }
  value *= b_omega;             /* b_omega outside the sum */
  return (value);
}                               /* sfoba_stepforward */


/*============================================================================*/
int ghmm_cmodel_forwardgibbs (ghmm_cmodel * smo, double *O, int T, double ***b,
                   double **alpha, double *scale, double *log_p, double***pmats)
{
# define CUR_PROC "ghmm_cmodel_forward"
  int res = -1;
  int i, t = 0, osc = 0;
  double c_t;
  int pos;

  /* T is length of sequence; divide by dimension to represent the number of time points */
  T /= smo->dim;
  /* calculate alpha and scale for t = 0 */
  if (b == NULL)
    sfoba_initforward(smo, alpha[0], O, scale, NULL);
  else
    sfoba_initforward(smo, alpha[0], O, scale, b[0]);
  if (scale[0] <= DBL_MIN) {
    /* means f(O[0], mue, u) << 0, first symbol very unlikely */
    /* GHMM_LOG(LCONVERTED, "scale[0] == .0!\n"); */
    goto STOP;
  }
  else {
    *log_p = -log (1 / scale[0]);

    if (smo->cos == 1) {
      osc = 0;
    }
    else {
      if (!smo->class_change->get_class) {
        printf ("ERROR: get_class not initialized\n");
        return (-1);
      }
      /* printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t); */
      osc = smo->class_change->get_class (smo, O, smo->class_change->k, t);
      if (osc >= smo->cos){
        printf("ERROR: get_class returned index %d but model has only %d classes !\n",osc,smo->cos);
        goto STOP;
      }

    }


    for (t = 1; t < T; t++) {
      scale[t] = 0.0;
      pos = t * smo->dim;
      /* b not calculated yet */
      if (b == NULL) {
        for (i = 0; i < smo->N; i++) {
          alpha[t][i] = sfoba_stepforward_gibbs(smo->s+i, alpha[t-1], osc,
                                          ghmm_cmodel_calc_b(smo->s+i, O+pos),
                                          pmats[t][i],smo->N);
          scale[t] += alpha[t][i];
        }
      }
      /* b precalculated */
      else {
        for (i = 0; i < smo->N; i++) {
          alpha[t][i] = sfoba_stepforward_gibbs (smo->s+i, alpha[t - 1], osc,
                                           b[t][i][smo->M], pmats[t][i], smo->N);
          scale[t] += alpha[t][i];
        }
      }
      if (scale[t] <= DBL_MIN) {        /* seq. can't be build */
        goto STOP;
        break;
      }
      c_t = 1 / scale[t];
      /* scale alpha */
      for (i = 0; i < smo->N; i++)
        alpha[t][i] *= c_t;
      /* summation of log(c[t]) for calculation of log( P(O|lambda) ) */
      *log_p -= log (c_t);

      if (smo->cos == 1) {
        osc = 0;
      }
      else {
        if (!smo->class_change->get_class) {
          printf ("ERROR: get_class not initialized\n");
          return (-1);
        }
        /* printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t); */
        osc = smo->class_change->get_class (smo, O, smo->class_change->k, t);
        if (osc >= smo->cos){
          printf("ERROR: get_class returned index %d but model has only %d classes !\n",osc,smo->cos);
          goto STOP;
        }		
      }
    }
  }
  /* log_p should not be smaller than value used for seqs. that 
     can't be build ???
     if (*log_p < (double)PENALTY_LOGP)
     *log_p = (double)PENALTY_LOGP;
   */
  return 0;
STOP:
  *log_p = (double) -DBL_MAX;
  return (res);
#undef CUR_PROC
}                               /* ghmm_cmodel_forward */
//====================================================================================
//==================================end forwards======================================
//====================================================================================
