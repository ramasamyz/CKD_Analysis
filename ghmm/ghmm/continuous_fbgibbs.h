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

#ifndef GHMM_CONTINUOUS_FBGIBBS_H
#define GHMM_CONTINUOUS_FBGIBBS_H

#include "ghmm/smodel.h"
#include "bayesian_hmm.h"
/* uses fbgibbs to get a sampled model from the posterior 
 * @param bayes: the prior model
 * @param mo: sampled model from priors after function called is sampled from posterior
 * @param seq: observation sequence 
 * @param burnIn: number of times to iterate
 * @param seed: 0 for defualt seed non zero for new seed
 * @return sampled state sequence
 */
int** ghmm_cmodel_fbgibbs(ghmm_cmodel *mo, ghmm_bayes_hmm *bayes,  ghmm_cseq* seq,
         int burnIn, int seed);


/* based off of 
 * fast mcmc sampling for hidden markov models to determine copy number variations
 * md pavel mahmud, alexander schliep
 *
 * XXX NORMAL ONLY
 * uses a compressed sequence to approximate fbgibbs to get a sampled model from the posterior 
 * @param bayes: the prior model
 * @param mo: sampled model from priors after function called is sampled from posterior
 * @param seq: observation sequence 
 * @param burnIn: number of times to iterate
 * @param seed: 0 for defualt seed non zero for new seed
 * @param width: parameter that controls size of block
 * @param delta: shrink factor for the blocks
 * @param max_len_permitted: maximum length of block
 * @return sampled state sequence
 */
int** ghmm_cmodel_cfbgibbs( ghmm_cmodel *mo, ghmm_bayes_hmm *bayes, ghmm_cseq* seq,
         int burnIn, int seed, double width, double delta, int max_len_permitted);
#endif
