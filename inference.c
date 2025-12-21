#include "mongrail.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*============================================================================
 * CACHE INITIALIZATION AND MANAGEMENT
 *
 * Pre-computes values that are reused across likelihood calculations:
 * - Q(z) values for each ancestry configuration z (per chromosome)
 * - log P(hap|popA) and log P(hap|popB) for all haplotypes
 * - U(hap) values are cached on first computation
 *============================================================================*/

struct likelihood_cache* cache_init(int noChr, int* no_loci, double recomb_rate,
				    unsigned long** marker_positions,
				    int** popA_hap_counts, int** popB_hap_counts,
				    unsigned int** haplist, int* no_haps,
				    int noSamplesPopA, int noSamplesPopB)
{
  struct likelihood_cache* cache = malloc(sizeof(struct likelihood_cache));
  cache->initialized = 1;
  cache->noChr = noChr;
  cache->no_loci = no_loci;
  cache->log2 = log(2.0);

  /* Allocate and pre-compute log Q(z) values */
  cache->logQ_cache = malloc(noChr * sizeof(double*));
  for (int c = 0; c < noChr; c++) {
    int num_z = 1 << no_loci[c];
    cache->logQ_cache[c] = malloc(num_z * sizeof(double));
    for (int z = 0; z < num_z; z++) {
      double Q_z = prQ(z, no_loci[c], recomb_rate, marker_positions[c]);
      cache->logQ_cache[c][z] = (Q_z > 0) ? log(Q_z) : -1e300;
    }
  }

  /* Allocate U cache (initialized to -1 = not computed) */
  cache->U_cache = malloc(noChr * sizeof(double*));
  for (int c = 0; c < noChr; c++) {
    cache->U_cache[c] = malloc(MAXHAPS * sizeof(double));
    for (int h = 0; h < MAXHAPS; h++) {
      cache->U_cache[c][h] = -1.0;  /* sentinel for "not computed" */
    }
  }

  /* Pre-compute log P(hap|popA) and log P(hap|popB) for all haplotypes in haplist */
  cache->logP_A = malloc(noChr * sizeof(double*));
  cache->logP_B = malloc(noChr * sizeof(double*));
  for (int c = 0; c < noChr; c++) {
    cache->logP_A[c] = malloc(MAXHAPS * sizeof(double));
    cache->logP_B[c] = malloc(MAXHAPS * sizeof(double));
    /* Initialize to sentinel */
    for (int h = 0; h < MAXHAPS; h++) {
      cache->logP_A[c][h] = -1e300;
      cache->logP_B[c][h] = -1e300;
    }
    /* Compute for haplotypes in haplist */
    for (int i = 0; i < no_haps[c]; i++) {
      unsigned int hap = haplist[c][i];
      cache->logP_A[c][hap] = log_prob_single_hap(hap, popA_hap_counts, haplist,
						  no_haps, noSamplesPopA, c);
      cache->logP_B[c][hap] = log_prob_single_hap(hap, popB_hap_counts, haplist,
						  no_haps, noSamplesPopB, c);
    }
  }

  return cache;
}

void cache_free(struct likelihood_cache* cache)
{
  if (!cache) return;
  for (int c = 0; c < cache->noChr; c++) {
    free(cache->logQ_cache[c]);
    free(cache->U_cache[c]);
    free(cache->logP_A[c]);
    free(cache->logP_B[c]);
  }
  free(cache->logQ_cache);
  free(cache->U_cache);
  free(cache->logP_A);
  free(cache->logP_B);
  free(cache);
}

/* Cached version of log_U_given_z - uses precomputed logP values */
static double log_U_given_z_cached(struct likelihood_cache* cache,
				   unsigned int hap, unsigned int z, int noloci,
				   int** popA_hap_counts, int** popB_hap_counts,
				   unsigned int** haplist, int* no_haps,
				   int noSamplesPopA, int noSamplesPopB, int chrom)
{
  unsigned int all_A = 0;
  unsigned int all_B = (1U << noloci) - 1;

  if (z == all_A) {
    /* All loci from pop A - use cached value if available */
    if (cache->logP_A[chrom][hap] > -1e299)
      return cache->logP_A[chrom][hap];
    return log_prob_single_hap(hap, popA_hap_counts, haplist, no_haps, noSamplesPopA, chrom);
  } else if (z == all_B) {
    /* All loci from pop B - use cached value if available */
    if (cache->logP_B[chrom][hap] > -1e299)
      return cache->logP_B[chrom][hap];
    return log_prob_single_hap(hap, popB_hap_counts, haplist, no_haps, noSamplesPopB, chrom);
  } else {
    /* Recombinant - split into sub-haplotypes */
    unsigned int ancvec_B = z;
    unsigned int ancvec_A = get_anc_complement(z, noloci);
    double log_prob_A = log_prob_subhap(hap, ancvec_A, popA_hap_counts, haplist,
					no_haps, noSamplesPopA, chrom);
    double log_prob_B = log_prob_subhap(hap, ancvec_B, popB_hap_counts, haplist,
					no_haps, noSamplesPopB, chrom);
    return log_prob_A + log_prob_B;
  }
}

/* Compute U(hap) with caching - returns cached value or computes and caches */
double compute_U_cached(struct likelihood_cache* cache, unsigned int hap, int chrom,
			int** popA_hap_counts, int** popB_hap_counts,
			unsigned int** haplist, int* no_haps,
			int noSamplesPopA, int noSamplesPopB)
{
  /* Check cache first */
  if (cache->U_cache[chrom][hap] >= 0) {
    return cache->U_cache[chrom][hap];
  }

  int noloci = cache->no_loci[chrom];
  int num_z = 1 << noloci;

  /* Use static buffer to avoid malloc - thread-unsafe but fast */
  static double log_terms[MAXANCESTRY];

  /* Compute log(U(hap|z) * Q(z)) for each z using cached log Q values */
  for (int z = 0; z < num_z; z++) {
    double log_U_z = log_U_given_z_cached(cache, hap, z, noloci,
					  popA_hap_counts, popB_hap_counts,
					  haplist, no_haps,
					  noSamplesPopA, noSamplesPopB, chrom);
    log_terms[z] = cache->logQ_cache[chrom][z] + log_U_z;
  }

  /* Log-sum-exp for numerical stability */
  double max_log = log_terms[0];
  for (int i = 1; i < num_z; i++) {
    if (log_terms[i] > max_log) max_log = log_terms[i];
  }

  double sum = 0.0;
  for (int i = 0; i < num_z; i++) {
    sum += exp(log_terms[i] - max_log);
  }

  double result = exp(max_log + log(sum));

  /* Cache the result */
  cache->U_cache[chrom][hap] = result;

  return result;
}

/*============================================================================
 * OPTIMIZED LIKELIHOOD FUNCTIONS FOR MODELS B, E, F
 * These use the cache for faster computation
 *============================================================================*/

double lik_b_e_cached(struct likelihood_cache* cache, int indivIndex,
		      struct indiv** hybrid_indiv,
		      int** popB_hap_counts, int** popA_hap_counts,
		      unsigned int** haplist, int* no_haps,
		      int noSamplesPopB, int noSamplesPopA, int noChr,
		      enum modelType model)
{
  unsigned int hap1, hap2;
  double logL = 0.0;
  static double probV[MAXHAPS];  /* static to avoid malloc */

  for (int i = 0; i < noChr; i++) {
    double probL = 0.0;

    for (int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k += 2) {
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      /* Get cached U values */
      double U1 = compute_U_cached(cache, hap1, i, popA_hap_counts, popB_hap_counts,
				   haplist, no_haps, noSamplesPopA, noSamplesPopB);
      double U2 = compute_U_cached(cache, hap2, i, popA_hap_counts, popB_hap_counts,
				   haplist, no_haps, noSamplesPopA, noSamplesPopB);

      /* Get cached log P values for pure haplotype */
      double p1_pure, p2_pure;
      if (model == MODEL_B) {
	p1_pure = exp(cache->logP_A[i][hap1]);
	p2_pure = exp(cache->logP_A[i][hap2]);
      } else {
	p1_pure = exp(cache->logP_B[i][hap1]);
	p2_pure = exp(cache->logP_B[i][hap2]);
      }

      if (identity2_hap(hap1, hap2)) {
	probV[k/2] = p1_pure * U2 + p2_pure * U1;
      } else {
	probV[k/2] = p1_pure * U2;
      }
    }

    int numDiplotypes = hybrid_indiv[i][indivIndex].numHaps / 2;
    double norm_factor = max_element(probV, numDiplotypes);
    if (norm_factor > 0) {
      double log_norm = log(norm_factor);
      for (int a = 0; a < numDiplotypes; a++) {
	probV[a] = exp(log(probV[a]) - cache->log2 - log_norm);
      }
      probL = log(kahanSum(probV, numDiplotypes)) + cache->log2 + log_norm;
    } else {
      probL = -1e300;
    }
    logL += probL;
  }

  return logL;
}

double lik_f_cached(struct likelihood_cache* cache, int indivIndex,
		    struct indiv** hybrid_indiv,
		    int** popB_hap_counts, int** popA_hap_counts,
		    unsigned int** haplist, int* no_haps,
		    int noSamplesPopB, int noSamplesPopA, int noChr)
{
  unsigned int hap1, hap2;
  double logL = 0.0;
  static double probV[MAXHAPS];

  for (int i = 0; i < noChr; i++) {
    double probL = 0.0;

    for (int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k += 2) {
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      double U1 = compute_U_cached(cache, hap1, i, popA_hap_counts, popB_hap_counts,
				   haplist, no_haps, noSamplesPopA, noSamplesPopB);
      double U2 = compute_U_cached(cache, hap2, i, popA_hap_counts, popB_hap_counts,
				   haplist, no_haps, noSamplesPopA, noSamplesPopB);

      if (identity2_hap(hap1, hap2)) {
	probV[k/2] = 2.0 * U1 * U2;
      } else {
	probV[k/2] = U1 * U2;
      }
    }

    int numDiplotypes = hybrid_indiv[i][indivIndex].numHaps / 2;
    double norm_factor = max_element(probV, numDiplotypes);
    if (norm_factor > 0) {
      double log_norm = log(norm_factor);
      for (int a = 0; a < numDiplotypes; a++) {
	probV[a] = exp(log(probV[a]) - cache->log2 - log_norm);
      }
      probL = log(kahanSum(probV, numDiplotypes)) + cache->log2 + log_norm;
    } else {
      probL = -1e300;
    }
    logL += probL;
  }

  return logL;
}

/* likelihood for models a and d */
double lik_a_d(int indivIndex, struct indiv** hybrid_indiv, int** popY_hap_counts,
	       unsigned int** haplist, int* no_haps, int noSamplesPopA,
	       int noSamplesPopB, int noChr, enum modelType model)
{
  unsigned int hap1, hap2;
  int phi = 0;
  double norm_factor;
  double logL = 0.0;
  double probL = 0.0;
  double *probV = (double *) malloc(MAXHAPS*sizeof(double));
  double t1, t3, t4, t5;

  if(model == MODEL_A)
    {
      t1 = T1B;
      t4 = T4B;
    }
  else
    {
      t1 = T1A;
      t4 = T4A;
    }
  assert(t1!=0.0); /* not MODEL_A or MODEL_D */
  
  /* loop over chromosomes */  
  for(int i=0; i<noChr; i++)
    {
      probL = 0.0;
      for(int k=0; k<hybrid_indiv[i][indivIndex].numHaps; k=k+2)
	{
	  hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
	  hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];
	  /* likelihood calculations */
	  t3 = identity2_hap(hap1, hap2)*L2; 
	  t5=0.0;
	  for(int j=0; j<no_haps[i]; j++)
	    {
	      phi = identity1_hap(haplist[i][j], hap1) + identity1_hap(haplist[i][j], hap2);
	      t5 += gammln(phi + popY_hap_counts[i][haplist[i][j]] + 1.0/no_haps[i]) -
		gammln(popY_hap_counts[i][haplist[i][j]]+1.0/no_haps[i]);
	    }
	  probV[k/2] = exp(t1 + t3 - t4 + t5);
	}
      /* normalize probabilities to avoid underflow */
      norm_factor = max_element(probV,hybrid_indiv[i][indivIndex].numHaps/2);
      for(int a=0; a<hybrid_indiv[i][indivIndex].numHaps/2; a++)
	probV[a]=exp(log(probV[a])-log(2.0)-log(norm_factor));
      probL = log(kahanSum(probV, hybrid_indiv[i][indivIndex].numHaps/2))+log(2.0)+log(norm_factor);
      logL += probL;
    }
  free(probV);
  return logL;
}

/* likelihood for model c */
double lik_c(int indivIndex, struct indiv** hybrid_indiv, int** popB_hap_counts, int** popA_hap_counts,
	     unsigned int** haplist, int* no_haps, int noSamplesPopB,  int noSamplesPopA, int noChr)
{
  unsigned int hap1, hap2;
  int phi1 = 0;
  int phi2 = 0;
  double logL = 0.0;
  double probL = 0.0;
  double norm_factor;
  double *probV = (double *) malloc(MAXHAPS*sizeof(double));
  double t2A, t2B;
  double t4h1A, t4h2A, t4h1B, t4h2B;
  double p1A, p1B, p2A, p2B;

  for(int i=0; i<noChr; i++)
    {
      probL = 0.0;
      for(int k=0; k<hybrid_indiv[i][indivIndex].numHaps; k=k+2) 
	{
	  hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
	  hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];
	  /* likelihood calculations */
	  t2A = 0.0;
	  t2B = 0.0;
	  t4h1A = 0.0;
	  t4h2A = 0.0;
	  t4h1B = 0.0;
	  t4h2B = 0.0;
	  for(int j=0; j<no_haps[i]; j++)
	    {
	      t2A += gammln(popA_hap_counts[i][haplist[i][j]]+1.0/no_haps[i]);
	      t2B += gammln(popB_hap_counts[i][haplist[i][j]]+1.0/no_haps[i]);
	    }
	  for(int j=0; j<no_haps[i]; j++)
	    {
	      phi1 = identity1_hap(haplist[i][j], hap1);
	      phi2 = identity1_hap(haplist[i][j], hap2);
	      t4h1A += gammln(phi1 + popA_hap_counts[i][haplist[i][j]] + 1.0/no_haps[i]);
	      t4h2A += gammln(phi2 + popA_hap_counts[i][haplist[i][j]] + 1.0/no_haps[i]);
	      t4h1B += gammln(phi1 + popB_hap_counts[i][haplist[i][j]] + 1.0/no_haps[i]);
	      t4h2B += gammln(phi2 + popB_hap_counts[i][haplist[i][j]] + 1.0/no_haps[i]);
	    }
	  p1A = exp(T1A - t2A - TC4A + t4h1A);
	  p1B = exp(T1B - t2B - TC4B + t4h1B);
	  p2A = exp(T1A - t2A - TC4A + t4h2A);
	  p2B = exp(T1B - t2B - TC4B + t4h2B);
	  probV[k/2] = identity2_hap(hap1, hap2)*(p1A*p2B+p2A*p1B)+(1 - identity2_hap(hap1, hap2))*(p1A*p2B);
	}
      norm_factor = max_element(probV,hybrid_indiv[i][indivIndex].numHaps/2);
      /* normalize probabilities to avoid underflow */
      for(int a=0; a<hybrid_indiv[i][indivIndex].numHaps/2; a++)
	probV[a]=exp(log(probV[a])-log(2.0)-log(norm_factor));
      probL = log(kahanSum(probV, hybrid_indiv[i][indivIndex].numHaps/2))+log(2.0)+log(norm_factor);
      logL += probL;
    }
  free(probV);
  return logL;
}

/*============================================================================
 * HELPER FUNCTIONS FOR MODELS B AND E (BACKCROSS MODELS)
 *
 * These functions implement the recombinant haplotype probability U(x^c_i)
 * as described in the Mongrail 2.0 paper (Appendix A).
 *
 * U(x^c_i) = Σ_z U(x^c_i|z) · Q(z|d,r)
 *
 * Where:
 *   z = ancestry vector (which loci come from pop A vs B)
 *   Q(z|d,r) = probability of ancestry configuration given recombination
 *   U(x^c_i|z) = probability of haplotype given ancestry state
 *============================================================================*/

/* Get distinct sub-haplotypes and their marginal counts from a population.
 * ancvec specifies which loci to consider (1 = include, 0 = exclude).
 * Returns the number of distinct sub-haplotypes found.
 */
int get_distinct_subhaps(unsigned int ancvec, int** pop_hap_counts,
			 unsigned int** haplist, int* no_haps, int chrom,
			 unsigned int* distinct_subhaps, int* marginal_counts)
{
  int num_distinct = 0;
  for(int i = 0; i < no_haps[chrom]; i++) {
    unsigned int subhap = haplist[chrom][i] & ancvec;
    int found = 0;
    for(int j = 0; j < num_distinct; j++) {
      if(distinct_subhaps[j] == subhap) {
	marginal_counts[j] += pop_hap_counts[chrom][haplist[chrom][i]];
	found = 1;
	break;
      }
    }
    if(!found) {
      distinct_subhaps[num_distinct] = subhap;
      marginal_counts[num_distinct] = pop_hap_counts[chrom][haplist[chrom][i]];
      num_distinct++;
    }
  }
  return num_distinct;
}

/* Compute log probability of a single haplotype from population k.
 * This is log Pr(φ^c(x_i)|n_k) from the paper.
 * Uses posterior predictive distribution with Dirichlet prior.
 */
double log_prob_single_hap(unsigned int hap, int** pop_hap_counts,
			   unsigned int** haplist, int* no_haps,
			   int noSamplesPop, int chrom)
{
  double theta = 1.0 + 2.0 * noSamplesPop;
  double H = (double)no_haps[chrom];
  double t1 = gammln(theta);
  double t2 = 0.0;
  double t3 = gammln(1.0 + theta);
  double t4 = 0.0;

  for(int j = 0; j < no_haps[chrom]; j++) {
    t2 += gammln(pop_hap_counts[chrom][haplist[chrom][j]] + 1.0/H);
    int phi = identity1_hap(haplist[chrom][j], hap);
    t4 += gammln(phi + pop_hap_counts[chrom][haplist[chrom][j]] + 1.0/H);
  }

  return t1 - t2 - t3 + t4;
}

/* Compute log probability of a sub-haplotype given ancestry mask.
 * ancvec specifies which loci are included in this sub-haplotype.
 * Uses marginal counts from the reference population.
 */
double log_prob_subhap(unsigned int hap, unsigned int ancvec, int** pop_hap_counts,
		       unsigned int** haplist, int* no_haps, int noSamplesPop, int chrom)
{
  unsigned int distinct_subhaps[MAXHAPS];
  int marginal_counts[MAXHAPS];

  /* Initialize marginal_counts to zero */
  for(int i = 0; i < MAXHAPS; i++) {
    marginal_counts[i] = 0;
  }

  int m = get_distinct_subhaps(ancvec, pop_hap_counts, haplist, no_haps, chrom,
			       distinct_subhaps, marginal_counts);

  if(m == 0) return -1e300;  /* No data - return very negative log prob */

  unsigned int target_subhap = hap & ancvec;
  double theta_prime = 1.0 + 2.0 * noSamplesPop;

  double t1 = gammln(theta_prime);
  double t2 = 0.0;
  double t3 = gammln(1.0 + theta_prime);
  double t4 = 0.0;

  for(int i = 0; i < m; i++) {
    t2 += gammln(marginal_counts[i] + 1.0/m);
    int indicator = (distinct_subhaps[i] == target_subhap) ? 1 : 0;
    t4 += gammln(indicator + marginal_counts[i] + 1.0/m);
  }

  return t1 - t2 - t3 + t4;
}

/* Compute log U(x^c_i|z) - probability of haplotype given ancestry vector z.
 * If z is all 0s (all from A) or all 1s (all from B), use full haplotype prob.
 * Otherwise, split into sub-haplotypes and multiply their probabilities.
 */
double log_U_given_z(unsigned int hap, unsigned int z, int noloci,
		     int** popA_hap_counts, int** popB_hap_counts,
		     unsigned int** haplist, int* no_haps,
		     int noSamplesPopA, int noSamplesPopB, int chrom)
{
  unsigned int all_A = 0;  /* z = 0 means all from A */
  unsigned int all_B = (1U << noloci) - 1;  /* z = 2^L - 1 means all from B */

  if(z == all_A) {
    /* All loci from pop A */
    return log_prob_single_hap(hap, popA_hap_counts, haplist, no_haps, noSamplesPopA, chrom);
  } else if(z == all_B) {
    /* All loci from pop B */
    return log_prob_single_hap(hap, popB_hap_counts, haplist, no_haps, noSamplesPopB, chrom);
  } else {
    /* Recombinant - split into sub-haplotypes from A and B */
    unsigned int ancvec_B = z;  /* Where z=1 (from B) */
    unsigned int ancvec_A = get_anc_complement(z, noloci);  /* Where z=0 (from A) */

    double log_prob_A = log_prob_subhap(hap, ancvec_A, popA_hap_counts, haplist, no_haps, noSamplesPopA, chrom);
    double log_prob_B = log_prob_subhap(hap, ancvec_B, popB_hap_counts, haplist, no_haps, noSamplesPopB, chrom);

    return log_prob_A + log_prob_B;
  }
}

/* Compute U(x^c_i) - the recombinant haplotype probability.
 * Sums over all possible ancestry configurations z, weighted by Q(z|d,r).
 * Returns the probability (not log probability) for numerical stability.
 */
double compute_U(unsigned int hap, int noloci, double recomb_rate,
		 unsigned long* positions, int** popA_hap_counts, int** popB_hap_counts,
		 unsigned int** haplist, int* no_haps,
		 int noSamplesPopA, int noSamplesPopB, int chrom)
{
  int num_z = 1 << noloci;  /* 2^L ancestry configurations */
  double* log_terms = (double*)malloc(num_z * sizeof(double));

  /* Compute log(U(hap|z) * Q(z)) for each z */
  for(unsigned int z = 0; z < (unsigned int)num_z; z++) {
    double Q_z = prQ(z, noloci, recomb_rate, positions);
    double log_Q = (Q_z > 0) ? log(Q_z) : -1e300;
    double log_U_z = log_U_given_z(hap, z, noloci, popA_hap_counts, popB_hap_counts,
				   haplist, no_haps, noSamplesPopA, noSamplesPopB, chrom);
    log_terms[z] = log_Q + log_U_z;
  }

  /* Sum in log space using log-sum-exp trick for numerical stability */
  double max_log = log_terms[0];
  for(int i = 1; i < num_z; i++) {
    if(log_terms[i] > max_log) max_log = log_terms[i];
  }

  double sum = 0.0;
  for(int i = 0; i < num_z; i++) {
    sum += exp(log_terms[i] - max_log);
  }

  free(log_terms);

  /* Return U as probability (not log) */
  return exp(max_log + log(sum));
}

/*============================================================================
 * LIKELIHOOD FOR MODELS B AND E (BACKCROSS MODELS)
 *
 * Model B: Backcross with Population A
 *   - One haplotype comes directly from population A (pure)
 *   - Other haplotype is a recombinant from F1 parent (mosaic of A and B)
 *
 * Model E: Backcross with Population B
 *   - One haplotype comes directly from population B (pure)
 *   - Other haplotype is a recombinant from F1 parent (mosaic of A and B)
 *
 * Log-likelihood formula from paper:
 *   log P(x_i|g=b,n) = I(x_i)*log[Pr(φ^M|n_A)*U(x^P) + Pr(φ^P|n_A)*U(x^M)]
 *                    + [1-I(x_i)]*log[Pr(φ^M|n_A)*U(x^P)]
 *============================================================================*/
double lik_b_e(int indivIndex, struct indiv** hybrid_indiv,
	       int** popB_hap_counts, int** popA_hap_counts,
	       unsigned int** haplist, int* no_haps,
	       int noSamplesPopB, int noSamplesPopA, int noChr,
	       int* no_loci, unsigned long** marker_positions, double recomb_rate,
	       enum modelType model)
{
  unsigned int hap1, hap2;
  double logL = 0.0;
  double norm_factor;
  double *probV = (double *)malloc(MAXHAPS * sizeof(double));
  double probL = 0.0;

  /* Loop over chromosomes */
  for(int i = 0; i < noChr; i++) {
    probL = 0.0;

    /* Loop over compatible haplotype pairs (diplotypes) */
    for(int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k = k + 2) {
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      /* Compute U(hap1) and U(hap2) - recombinant haplotype probabilities */
      double U1 = compute_U(hap1, no_loci[i], recomb_rate, marker_positions[i],
			    popA_hap_counts, popB_hap_counts, haplist, no_haps,
			    noSamplesPopA, noSamplesPopB, i);
      double U2 = compute_U(hap2, no_loci[i], recomb_rate, marker_positions[i],
			    popA_hap_counts, popB_hap_counts, haplist, no_haps,
			    noSamplesPopA, noSamplesPopB, i);

      /* Compute Pr(hap|pop) for the "pure" haplotype */
      double p1_pure, p2_pure;
      if(model == MODEL_B) {
	/* Model B: pure haplotype from pop A */
	p1_pure = exp(log_prob_single_hap(hap1, popA_hap_counts, haplist, no_haps, noSamplesPopA, i));
	p2_pure = exp(log_prob_single_hap(hap2, popA_hap_counts, haplist, no_haps, noSamplesPopA, i));
      } else {
	/* Model E: pure haplotype from pop B */
	p1_pure = exp(log_prob_single_hap(hap1, popB_hap_counts, haplist, no_haps, noSamplesPopB, i));
	p2_pure = exp(log_prob_single_hap(hap2, popB_hap_counts, haplist, no_haps, noSamplesPopB, i));
      }

      /* Combine: one haplotype pure, other recombinant */
      if(identity2_hap(hap1, hap2)) {
	/* Haplotypes are different - sum over both assignments */
	probV[k/2] = p1_pure * U2 + p2_pure * U1;
      } else {
	/* Haplotypes are identical */
	probV[k/2] = p1_pure * U2;
      }
    }

    /* Normalize probabilities to avoid underflow */
    norm_factor = max_element(probV, hybrid_indiv[i][indivIndex].numHaps/2);
    if(norm_factor > 0) {
      for(int a = 0; a < hybrid_indiv[i][indivIndex].numHaps/2; a++) {
	probV[a] = exp(log(probV[a]) - log(2.0) - log(norm_factor));
      }
      probL = log(kahanSum(probV, hybrid_indiv[i][indivIndex].numHaps/2)) + log(2.0) + log(norm_factor);
    } else {
      probL = -1e300;  /* Very negative log probability if all probs are zero */
    }
    logL += probL;
  }

  free(probV);
  return logL;
}

/*============================================================================
 * LIKELIHOOD FOR MODEL F (F2 HYBRID)
 *
 * Model F: Second-generation hybrid (F2)
 *   - Both parents were F1 hybrids
 *   - Both haplotypes are recombinant (mosaic of A and B ancestry)
 *
 * Log-likelihood formula:
 *   log P(x_i|g=f,n) = I(x_i)*log[U(x^M)*U(x^P) + U(x^P)*U(x^M)]
 *                    + [1-I(x_i)]*log[U(x^M)*U(x^P)]
 *
 * Where U(x) is the recombinant haplotype probability from compute_U().
 *============================================================================*/
double lik_f(int indivIndex, struct indiv** hybrid_indiv,
	     int** popB_hap_counts, int** popA_hap_counts,
	     unsigned int** haplist, int* no_haps,
	     int noSamplesPopB, int noSamplesPopA, int noChr,
	     int* no_loci, unsigned long** marker_positions, double recomb_rate)
{
  unsigned int hap1, hap2;
  double logL = 0.0;
  double norm_factor;
  double *probV = (double *)malloc(MAXHAPS * sizeof(double));
  double probL = 0.0;

  /* Loop over chromosomes */
  for(int i = 0; i < noChr; i++) {
    probL = 0.0;

    /* Loop over compatible haplotype pairs (diplotypes) */
    for(int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k = k + 2) {
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      /* Compute U(hap1) and U(hap2) - recombinant haplotype probabilities */
      double U1 = compute_U(hap1, no_loci[i], recomb_rate, marker_positions[i],
			    popA_hap_counts, popB_hap_counts, haplist, no_haps,
			    noSamplesPopA, noSamplesPopB, i);
      double U2 = compute_U(hap2, no_loci[i], recomb_rate, marker_positions[i],
			    popA_hap_counts, popB_hap_counts, haplist, no_haps,
			    noSamplesPopA, noSamplesPopB, i);

      /* Both haplotypes are recombinant */
      if(identity2_hap(hap1, hap2)) {
	/* Haplotypes are different - sum over both assignments */
	probV[k/2] = 2.0 * U1 * U2;
      } else {
	/* Haplotypes are identical */
	probV[k/2] = U1 * U2;
      }
    }

    /* Normalize probabilities to avoid underflow */
    norm_factor = max_element(probV, hybrid_indiv[i][indivIndex].numHaps/2);
    if(norm_factor > 0) {
      for(int a = 0; a < hybrid_indiv[i][indivIndex].numHaps/2; a++) {
	probV[a] = exp(log(probV[a]) - log(2.0) - log(norm_factor));
      }
      probL = log(kahanSum(probV, hybrid_indiv[i][indivIndex].numHaps/2)) + log(2.0) + log(norm_factor);
    } else {
      probL = -1e300;  /* Very negative log probability if all probs are zero */
    }
    logL += probL;
  }

  free(probV);
  return logL;
}

/* Kahan summation avoids numerical error */
double kahanSum(double arr[], int n) {
    double sum = 0.0;
    double c = 0.0; // A running compensation for lost low-order bits

    for (int i = 0; i < n; i++) {
        double y = arr[i] - c;
        double t = sum + y;
        c = (t - sum) - y; // Calculate the lost low-order bits
        sum = t;
    }
    return sum;
}

/* find maximum element in an array */
double max_element(double arr[], int n)
{
  double currmax=0.0;
  for(int i=0; i<n; i++)
    currmax = fmax(currmax,arr[i]);  /* use fmax for double, not fmaxf (float) */
  return(currmax);
}

/* add haplotype to haplotype list if not already present */
void add_hap(unsigned int hap, unsigned int** haplist, int* no_haps, int chrom)
{
  int hap_found = 0;
  if(no_haps[chrom] == 0)
    {
      haplist[chrom][0] = hap;
      no_haps[chrom] = 1;
    }
  else
    {
      for(int i=0; i < no_haps[chrom]; i++) 
	{
	  if(haplist[chrom][i] == hap)
	    {
	      hap_found = 1;
	      break;
	    }
	}
      if(hap_found == 0)
	{
	  haplist[chrom][no_haps[chrom]] = hap;
	  no_haps[chrom] = no_haps[chrom] + 1;
	}
    }
}

/* add haplotype to local haplotype list if not already present */
/* void add_hap_lcopy(unsigned int hap, unsigned int* hlist, int* nhaps) */
/* { */
/*   int hap_found = 0; */
/*   if(*nhaps == 0) */
/*     { */
/*       hlist[0] = hap; */
/*       *nhaps = 1; */
/*     } */
/*   else */
/*     { */
/*       for(int i=0; i < *nhaps; i++)  */
/* 	{ */
/* 	  if(hlist[i] == hap) */
/* 	    { */
/* 	      hap_found = 1; */
/* 	      break; */
/* 	    } */
/* 	} */
/*       if(hap_found == 0) */
/* 	{ */
/* 	  hlist[*nhaps] = hap; */
/* 	  *nhaps = *nhaps + 1; */
/* 	} */
/*     } */
/* } */

/* check if two haplotypes are different */
int identity2_hap(unsigned int hap1, unsigned int hap2)
{
  if(hap1==hap2)
    return 0;
  else
    return 1;
}

/* check if two haplotypes are identical */
int identity1_hap(unsigned int hap1, unsigned int hap2)
{
  if(hap1==hap2)
    return 1;
  else
    return 0;
}

/* check if sub-haplotypes match given ancestry vector */
int match_sub_hap(unsigned int hyb_hap, unsigned int pop_hap, unsigned int ancvec)
{
  if((hyb_hap & ancvec)==(pop_hap & ancvec))
    return 1;
  else
    return 0;
}

/* get marginal haplotype counts */
int get_marginal_hap_counts(unsigned int hyb_hap, int population, unsigned int ancvec, int** popY_hap_counts, unsigned int** haplist, int* no_haps, int chrom)
{
  int mhcount=0;
  for(int i=0; i<no_haps[chrom]; i++)
    if(match_sub_hap(hyb_hap,haplist[chrom][i],ancvec))
      mhcount += popY_hap_counts[chrom][haplist[chrom][i]];
  return mhcount;
}

/* get complement of ancestry vector */
/* Returns the bitwise complement of ancvec, masked to noloci bits */
unsigned int get_anc_complement(unsigned int ancvec, int noloci)
{
  unsigned int mask = (1U << noloci) - 1;  /* mask for lower noloci bits */
  return (~ancvec) & mask;
}

/* probability of switching populations */
double prPopSwitch(double r, unsigned long d)
{
  return 1.0 - cosh(r*d)*exp(-r*d); 
}

/* probability of ancestry configuration */
double prQ(unsigned int anc, int noloci, double rec_rate_bp, unsigned long* position)
{
  double t1 = 0.0;
  double t2 = 0.0;
  double Q = 0.0;
  double p = 0.0;
  double pprod = 1.0;
  unsigned int mask1 = intpow(2,noloci) - 2;
  p = prPopSwitch(rec_rate_bp,position[0]);
  if((anc & mask1) == 0)
    {
      t1 = 0.5*(1.0 - p);
      t2 = 0.5*p;
    }
  else
    {
      t1 = 0.5*p;
      t2 = 0.5*(1.0 - p);
    }
  int curr_anc = 0;
  int prev_anc = ((anc >> 0)&1);
  for(int i=1; i<noloci; i++)
    {
      curr_anc = ((anc >> i)&1);
      p = prPopSwitch(rec_rate_bp,position[i]-position[i-1]); 
      if(curr_anc != prev_anc)
	pprod *= p;
      else
	pprod *= 1.0 - p;
      prev_anc = ((anc >> i)&1);
    }
  Q = pprod*(t1 + t2);
  return Q;
}

