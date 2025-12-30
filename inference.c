#include "mongrail.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/*============================================================================
 * CACHE INITIALIZATION AND MANAGEMENT
 *
 * Pre-computes values that are reused across likelihood calculations:
 * - log P(hap|popA) and log P(hap|popB) for all haplotypes
 * - U(hap) values are cached on first computation
 * - Q(z) values computed on-the-fly with sparse enumeration (≤k switches)
 *============================================================================*/

/* Count the number of switches (transitions) in an ancestry vector.
 * A switch occurs when adjacent bits differ: z[i] != z[i-1] */
int count_switches(unsigned int z, int noloci)
{
  int switches = 0;
  for (int i = 1; i < noloci; i++) {
    int prev_bit = (z >> (i-1)) & 1;
    int curr_bit = (z >> i) & 1;
    if (prev_bit != curr_bit) {
      switches++;
    }
  }
  return switches;
}

struct likelihood_cache* cache_init(int noChr, int* no_loci, double recomb_rate,
				    unsigned long** marker_positions,
				    int** popA_hap_counts, int** popB_hap_counts,
				    unsigned int** haplist, int* no_haps,
				    int noSamplesPopA, int noSamplesPopB,
				    int max_recomb)
{
  struct likelihood_cache* cache = malloc(sizeof(struct likelihood_cache));
  cache->initialized = 1;
  cache->noChr = noChr;
  cache->no_loci = no_loci;
  cache->max_recomb = max_recomb;
  cache->recomb_rate = recomb_rate;
  cache->marker_positions = marker_positions;
  cache->log2 = log(2.0);

  /* Allocate U cache (initialized to -1 = not computed) */
  cache->U_cache = malloc(noChr * sizeof(double*));
  for (int c = 0; c < noChr; c++) {
    cache->U_cache[c] = malloc(MAXHAPS * sizeof(double));
    for (int h = 0; h < MAXHAPS; h++) {
      cache->U_cache[c][h] = -1.0;  /* sentinel for "not computed" */
    }
  }

  /* Pre-compute log P(hap|popA) and log P(hap|popB) for all haplotypes in haplist
   * Indexed by haplist position, not by haplotype value */
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
    /* Compute for haplotypes in haplist - stored at index i, not at hap value */
    for (int i = 0; i < no_haps[c]; i++) {
      unsigned int hap = haplist[c][i];
      cache->logP_A[c][i] = log_prob_single_hap(hap, popA_hap_counts, haplist,
						  no_haps, noSamplesPopA, c);
      cache->logP_B[c][i] = log_prob_single_hap(hap, popB_hap_counts, haplist,
						  no_haps, noSamplesPopB, c);
    }
  }

  /* Pre-compute U_k1 for all reference haplotypes (for scoring in preranked pruning) */
  cache->U_k1_cache = malloc(noChr * sizeof(double*));
  for (int c = 0; c < noChr; c++) {
    cache->U_k1_cache[c] = malloc(no_haps[c] * sizeof(double));
    int noloci = no_loci[c];
    unsigned long* positions = marker_positions[c];

    /* Precompute Q values for k=1 configurations once per chromosome */
    int max_configs = 2 + 2 * (noloci - 1);
    unsigned int* z_values = malloc(max_configs * sizeof(unsigned int));
    double* log_Q_values = malloc(max_configs * sizeof(double));
    int n_configs = 0;

    unsigned int all_A = 0;
    unsigned int all_B = (1U << noloci) - 1;

    /* z = all 0s (all from A) */
    double Q_z = prQ(all_A, noloci, recomb_rate, positions);
    z_values[n_configs] = all_A;
    log_Q_values[n_configs] = (Q_z > 0) ? log(Q_z) : -1e300;
    n_configs++;

    /* z = all 1s (all from B) */
    Q_z = prQ(all_B, noloci, recomb_rate, positions);
    z_values[n_configs] = all_B;
    log_Q_values[n_configs] = (Q_z > 0) ? log(Q_z) : -1e300;
    n_configs++;

    /* z with exactly 1 switch */
    for (int switch_pos = 1; switch_pos < noloci; switch_pos++) {
      unsigned int z1 = all_B << switch_pos;
      z1 &= all_B;
      Q_z = prQ(z1, noloci, recomb_rate, positions);
      z_values[n_configs] = z1;
      log_Q_values[n_configs] = (Q_z > 0) ? log(Q_z) : -1e300;
      n_configs++;

      unsigned int z2 = all_B >> (noloci - switch_pos);
      Q_z = prQ(z2, noloci, recomb_rate, positions);
      z_values[n_configs] = z2;
      log_Q_values[n_configs] = (Q_z > 0) ? log(Q_z) : -1e300;
      n_configs++;
    }

    /* Compute U_k1 for each reference haplotype using precomputed Q */
    for (int i = 0; i < no_haps[c]; i++) {
      unsigned int hap = haplist[c][i];
      double* log_terms = malloc(n_configs * sizeof(double));

      for (int j = 0; j < n_configs; j++) {
        double log_U_z = log_U_given_z(hap, z_values[j], noloci,
                                       popA_hap_counts, popB_hap_counts,
                                       haplist, no_haps, noSamplesPopA, noSamplesPopB, c);
        log_terms[j] = log_Q_values[j] + log_U_z;
      }

      /* Log-sum-exp */
      double max_log = log_terms[0];
      for (int j = 1; j < n_configs; j++) {
        if (log_terms[j] > max_log) max_log = log_terms[j];
      }
      double sum = 0.0;
      for (int j = 0; j < n_configs; j++) {
        sum += exp(log_terms[j] - max_log);
      }
      cache->U_k1_cache[c][i] = exp(max_log) * sum;

      free(log_terms);
    }

    free(z_values);
    free(log_Q_values);
  }

  return cache;
}

void cache_free(struct likelihood_cache* cache)
{
  if (!cache) return;
  for (int c = 0; c < cache->noChr; c++) {
    free(cache->U_cache[c]);
    free(cache->logP_A[c]);
    free(cache->logP_B[c]);
    free(cache->U_k1_cache[c]);
  }
  free(cache->U_cache);
  free(cache->logP_A);
  free(cache->logP_B);
  free(cache->U_k1_cache);
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
    int idx = find_hap_index(hap, haplist[chrom], no_haps[chrom]);
    if (idx >= 0 && cache->logP_A[chrom][idx] > -1e299)
      return cache->logP_A[chrom][idx];
    return log_prob_single_hap(hap, popA_hap_counts, haplist, no_haps, noSamplesPopA, chrom);
  } else if (z == all_B) {
    /* All loci from pop B - use cached value if available */
    int idx = find_hap_index(hap, haplist[chrom], no_haps[chrom]);
    if (idx >= 0 && cache->logP_B[chrom][idx] > -1e299)
      return cache->logP_B[chrom][idx];
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

/* Compute U(hap) with caching using sparse enumeration.
 * Only considers ancestry configurations with ≤ max_recomb switches.
 * Uses two-pass log-sum-exp for numerical stability. */
/* Generate z value from switch positions.
 * switch_positions[0..n_switches-1] are the positions where state changes.
 * start_state is 0 or 1 (starting ancestry).
 * Returns z with bits set according to ancestry pattern. */
static unsigned int generate_z_from_switches(int* switch_positions, int n_switches,
                                              int start_state, int noloci)
{
  unsigned int z = 0;
  int current_state = start_state;
  int switch_idx = 0;

  for (int i = 0; i < noloci; i++) {
    /* Switch happens BEFORE locus i if switch_position equals i */
    if (switch_idx < n_switches && i == switch_positions[switch_idx]) {
      current_state = 1 - current_state;
      switch_idx++;
    }
    if (current_state) z |= (1U << i);
  }
  return z;
}

/* Iterate through all z values with exactly n_switches switches.
 * Calls callback for each z value generated.
 * Returns count of z values generated. */
static int enumerate_z_with_switches(int n_switches, int noloci, int start_state,
                                      void (*callback)(unsigned int z, void* ctx),
                                      void* ctx)
{
  if (n_switches == 0) {
    unsigned int z = start_state ? ((1U << noloci) - 1) : 0;
    callback(z, ctx);
    return 1;
  }

  /* Generate all combinations of n_switches positions from {1, 2, ..., noloci-1} */
  int positions[32];  /* max switches */
  int count = 0;

  /* Initialize to first combination: {1, 2, ..., n_switches} */
  for (int i = 0; i < n_switches; i++) {
    positions[i] = i + 1;  /* switch positions are 1 to noloci-1 */
  }

  while (1) {
    unsigned int z = generate_z_from_switches(positions, n_switches, start_state, noloci);
    callback(z, ctx);
    count++;

    /* Generate next combination */
    int i = n_switches - 1;
    while (i >= 0 && positions[i] == noloci - n_switches + i) {
      i--;
    }
    if (i < 0) break;
    positions[i]++;
    for (int j = i + 1; j < n_switches; j++) {
      positions[j] = positions[j-1] + 1;
    }
  }

  return count;
}

/* Context for compute_U sparse enumeration */
struct compute_U_ctx {
  struct likelihood_cache* cache;
  unsigned int hap;
  int noloci;
  int chrom;
  double recomb_rate;
  unsigned long* positions;
  int** popA_hap_counts;
  int** popB_hap_counts;
  unsigned int** haplist;
  int* no_haps;
  int noSamplesPopA;
  int noSamplesPopB;
  double* log_terms;
  int term_count;
};

static void compute_U_callback(unsigned int z, void* ctx_ptr)
{
  struct compute_U_ctx* ctx = (struct compute_U_ctx*)ctx_ptr;
  double Q_z = prQ(z, ctx->noloci, ctx->recomb_rate, ctx->positions);
  double log_Q = (Q_z > 0) ? log(Q_z) : -1e300;
  double log_U_z = log_U_given_z_cached(ctx->cache, ctx->hap, z, ctx->noloci,
                                         ctx->popA_hap_counts, ctx->popB_hap_counts,
                                         ctx->haplist, ctx->no_haps,
                                         ctx->noSamplesPopA, ctx->noSamplesPopB,
                                         ctx->chrom);
  ctx->log_terms[ctx->term_count++] = log_Q + log_U_z;
}

double compute_U_cached(struct likelihood_cache* cache, unsigned int hap, int chrom,
			int** popA_hap_counts, int** popB_hap_counts,
			unsigned int** haplist, int* no_haps,
			int noSamplesPopA, int noSamplesPopB)
{
  /* Check cache first - indexed by haplist position */
  int hap_idx = find_hap_index(hap, haplist[chrom], no_haps[chrom]);
  if (hap_idx >= 0 && cache->U_cache[chrom][hap_idx] >= 0) {
    return cache->U_cache[chrom][hap_idx];
  }

  int noloci = cache->no_loci[chrom];
  int max_recomb = cache->max_recomb;
  double recomb_rate = cache->recomb_rate;
  unsigned long* positions = cache->marker_positions[chrom];

  /* Use full enumeration when all configs have ≤ max_recomb switches */
  int use_full_enum = (noloci <= max_recomb + 1);

  double result;

  if (use_full_enum) {
    /* Full enumeration for small noloci */
    unsigned int num_z = 1U << noloci;
    double* log_terms = malloc(num_z * sizeof(double));
    int n = 0;

    for (unsigned int z = 0; z < num_z; z++) {
      double Q_z = prQ(z, noloci, recomb_rate, positions);
      double log_Q = (Q_z > 0) ? log(Q_z) : -1e300;
      double log_U_z = log_U_given_z_cached(cache, hap, z, noloci,
					    popA_hap_counts, popB_hap_counts,
					    haplist, no_haps,
					    noSamplesPopA, noSamplesPopB, chrom);
      log_terms[n++] = log_Q + log_U_z;
    }

    /* Log-sum-exp */
    double max_log = log_terms[0];
    for (int i = 1; i < n; i++)
      if (log_terms[i] > max_log) max_log = log_terms[i];

    double sum = 0.0;
    for (int i = 0; i < n; i++)
      sum += exp(log_terms[i] - max_log);

    result = exp(max_log + log(sum));
    free(log_terms);
  } else {
    /* Sparse enumeration: generate only configs with ≤ max_recomb switches */
    /* Estimate max terms: sum of C(noloci-1, k) for k=0..max_recomb, times 2 */
    int max_terms = 2;  /* for 0 switches */
    int binom = 1;
    for (int k = 1; k <= max_recomb && k < noloci; k++) {
      binom = binom * (noloci - k) / k;
      max_terms += 2 * binom;
    }

    struct compute_U_ctx ctx;
    ctx.cache = cache;
    ctx.hap = hap;
    ctx.noloci = noloci;
    ctx.chrom = chrom;
    ctx.recomb_rate = recomb_rate;
    ctx.positions = positions;
    ctx.popA_hap_counts = popA_hap_counts;
    ctx.popB_hap_counts = popB_hap_counts;
    ctx.haplist = haplist;
    ctx.no_haps = no_haps;
    ctx.noSamplesPopA = noSamplesPopA;
    ctx.noSamplesPopB = noSamplesPopB;
    ctx.log_terms = malloc(max_terms * sizeof(double));
    ctx.term_count = 0;

    /* Enumerate all z with ≤ max_recomb switches */
    for (int n_switches = 0; n_switches <= max_recomb; n_switches++) {
      enumerate_z_with_switches(n_switches, noloci, 0, compute_U_callback, &ctx);
      enumerate_z_with_switches(n_switches, noloci, 1, compute_U_callback, &ctx);
    }

    /* Log-sum-exp */
    double max_log = ctx.log_terms[0];
    for (int i = 1; i < ctx.term_count; i++)
      if (ctx.log_terms[i] > max_log) max_log = ctx.log_terms[i];

    double sum = 0.0;
    for (int i = 0; i < ctx.term_count; i++)
      sum += exp(ctx.log_terms[i] - max_log);

    result = exp(max_log + log(sum));
    free(ctx.log_terms);
  }

  /* Cache the result - indexed by haplist position */
  if (hap_idx >= 0) {
    cache->U_cache[chrom][hap_idx] = result;
  }

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
      int idx1 = find_hap_index(hap1, haplist[i], no_haps[i]);
      int idx2 = find_hap_index(hap2, haplist[i], no_haps[i]);
      double p1_pure, p2_pure;
      if (model == MODEL_B) {
	p1_pure = (idx1 >= 0) ? exp(cache->logP_A[i][idx1]) : 0.0;
	p2_pure = (idx2 >= 0) ? exp(cache->logP_A[i][idx2]) : 0.0;
      } else {
	p1_pure = (idx1 >= 0) ? exp(cache->logP_B[i][idx1]) : 0.0;
	p2_pure = (idx2 >= 0) ? exp(cache->logP_B[i][idx2]) : 0.0;
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
	      t5 += gammln(phi + popY_hap_counts[i][j] + 1.0/no_haps[i]) -
		gammln(popY_hap_counts[i][j]+1.0/no_haps[i]);
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
	      t2A += gammln(popA_hap_counts[i][j]+1.0/no_haps[i]);
	      t2B += gammln(popB_hap_counts[i][j]+1.0/no_haps[i]);
	    }
	  for(int j=0; j<no_haps[i]; j++)
	    {
	      phi1 = identity1_hap(haplist[i][j], hap1);
	      phi2 = identity1_hap(haplist[i][j], hap2);
	      t4h1A += gammln(phi1 + popA_hap_counts[i][j] + 1.0/no_haps[i]);
	      t4h2A += gammln(phi2 + popA_hap_counts[i][j] + 1.0/no_haps[i]);
	      t4h1B += gammln(phi1 + popB_hap_counts[i][j] + 1.0/no_haps[i]);
	      t4h2B += gammln(phi2 + popB_hap_counts[i][j] + 1.0/no_haps[i]);
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
	marginal_counts[j] += pop_hap_counts[chrom][i];
	found = 1;
	break;
      }
    }
    if(!found) {
      distinct_subhaps[num_distinct] = subhap;
      marginal_counts[num_distinct] = pop_hap_counts[chrom][i];
      num_distinct++;
    }
  }
  return num_distinct;
}

/* Compute log probability of a single haplotype from population k.
 * This is log Pr(φ^c(x_i)|n_k) from the paper.
 * Uses posterior predictive distribution with Dirichlet prior.
 */
/* Cache for log_prob_single_hap gammln values */
#define SINGLE_HAP_CACHE_SIZE 256
static struct {
  int chrom;
  int** pop_hap_counts;
  int noSamplesPop;
  double cached_t1_t2_t3;
  double* gammln_base;
  double* gammln_plus1;
  int capacity;
  int valid;
} single_hap_cache[SINGLE_HAP_CACHE_SIZE] = {{0}};
static int single_hap_cache_next = 0;
#ifdef _OPENMP
#pragma omp threadprivate(single_hap_cache, single_hap_cache_next)
#endif

double log_prob_single_hap(unsigned int hap, int** pop_hap_counts,
			   unsigned int** haplist, int* no_haps,
			   int noSamplesPop, int chrom)
{
  int H = no_haps[chrom];
  double alpha = 1.0 / H;

  /* Look for cache hit */
  int slot = -1;
  for (int i = 0; i < SINGLE_HAP_CACHE_SIZE; i++) {
    if (single_hap_cache[i].valid &&
        single_hap_cache[i].pop_hap_counts == pop_hap_counts &&
        single_hap_cache[i].chrom == chrom &&
        single_hap_cache[i].noSamplesPop == noSamplesPop) {
      slot = i;
      break;
    }
  }

  /* Cache miss - compute and store gammln values */
  if (slot < 0) {
    slot = single_hap_cache_next;
    single_hap_cache_next = (single_hap_cache_next + 1) % SINGLE_HAP_CACHE_SIZE;

    /* Allocate/reallocate if needed */
    if (single_hap_cache[slot].capacity < H) {
      free(single_hap_cache[slot].gammln_base);
      free(single_hap_cache[slot].gammln_plus1);
      single_hap_cache[slot].capacity = H;
      single_hap_cache[slot].gammln_base = malloc(H * sizeof(double));
      single_hap_cache[slot].gammln_plus1 = malloc(H * sizeof(double));
    }

    double theta = 1.0 + 2.0 * noSamplesPop;
    double t1 = gammln(theta);
    double t2 = 0.0;
    double t3 = gammln(1.0 + theta);

    for (int j = 0; j < H; j++) {
      double base = pop_hap_counts[chrom][j] + alpha;
      single_hap_cache[slot].gammln_base[j] = gammln(base);
      single_hap_cache[slot].gammln_plus1[j] = gammln(base + 1.0);
      t2 += single_hap_cache[slot].gammln_base[j];
    }

    single_hap_cache[slot].cached_t1_t2_t3 = t1 - t2 - t3;
    single_hap_cache[slot].chrom = chrom;
    single_hap_cache[slot].pop_hap_counts = pop_hap_counts;
    single_hap_cache[slot].noSamplesPop = noSamplesPop;
    single_hap_cache[slot].valid = 1;
  }

  /* Use cached values to compute t4 */
  double t4 = 0.0;
  for (int j = 0; j < H; j++) {
    int phi = identity1_hap(haplist[chrom][j], hap);
    t4 += phi ? single_hap_cache[slot].gammln_plus1[j] : single_hap_cache[slot].gammln_base[j];
  }

  return single_hap_cache[slot].cached_t1_t2_t3 + t4;
}

/* LRU cache for get_distinct_subhaps results - needs to be large enough
 * to hold all unique (chrom, ancvec, population) combinations we encounter.
 * With k switches and L loci, we have O(L^k) ancestry vectors per chromosome.
 * For k=4, L=12: ~2500 configs * 2 pops = 5000 entries needed */
#define SUBHAP_CACHE_SIZE 32768
static struct {
  unsigned int ancvec;
  int chrom;
  int** pop_hap_counts;
  unsigned int* distinct_subhaps;
  int* marginal_counts;
  int num_distinct;
  int capacity;
  int valid;
  /* Cached gammln values for this entry */
  int noSamplesPop;           /* Sample count used to compute cached values */
  double cached_t1_t2_t3;     /* Precomputed t1 - t2 - t3 constant */
  double* gammln_base;        /* gammln(marginal_counts[i] + 1/m) */
  double* gammln_plus1;       /* gammln(marginal_counts[i] + 1 + 1/m) */
  int gammln_valid;           /* 1 if gammln caches are valid */
} subhap_cache[SUBHAP_CACHE_SIZE] = {{0}};
static int subhap_cache_next = 0;
#ifdef _OPENMP
#pragma omp threadprivate(subhap_cache, subhap_cache_next)
#endif

/* Get cached or compute distinct subhaps for a given ancvec */
static int get_distinct_subhaps_cached(unsigned int ancvec, int** pop_hap_counts,
                                        unsigned int** haplist, int* no_haps, int chrom,
                                        unsigned int** out_distinct, int** out_counts)
{
  /* Check cache for hit */
  for (int i = 0; i < SUBHAP_CACHE_SIZE; i++) {
    if (subhap_cache[i].valid &&
        subhap_cache[i].pop_hap_counts == pop_hap_counts &&
        subhap_cache[i].chrom == chrom &&
        subhap_cache[i].ancvec == ancvec) {
      *out_distinct = subhap_cache[i].distinct_subhaps;
      *out_counts = subhap_cache[i].marginal_counts;
      return subhap_cache[i].num_distinct;
    }
  }

  /* Cache miss - use next slot (simple round-robin) */
  int slot = subhap_cache_next;
  subhap_cache_next = (subhap_cache_next + 1) % SUBHAP_CACHE_SIZE;

  /* Invalidate and free gammln cache when evicting/reusing slot.
   * We must free the arrays because they're sized for the old num_distinct,
   * and that value will be overwritten below before log_prob_subhap checks it. */
  subhap_cache[slot].gammln_valid = 0;
  free(subhap_cache[slot].gammln_base);
  free(subhap_cache[slot].gammln_plus1);
  subhap_cache[slot].gammln_base = NULL;
  subhap_cache[slot].gammln_plus1 = NULL;

  /* Allocate/reallocate if needed */
  if (subhap_cache[slot].capacity < no_haps[chrom]) {
    free(subhap_cache[slot].distinct_subhaps);
    free(subhap_cache[slot].marginal_counts);
    subhap_cache[slot].capacity = no_haps[chrom];
    subhap_cache[slot].distinct_subhaps = malloc(subhap_cache[slot].capacity * sizeof(unsigned int));
    subhap_cache[slot].marginal_counts = malloc(subhap_cache[slot].capacity * sizeof(int));
  }

  /* Clear marginal counts */
  memset(subhap_cache[slot].marginal_counts, 0, no_haps[chrom] * sizeof(int));

  /* Compute */
  subhap_cache[slot].num_distinct = get_distinct_subhaps(ancvec, pop_hap_counts, haplist, no_haps, chrom,
                                                          subhap_cache[slot].distinct_subhaps,
                                                          subhap_cache[slot].marginal_counts);
  subhap_cache[slot].ancvec = ancvec;
  subhap_cache[slot].chrom = chrom;
  subhap_cache[slot].pop_hap_counts = pop_hap_counts;
  subhap_cache[slot].valid = 1;

  *out_distinct = subhap_cache[slot].distinct_subhaps;
  *out_counts = subhap_cache[slot].marginal_counts;
  return subhap_cache[slot].num_distinct;
}

/* Find cache slot for a given (ancvec, chrom, pop_hap_counts) and return its index */
static int find_subhap_cache_slot(unsigned int ancvec, int** pop_hap_counts, int chrom)
{
  for (int i = 0; i < SUBHAP_CACHE_SIZE; i++) {
    if (subhap_cache[i].valid &&
        subhap_cache[i].pop_hap_counts == pop_hap_counts &&
        subhap_cache[i].chrom == chrom &&
        subhap_cache[i].ancvec == ancvec) {
      return i;
    }
  }
  return -1;
}

/* Compute log probability of a sub-haplotype given ancestry mask.
 * ancvec specifies which loci are included in this sub-haplotype.
 * Uses marginal counts from the reference population.
 * Uses cached gammln values when available.
 */
double log_prob_subhap(unsigned int hap, unsigned int ancvec, int** pop_hap_counts,
		       unsigned int** haplist, int* no_haps, int noSamplesPop, int chrom)
{
  unsigned int *distinct_subhaps;
  int *marginal_counts;

  int m = get_distinct_subhaps_cached(ancvec, pop_hap_counts, haplist, no_haps, chrom,
                                       &distinct_subhaps, &marginal_counts);

  if(m == 0) {
    return -1e300;  /* No data - return very negative log prob */
  }

  /* Find the cache slot we just used/created */
  int slot = find_subhap_cache_slot(ancvec, pop_hap_counts, chrom);
  assert(slot >= 0 && "subhap cache slot lookup failed immediately after caching - cache may be too small");
  if (slot < 0) {
    /* Should not happen, but fallback to uncached computation */
    goto uncached;
  }

  /* Check if gammln values are cached for this noSamplesPop */
  if (!subhap_cache[slot].gammln_valid || subhap_cache[slot].noSamplesPop != noSamplesPop) {
    /* Compute and cache gammln values */
    double theta_prime = 1.0 + 2.0 * noSamplesPop;
    double alpha = 1.0 / m;
    double t1 = gammln(theta_prime);
    double t2 = 0.0;
    double t3 = gammln(1.0 + theta_prime);

    /* Allocate/reallocate gammln caches if needed - use num_distinct not capacity */
    if (subhap_cache[slot].gammln_base == NULL || subhap_cache[slot].num_distinct != m) {
      free(subhap_cache[slot].gammln_base);
      free(subhap_cache[slot].gammln_plus1);
      subhap_cache[slot].gammln_base = malloc(m * sizeof(double));
      subhap_cache[slot].gammln_plus1 = malloc(m * sizeof(double));
      assert(subhap_cache[slot].gammln_base != NULL && "gammln_base allocation failed");
      assert(subhap_cache[slot].gammln_plus1 != NULL && "gammln_plus1 allocation failed");
    }

    for (int i = 0; i < m; i++) {
      double base = marginal_counts[i] + alpha;
      subhap_cache[slot].gammln_base[i] = gammln(base);
      subhap_cache[slot].gammln_plus1[i] = gammln(base + 1.0);
      t2 += subhap_cache[slot].gammln_base[i];
    }

    subhap_cache[slot].cached_t1_t2_t3 = t1 - t2 - t3;
    subhap_cache[slot].noSamplesPop = noSamplesPop;
    subhap_cache[slot].gammln_valid = 1;
  }

  /* Use cached values to compute t4 */
  assert(subhap_cache[slot].gammln_valid && "Using gammln cache that was not properly initialized");
  assert(subhap_cache[slot].noSamplesPop == noSamplesPop && "gammln cache has wrong noSamplesPop");
  assert(subhap_cache[slot].gammln_base != NULL && "gammln_base is NULL when using cache");
  assert(subhap_cache[slot].gammln_plus1 != NULL && "gammln_plus1 is NULL when using cache");
  assert(subhap_cache[slot].num_distinct == m && "gammln cache size mismatch");

  unsigned int target_subhap = hap & ancvec;
  double t4 = 0.0;
  for (int i = 0; i < m; i++) {
    if (distinct_subhaps[i] == target_subhap) {
      t4 += subhap_cache[slot].gammln_plus1[i];
    } else {
      t4 += subhap_cache[slot].gammln_base[i];
    }
  }

  return subhap_cache[slot].cached_t1_t2_t3 + t4;

uncached:
  {
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
      mhcount += popY_hap_counts[chrom][i];
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

/*============================================================================
 * PRE-RANKED HAPLOTYPE PAIR PRUNING
 *
 * Instead of enumerating all 2^H compatible haplotype pairs for a hybrid
 * (where H = number of heterozygous sites), we:
 * 1. Pre-rank all haplotype pairs from reference populations by probability
 * 2. For each hybrid, iterate through pre-ranked pairs checking compatibility
 * 3. Stop when cumulative probability reaches threshold
 *
 * This is O(R^2) where R = reference haplotypes, instead of O(2^H).
 *============================================================================*/

/* Structure for a ranked haplotype pair */
struct ranked_pair {
  unsigned int hap1;
  unsigned int hap2;
  int idx1;  /* index in haplist for fast lookup */
  int idx2;
  double score;
};

/* Check if haplotype pair (h1, h2) is compatible with genotype (g1, g2).
 * Compatible means the unordered pair {h1, h2} produces the same genotype.
 * Condition: (h1 OR h2) == (g1 OR g2) AND (h1 AND h2) == (g1 AND g2)
 */
static int is_pair_compatible(unsigned int h1, unsigned int h2,
                              unsigned int g1, unsigned int g2)
{
  unsigned int g_or = g1 | g2;
  unsigned int g_and = g1 & g2;
  unsigned int h_or = h1 | h2;
  unsigned int h_and = h1 & h2;
  return (h_or == g_or) && (h_and == g_and);
}

/* Comparison function for qsort (descending order by score) */
static int compare_ranked_pairs(const void* a, const void* b) {
  const struct ranked_pair* pa = (const struct ranked_pair*)a;
  const struct ranked_pair* pb = (const struct ranked_pair*)b;
  if (pb->score > pa->score) return 1;
  if (pb->score < pa->score) return -1;
  return 0;
}

/*============================================================================
 * MAX-HEAP FOR PARTIAL SORTING
 * Instead of fully sorting all pairs, we use a max-heap to extract pairs
 * in score order until we reach the threshold. This avoids O(n log n) sort
 * when we only need a fraction of the pairs.
 *============================================================================*/

/* Heap helper: swap two pairs */
static inline void heap_swap(struct ranked_pair* a, struct ranked_pair* b)
{
  struct ranked_pair tmp = *a;
  *a = *b;
  *b = tmp;
}

/* Heap helper: sift down to maintain max-heap property */
static void heap_sift_down(struct ranked_pair* heap, int n, int i)
{
  while (1) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && heap[left].score > heap[largest].score)
      largest = left;
    if (right < n && heap[right].score > heap[largest].score)
      largest = right;

    if (largest == i) break;

    heap_swap(&heap[i], &heap[largest]);
    i = largest;
  }
}

/* Build max-heap in O(n) time */
static void heap_build(struct ranked_pair* heap, int n)
{
  for (int i = n / 2 - 1; i >= 0; i--) {
    heap_sift_down(heap, n, i);
  }
}

/* Extract max element from heap, returns the extracted pair */
static struct ranked_pair heap_extract_max(struct ranked_pair* heap, int* n)
{
  struct ranked_pair max = heap[0];
  (*n)--;
  heap[0] = heap[*n];
  heap_sift_down(heap, *n, 0);
  return max;
}

/*
 * Pre-ranked likelihood for pure models (Ma, Md).
 * Iterates through reference haplotype pairs in probability order,
 * checking compatibility with the hybrid's genotype.
 */
double lik_a_d_preranked(int indivIndex, struct indiv** hybrid_indiv,
                         int** popY_hap_counts,
                         unsigned int** haplist, int* no_haps,
                         int noSamplesPopA, int noSamplesPopB, int noChr,
                         enum modelType model, double threshold)
{
  double logL = 0.0;
  double t1, t4;
  int noSamplesPop;

  if (model == MODEL_A) {
    t1 = T1B;
    t4 = T4B;
    noSamplesPop = noSamplesPopB;
  } else {
    t1 = T1A;
    t4 = T4A;
    noSamplesPop = noSamplesPopA;
  }

  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:logL)
  #endif
  for (int c = 0; c < noChr; c++) {
    unsigned int g1 = hybrid_indiv[c][indivIndex].genotype1;
    unsigned int g2 = hybrid_indiv[c][indivIndex].genotype2;
    int numRefHaps = no_haps[c];

    /* Precompute gammln values for this chromosome */
    double alpha = 1.0 / numRefHaps;
    double* gammln_0 = malloc(numRefHaps * sizeof(double));  /* gammln(count + alpha) */
    double* gammln_1 = malloc(numRefHaps * sizeof(double));  /* gammln(count + 1 + alpha) */
    double* gammln_2 = malloc(numRefHaps * sizeof(double));  /* gammln(count + 2 + alpha) */
    for (int j = 0; j < numRefHaps; j++) {
      double base = popY_hap_counts[c][j] + alpha;
      gammln_0[j] = gammln(base);
      gammln_1[j] = gammln(base + 1.0);
      gammln_2[j] = gammln(base + 2.0);
    }

    /* Build ranked list of COMPATIBLE pairs only */
    int maxPairs = numRefHaps * numRefHaps;
    struct ranked_pair* pairs = malloc(maxPairs * sizeof(struct ranked_pair));
    double* probV = malloc(maxPairs * sizeof(double));

    double denom = 2.0 * noSamplesPop + numRefHaps;
    int numPairs = 0;
    double total_score = 0.0;

    for (int i = 0; i < numRefHaps; i++) {
      for (int j = i; j < numRefHaps; j++) {
        unsigned int h1 = haplist[c][i];
        unsigned int h2 = haplist[c][j];

        /* Skip incompatible pairs early */
        if (!is_pair_compatible(h1, h2, g1, g2)) continue;

        double p1 = (popY_hap_counts[c][i] + alpha) / denom;
        double p2 = (popY_hap_counts[c][j] + alpha) / denom;
        double score = p1 * p2;
        if (i != j) score *= 2.0;  /* count both orderings */

        pairs[numPairs].hap1 = h1;
        pairs[numPairs].hap2 = h2;
        pairs[numPairs].idx1 = i;
        pairs[numPairs].idx2 = j;
        pairs[numPairs].score = score;
        total_score += score;
        numPairs++;
      }
    }

    /* Sort only compatible pairs */
    qsort(pairs, numPairs, sizeof(struct ranked_pair), compare_ranked_pairs);

    /* Iterate through sorted compatible pairs */
    double cumulative_score = 0.0;
    int pairsUsed = 0;

    for (int p = 0; p < numPairs; p++) {
      unsigned int h1 = pairs[p].hap1;
      unsigned int h2 = pairs[p].hap2;

      /* Compute likelihood for this compatible pair using cached gammln */
      double t3 = identity2_hap(h1, h2) * L2;
      double t5 = 0.0;
      for (int j = 0; j < numRefHaps; j++) {
        int phi = identity1_hap(haplist[c][j], h1) + identity1_hap(haplist[c][j], h2);
        double g_phi = (phi == 0) ? gammln_0[j] : (phi == 1) ? gammln_1[j] : gammln_2[j];
        t5 += g_phi - gammln_0[j];
      }
      probV[pairsUsed] = exp(t1 + t3 - t4 + t5);
      pairsUsed++;

      cumulative_score += pairs[p].score;
      if (cumulative_score / total_score >= threshold) break;
    }

    /* Compute log-likelihood for this chromosome */
    double probL;
    if (pairsUsed > 0) {
      double norm_factor = max_element(probV, pairsUsed);
      for (int a = 0; a < pairsUsed; a++)
        probV[a] = exp(log(probV[a]) - log(2.0) - log(norm_factor));
      probL = log(kahanSum(probV, pairsUsed)) + log(2.0) + log(norm_factor);
    } else {
      probL = -1e300;
    }
    logL += probL;

    free(gammln_0);
    free(gammln_1);
    free(gammln_2);
    free(pairs);
    free(probV);
  }

  return logL;
}

/*
 * Pre-ranked likelihood for F1 model (Mc).
 * Iterates through pairs of (popA haplotype, popB haplotype) in probability order.
 * Filters compatible pairs first, then sorts only those.
 */
double lik_c_preranked(int indivIndex, struct indiv** hybrid_indiv,
                       int** popB_hap_counts, int** popA_hap_counts,
                       unsigned int** haplist, int* no_haps,
                       int noSamplesPopB, int noSamplesPopA, int noChr,
                       double threshold)
{
  double logL = 0.0;

  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:logL)
  #endif
  for (int c = 0; c < noChr; c++) {
    unsigned int g1 = hybrid_indiv[c][indivIndex].genotype1;
    unsigned int g2 = hybrid_indiv[c][indivIndex].genotype2;
    int numRefHaps = no_haps[c];
    double alpha = 1.0 / numRefHaps;

    /* Precompute gammln values for this chromosome */
    double* gammln_A_0 = malloc(numRefHaps * sizeof(double));
    double* gammln_A_1 = malloc(numRefHaps * sizeof(double));
    double* gammln_B_0 = malloc(numRefHaps * sizeof(double));
    double* gammln_B_1 = malloc(numRefHaps * sizeof(double));
    double t2A = 0.0, t2B = 0.0;  /* constant sums for this chromosome */
    for (int k = 0; k < numRefHaps; k++) {
      double baseA = popA_hap_counts[c][k] + alpha;
      double baseB = popB_hap_counts[c][k] + alpha;
      gammln_A_0[k] = gammln(baseA);
      gammln_A_1[k] = gammln(baseA + 1.0);
      gammln_B_0[k] = gammln(baseB);
      gammln_B_1[k] = gammln(baseB + 1.0);
      t2A += gammln_A_0[k];
      t2B += gammln_B_0[k];
    }

    /* Build ranked list of COMPATIBLE pairs only */
    int maxPairs = numRefHaps * numRefHaps;
    struct ranked_pair* pairs = malloc(maxPairs * sizeof(struct ranked_pair));
    double* probV = malloc(maxPairs * sizeof(double));

    double denomA = 2.0 * noSamplesPopA + numRefHaps;
    double denomB = 2.0 * noSamplesPopB + numRefHaps;
    int numPairs = 0;
    double total_score = 0.0;

    /* For F1: score = P(h1|A)*P(h2|B) + P(h1|B)*P(h2|A) */
    for (int i = 0; i < numRefHaps; i++) {
      for (int j = i; j < numRefHaps; j++) {
        unsigned int h1 = haplist[c][i];
        unsigned int h2 = haplist[c][j];

        /* Skip incompatible pairs early */
        if (!is_pair_compatible(h1, h2, g1, g2)) continue;

        double p1_A = (popA_hap_counts[c][i] + alpha) / denomA;
        double p2_A = (popA_hap_counts[c][j] + alpha) / denomA;
        double p1_B = (popB_hap_counts[c][i] + alpha) / denomB;
        double p2_B = (popB_hap_counts[c][j] + alpha) / denomB;

        double score = p1_A * p2_B + p1_B * p2_A;
        if (i != j) score *= 2.0;  /* count both orderings for het pairs */

        pairs[numPairs].hap1 = h1;
        pairs[numPairs].hap2 = h2;
        pairs[numPairs].idx1 = i;
        pairs[numPairs].idx2 = j;
        pairs[numPairs].score = score;
        total_score += score;
        numPairs++;
      }
    }

    /* Sort only compatible pairs */
    qsort(pairs, numPairs, sizeof(struct ranked_pair), compare_ranked_pairs);

    /* Iterate through sorted compatible pairs */
    double cumulative_score = 0.0;
    int pairsUsed = 0;

    for (int p = 0; p < numPairs; p++) {
      unsigned int h1 = pairs[p].hap1;
      unsigned int h2 = pairs[p].hap2;

      /* Compute F1 likelihood using cached gammln values */
      double t4h1A = 0.0, t4h2A = 0.0, t4h1B = 0.0, t4h2B = 0.0;
      for (int k = 0; k < numRefHaps; k++) {
        int phi1 = identity1_hap(haplist[c][k], h1);
        int phi2 = identity1_hap(haplist[c][k], h2);
        t4h1A += phi1 ? gammln_A_1[k] : gammln_A_0[k];
        t4h2A += phi2 ? gammln_A_1[k] : gammln_A_0[k];
        t4h1B += phi1 ? gammln_B_1[k] : gammln_B_0[k];
        t4h2B += phi2 ? gammln_B_1[k] : gammln_B_0[k];
      }

      double p1A = exp(T1A - t2A - TC4A + t4h1A);
      double p1B = exp(T1B - t2B - TC4B + t4h1B);
      double p2A = exp(T1A - t2A - TC4A + t4h2A);
      double p2B = exp(T1B - t2B - TC4B + t4h2B);

      probV[pairsUsed] = identity2_hap(h1, h2) * (p1A*p2B + p2A*p1B) +
                         (1 - identity2_hap(h1, h2)) * (p1A*p2B);
      pairsUsed++;

      cumulative_score += pairs[p].score;
      if (cumulative_score / total_score >= threshold) break;
    }

    /* Compute log-likelihood for this chromosome */
    double probL;
    if (pairsUsed > 0) {
      double norm_factor = max_element(probV, pairsUsed);
      for (int a = 0; a < pairsUsed; a++)
        probV[a] = exp(log(probV[a]) - log(2.0) - log(norm_factor));
      probL = log(kahanSum(probV, pairsUsed)) + log(2.0) + log(norm_factor);
    } else {
      probL = -1e300;
    }
    logL += probL;

    free(gammln_A_0);
    free(gammln_A_1);
    free(gammln_B_0);
    free(gammln_B_1);
    free(pairs);
    free(probV);
  }

  return logL;
}

/*
 * Structure to hold precomputed Q values for k=1 ancestry configurations.
 * This avoids recomputing Q(z) for each haplotype.
 */
struct precomputed_Q_k1 {
  int n_configs;
  unsigned int* z_values;   /* ancestry configurations */
  double* log_Q_values;     /* log Q(z) for each configuration */
};

/*
 * Precompute Q values for all k=1 ancestry configurations.
 * These only depend on marker positions and recombination rate, not on haplotypes.
 */
static struct precomputed_Q_k1* precompute_Q_k1(int noloci, double recomb_rate,
                                                 unsigned long* positions)
{
  struct precomputed_Q_k1* pq = malloc(sizeof(struct precomputed_Q_k1));
  int max_configs = 2 + 2 * (noloci - 1);
  pq->z_values = malloc(max_configs * sizeof(unsigned int));
  pq->log_Q_values = malloc(max_configs * sizeof(double));
  pq->n_configs = 0;

  unsigned int all_A = 0;
  unsigned int all_B = (1U << noloci) - 1;

  /* z = all 0s (all from A) */
  double Q_z = prQ(all_A, noloci, recomb_rate, positions);
  pq->z_values[pq->n_configs] = all_A;
  pq->log_Q_values[pq->n_configs] = (Q_z > 0) ? log(Q_z) : -1e300;
  pq->n_configs++;

  /* z = all 1s (all from B) */
  Q_z = prQ(all_B, noloci, recomb_rate, positions);
  pq->z_values[pq->n_configs] = all_B;
  pq->log_Q_values[pq->n_configs] = (Q_z > 0) ? log(Q_z) : -1e300;
  pq->n_configs++;

  /* z with exactly 1 switch */
  for (int switch_pos = 1; switch_pos < noloci; switch_pos++) {
    /* Start with A, switch to B */
    unsigned int z1 = all_B << switch_pos;
    z1 &= all_B;
    Q_z = prQ(z1, noloci, recomb_rate, positions);
    pq->z_values[pq->n_configs] = z1;
    pq->log_Q_values[pq->n_configs] = (Q_z > 0) ? log(Q_z) : -1e300;
    pq->n_configs++;

    /* Start with B, switch to A */
    unsigned int z2 = all_B >> (noloci - switch_pos);
    Q_z = prQ(z2, noloci, recomb_rate, positions);
    pq->z_values[pq->n_configs] = z2;
    pq->log_Q_values[pq->n_configs] = (Q_z > 0) ? log(Q_z) : -1e300;
    pq->n_configs++;
  }

  return pq;
}

static void free_precomputed_Q_k1(struct precomputed_Q_k1* pq)
{
  if (pq) {
    free(pq->z_values);
    free(pq->log_Q_values);
    free(pq);
  }
}

/*
 * Compute U_k1 for a haplotype using precomputed Q values.
 * This is much faster than compute_U_k1_simple when computing for many haplotypes.
 */
static double compute_U_k1_with_precomputed_Q(unsigned int hap, int noloci,
                                               struct precomputed_Q_k1* pq,
                                               int** popA_hap_counts, int** popB_hap_counts,
                                               unsigned int** haplist, int* no_haps,
                                               int noSamplesPopA, int noSamplesPopB, int chrom)
{
  double* log_terms = malloc(pq->n_configs * sizeof(double));

  for (int i = 0; i < pq->n_configs; i++) {
    unsigned int z = pq->z_values[i];
    double log_U_z = log_U_given_z(hap, z, noloci, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB, chrom);
    log_terms[i] = pq->log_Q_values[i] + log_U_z;
  }

  /* Log-sum-exp */
  double max_log = log_terms[0];
  for (int i = 1; i < pq->n_configs; i++) {
    if (log_terms[i] > max_log) max_log = log_terms[i];
  }

  double sum = 0.0;
  for (int i = 0; i < pq->n_configs; i++) {
    sum += exp(log_terms[i] - max_log);
  }

  free(log_terms);
  return exp(max_log) * sum;
}

/*
 * Pre-ranked likelihood for backcross models (Mb, Me).
 * Iterates through reference haplotype pairs in probability order.
 * For backcross: one haplotype is pure, other is recombinant.
 */
double lik_b_e_preranked(struct likelihood_cache* cache, int indivIndex,
                         struct indiv** hybrid_indiv,
                         int** popB_hap_counts, int** popA_hap_counts,
                         unsigned int** haplist, int* no_haps,
                         int noSamplesPopB, int noSamplesPopA, int noChr,
                         enum modelType model, double threshold)
{
  double logL = 0.0;

  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:logL)
  #endif
  for (int c = 0; c < noChr; c++) {
    unsigned int g1 = hybrid_indiv[c][indivIndex].genotype1;
    unsigned int g2 = hybrid_indiv[c][indivIndex].genotype2;
    int numRefHaps = no_haps[c];

    /* Use cached U_k1 and P_pure values - no recomputation needed */
    double* U_k1_values = cache->U_k1_cache[c];
    double* P_pure_values = malloc(numRefHaps * sizeof(double));

    for (int i = 0; i < numRefHaps; i++) {
      /* Use cached logP values (proper posterior predictive formula) */
      if (model == MODEL_B) {
        P_pure_values[i] = exp(cache->logP_A[c][i]);
      } else {
        P_pure_values[i] = exp(cache->logP_B[c][i]);
      }
    }

    /* Build ranked list of COMPATIBLE pairs only */
    int maxPairs = numRefHaps * numRefHaps;
    struct ranked_pair* pairs = malloc(maxPairs * sizeof(struct ranked_pair));
    double* probV = malloc(maxPairs * sizeof(double));

    int numPairs = 0;
    double total_score = 0.0;

    for (int i = 0; i < numRefHaps; i++) {
      for (int j = i; j < numRefHaps; j++) {
        unsigned int h1 = haplist[c][i];
        unsigned int h2 = haplist[c][j];

        /* Skip incompatible pairs early */
        if (!is_pair_compatible(h1, h2, g1, g2)) continue;

        /* Backcross: one pure, one recombinant */
        double score1 = P_pure_values[i] * U_k1_values[j];
        double score2 = P_pure_values[j] * U_k1_values[i];
        double score = (score1 > score2) ? score1 : score2;
        if (i != j) score *= 2.0;

        pairs[numPairs].hap1 = h1;
        pairs[numPairs].hap2 = h2;
        pairs[numPairs].idx1 = i;
        pairs[numPairs].idx2 = j;
        pairs[numPairs].score = score;
        total_score += score;
        numPairs++;
      }
    }

    /* Sort only compatible pairs */
    qsort(pairs, numPairs, sizeof(struct ranked_pair), compare_ranked_pairs);

    /* Iterate through sorted compatible pairs */
    double cumulative_score = 0.0;
    int pairsUsed = 0;

    for (int p = 0; p < numPairs; p++) {
      unsigned int h1 = pairs[p].hap1;
      unsigned int h2 = pairs[p].hap2;
      int idx1 = pairs[p].idx1;
      int idx2 = pairs[p].idx2;

      /* Compute actual U values using cache (full k recombinations) */
      double U1 = compute_U_cached(cache, h1, c, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB);
      double U2 = compute_U_cached(cache, h2, c, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB);

      /* Get P_pure using stored indices */
      double p1_pure = P_pure_values[idx1];
      double p2_pure = P_pure_values[idx2];

      /* Backcross likelihood: match original formula exactly */
      if (identity2_hap(h1, h2)) {
        probV[pairsUsed] = p1_pure * U2 + p2_pure * U1;
      } else {
        probV[pairsUsed] = p1_pure * U2;
      }
      pairsUsed++;

      cumulative_score += pairs[p].score;
      if (cumulative_score / total_score >= threshold) break;
    }

    /* Compute log-likelihood for this chromosome */
    double probL;
    if (pairsUsed > 0) {
      double norm_factor = max_element(probV, pairsUsed);
      if (norm_factor > 0) {
        double log_norm = log(norm_factor);
        for (int a = 0; a < pairsUsed; a++)
          probV[a] = exp(log(probV[a]) - cache->log2 - log_norm);
        probL = log(kahanSum(probV, pairsUsed)) + cache->log2 + log_norm;
      } else {
        probL = -1e300;
      }
    } else {
      probL = -1e300;
    }
    logL += probL;

    free(P_pure_values);
    free(pairs);
    free(probV);
  }

  return logL;
}

/*
 * Pre-ranked likelihood for F2 model (Mf).
 * Iterates through reference haplotype pairs in probability order.
 * For F2: both haplotypes are recombinant (have U values).
 */
double lik_f_preranked(struct likelihood_cache* cache, int indivIndex,
                       struct indiv** hybrid_indiv,
                       int** popB_hap_counts, int** popA_hap_counts,
                       unsigned int** haplist, int* no_haps,
                       int noSamplesPopB, int noSamplesPopA, int noChr,
                       double threshold)
{
  double logL = 0.0;

  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:logL)
  #endif
  for (int c = 0; c < noChr; c++) {
    unsigned int g1 = hybrid_indiv[c][indivIndex].genotype1;
    unsigned int g2 = hybrid_indiv[c][indivIndex].genotype2;
    int numRefHaps = no_haps[c];

    /* Use cached U_k1 values - no recomputation needed */
    double* U_k1_values = cache->U_k1_cache[c];

    /* Build ranked list of COMPATIBLE pairs only */
    int maxPairs = numRefHaps * numRefHaps;
    struct ranked_pair* pairs = malloc(maxPairs * sizeof(struct ranked_pair));
    double* probV = malloc(maxPairs * sizeof(double));

    int numPairs = 0;
    double total_score = 0.0;

    for (int i = 0; i < numRefHaps; i++) {
      for (int j = i; j < numRefHaps; j++) {
        unsigned int h1 = haplist[c][i];
        unsigned int h2 = haplist[c][j];

        /* Skip incompatible pairs early */
        if (!is_pair_compatible(h1, h2, g1, g2)) continue;

        double score = U_k1_values[i] * U_k1_values[j];
        if (i != j) score *= 2.0;

        pairs[numPairs].hap1 = h1;
        pairs[numPairs].hap2 = h2;
        pairs[numPairs].idx1 = i;
        pairs[numPairs].idx2 = j;
        pairs[numPairs].score = score;
        total_score += score;
        numPairs++;
      }
    }

    /* Sort only compatible pairs */
    qsort(pairs, numPairs, sizeof(struct ranked_pair), compare_ranked_pairs);

    /* Iterate through sorted compatible pairs */
    double cumulative_score = 0.0;
    int pairsUsed = 0;

    for (int p = 0; p < numPairs; p++) {
      unsigned int h1 = pairs[p].hap1;
      unsigned int h2 = pairs[p].hap2;

      /* Compute actual U values using cache */
      double U1 = compute_U_cached(cache, h1, c, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB);
      double U2 = compute_U_cached(cache, h2, c, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB);

      if (identity2_hap(h1, h2)) {
        probV[pairsUsed] = 2.0 * U1 * U2;
      } else {
        probV[pairsUsed] = U1 * U2;
      }
      pairsUsed++;

      cumulative_score += pairs[p].score;
      if (cumulative_score / total_score >= threshold) break;
    }

    /* Compute log-likelihood for this chromosome */
    double probL;
    if (pairsUsed > 0) {
      double norm_factor = max_element(probV, pairsUsed);
      if (norm_factor > 0) {
        double log_norm = log(norm_factor);
        for (int a = 0; a < pairsUsed; a++)
          probV[a] = exp(log(probV[a]) - cache->log2 - log_norm);
        probL = log(kahanSum(probV, pairsUsed)) + cache->log2 + log_norm;
      } else {
        probL = -1e300;
      }
    } else {
      probL = -1e300;
    }
    logL += probL;

    free(pairs);
    free(probV);
  }

  return logL;
}

/*============================================================================
 * HAPLOTYPE PAIR SCORING AND PRUNING (OLD APPROACH - KEPT FOR REFERENCE)
 *
 * For unphased data with many heterozygous sites, we can avoid computing
 * likelihoods for implausible haplotype phasings by:
 * 1. Scoring each compatible haplotype pair by posterior probability
 * 2. Sorting pairs by score
 * 3. Summing only pairs that contribute to cumulative probability threshold
 *============================================================================*/

/* Structure for sorting haplotype pairs by score */
struct scored_pair {
  int pair_index;      /* index into compHaps (k for pair k/2) */
  double score;        /* posterior probability score */
};

/* Comparison function for qsort (descending order by score) */
static int compare_scored_pairs(const void* a, const void* b) {
  const struct scored_pair* pa = (const struct scored_pair*)a;
  const struct scored_pair* pb = (const struct scored_pair*)b;
  if (pb->score > pa->score) return 1;
  if (pb->score < pa->score) return -1;
  return 0;
}

/*
 * Score a haplotype pair for pure model (Ma or Md).
 * Score = Pr(h1|pop) * Pr(h2|pop) using posterior predictive distribution.
 */
static double score_diplotype_pure(unsigned int hap1, unsigned int hap2,
                                   int* pop_hap_counts,
                                   unsigned int* haplist, int no_haps,
                                   int noSamplesPop)
{
  int idx1 = find_hap_index(hap1, haplist, no_haps);
  int idx2 = find_hap_index(hap2, haplist, no_haps);

  double denom = 2.0 * noSamplesPop + no_haps;
  double p1 = (idx1 >= 0 ? pop_hap_counts[idx1] + 1.0 : 1.0) / denom;
  double p2 = (idx2 >= 0 ? pop_hap_counts[idx2] + 1.0 : 1.0) / denom;

  return p1 * p2;
}

/*
 * Score a haplotype pair for F1 model (Mc).
 * Score = Pr(h1|A)*Pr(h2|B) + Pr(h1|B)*Pr(h2|A)
 */
static double score_diplotype_f1(unsigned int hap1, unsigned int hap2,
                                 int* popA_hap_counts, int* popB_hap_counts,
                                 unsigned int* haplist, int no_haps,
                                 int noSamplesPopA, int noSamplesPopB)
{
  int idx1 = find_hap_index(hap1, haplist, no_haps);
  int idx2 = find_hap_index(hap2, haplist, no_haps);

  double denomA = 2.0 * noSamplesPopA + no_haps;
  double denomB = 2.0 * noSamplesPopB + no_haps;

  double p1_A = (idx1 >= 0 ? popA_hap_counts[idx1] + 1.0 : 1.0) / denomA;
  double p2_A = (idx2 >= 0 ? popA_hap_counts[idx2] + 1.0 : 1.0) / denomA;
  double p1_B = (idx1 >= 0 ? popB_hap_counts[idx1] + 1.0 : 1.0) / denomB;
  double p2_B = (idx2 >= 0 ? popB_hap_counts[idx2] + 1.0 : 1.0) / denomB;

  /* Sum over both parental orderings */
  return p1_A * p2_B + p1_B * p2_A;
}

/*
 * Combined score for pruning: max over Ma, Md, Mc models.
 * This ensures we don't prune pairs that are important for any model.
 */
static double score_diplotype_combined(unsigned int hap1, unsigned int hap2,
                                       int* popA_hap_counts, int* popB_hap_counts,
                                       unsigned int* haplist, int no_haps,
                                       int noSamplesPopA, int noSamplesPopB)
{
  double score_Ma = score_diplotype_pure(hap1, hap2, popA_hap_counts,
                                         haplist, no_haps, noSamplesPopA);
  double score_Md = score_diplotype_pure(hap1, hap2, popB_hap_counts,
                                         haplist, no_haps, noSamplesPopB);
  double score_Mc = score_diplotype_f1(hap1, hap2, popA_hap_counts, popB_hap_counts,
                                       haplist, no_haps, noSamplesPopA, noSamplesPopB);

  double max_score = score_Ma;
  if (score_Md > max_score) max_score = score_Md;
  if (score_Mc > max_score) max_score = score_Mc;

  return max_score;
}

/*
 * Pruned likelihood for pure models (Ma, Md).
 * Only sums over haplotype pairs up to cumulative probability threshold.
 */
double lik_a_d_pruned(int indivIndex, struct indiv** hybrid_indiv,
                      int** popY_hap_counts,
                      int** popA_hap_counts, int** popB_hap_counts,
                      unsigned int** haplist, int* no_haps,
                      int noSamplesPopA, int noSamplesPopB, int noChr,
                      enum modelType model, double threshold)
{
  unsigned int hap1, hap2;
  double logL = 0.0;
  double probL = 0.0;
  double norm_factor;
  double *probV = (double *)malloc(MAXHAPS * sizeof(double));
  struct scored_pair *pairs = (struct scored_pair *)malloc(MAXHAPS * sizeof(struct scored_pair));

  double t1, t4;
  if (model == MODEL_A) {
    t1 = T1B;
    t4 = T4B;
  } else {
    t1 = T1A;
    t4 = T4A;
  }

  for (int i = 0; i < noChr; i++) {
    int numPairs = hybrid_indiv[i][indivIndex].numHaps / 2;

    /* If only 1-2 pairs, no point in pruning */
    if (numPairs <= 2) {
      /* Use original calculation */
      probL = 0.0;
      for (int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k += 2) {
        hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
        hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

        double t3 = identity2_hap(hap1, hap2) * L2;
        double t5 = 0.0;
        for (int j = 0; j < no_haps[i]; j++) {
          int phi = identity1_hap(haplist[i][j], hap1) + identity1_hap(haplist[i][j], hap2);
          t5 += gammln(phi + popY_hap_counts[i][j] + 1.0/no_haps[i]) -
                gammln(popY_hap_counts[i][j] + 1.0/no_haps[i]);
        }
        probV[k/2] = exp(t1 + t3 - t4 + t5);
      }
      norm_factor = max_element(probV, numPairs);
      for (int a = 0; a < numPairs; a++)
        probV[a] = exp(log(probV[a]) - log(2.0) - log(norm_factor));
      probL = log(kahanSum(probV, numPairs)) + log(2.0) + log(norm_factor);
      logL += probL;
      continue;
    }

    /* Score and sort haplotype pairs */
    double total_score = 0.0;
    for (int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k += 2) {
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      double score = score_diplotype_combined(hap1, hap2,
                                              popA_hap_counts[i], popB_hap_counts[i],
                                              haplist[i], no_haps[i],
                                              noSamplesPopA, noSamplesPopB);
      pairs[k/2].pair_index = k;
      pairs[k/2].score = score;
      total_score += score;
    }

    /* Sort by score descending */
    qsort(pairs, numPairs, sizeof(struct scored_pair), compare_scored_pairs);

    /* Compute likelihoods for pairs until cumulative score exceeds threshold */
    double cumulative_score = 0.0;
    int pairsUsed = 0;

    for (int p = 0; p < numPairs; p++) {
      int k = pairs[p].pair_index;
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      double t3 = identity2_hap(hap1, hap2) * L2;
      double t5 = 0.0;
      for (int j = 0; j < no_haps[i]; j++) {
        int phi = identity1_hap(haplist[i][j], hap1) + identity1_hap(haplist[i][j], hap2);
        t5 += gammln(phi + popY_hap_counts[i][j] + 1.0/no_haps[i]) -
              gammln(popY_hap_counts[i][j] + 1.0/no_haps[i]);
      }
      probV[pairsUsed] = exp(t1 + t3 - t4 + t5);
      pairsUsed++;

      cumulative_score += pairs[p].score;
      if (cumulative_score / total_score >= threshold) break;
    }

    norm_factor = max_element(probV, pairsUsed);
    for (int a = 0; a < pairsUsed; a++)
      probV[a] = exp(log(probV[a]) - log(2.0) - log(norm_factor));
    probL = log(kahanSum(probV, pairsUsed)) + log(2.0) + log(norm_factor);
    logL += probL;
  }

  free(probV);
  free(pairs);
  return logL;
}

/*
 * Pruned likelihood for F1 model (Mc).
 * Only sums over haplotype pairs up to cumulative probability threshold.
 */
double lik_c_pruned(int indivIndex, struct indiv** hybrid_indiv,
                    int** popB_hap_counts, int** popA_hap_counts,
                    unsigned int** haplist, int* no_haps,
                    int noSamplesPopB, int noSamplesPopA, int noChr,
                    double threshold)
{
  unsigned int hap1, hap2;
  double logL = 0.0;
  double probL = 0.0;
  double norm_factor;
  double *probV = (double *)malloc(MAXHAPS * sizeof(double));
  struct scored_pair *pairs = (struct scored_pair *)malloc(MAXHAPS * sizeof(struct scored_pair));
  double t2A, t2B;
  double t4h1A, t4h2A, t4h1B, t4h2B;
  double p1A, p1B, p2A, p2B;

  for (int i = 0; i < noChr; i++) {
    int numPairs = hybrid_indiv[i][indivIndex].numHaps / 2;

    /* If only 1-2 pairs, no point in pruning */
    if (numPairs <= 2) {
      probL = 0.0;
      for (int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k += 2) {
        hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
        hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

        /* First compute t2A, t2B (baseline sums) */
        t2A = 0.0;
        t2B = 0.0;
        for (int j = 0; j < no_haps[i]; j++) {
          t2A += gammln(popA_hap_counts[i][j] + 1.0/no_haps[i]);
          t2B += gammln(popB_hap_counts[i][j] + 1.0/no_haps[i]);
        }

        /* Then compute t4h1/t4h2 (sums including haplotype contributions) */
        t4h1A = 0.0; t4h2A = 0.0; t4h1B = 0.0; t4h2B = 0.0;
        for (int j = 0; j < no_haps[i]; j++) {
          int phi1 = identity1_hap(haplist[i][j], hap1);
          int phi2 = identity1_hap(haplist[i][j], hap2);
          t4h1A += gammln(phi1 + popA_hap_counts[i][j] + 1.0/no_haps[i]);
          t4h2A += gammln(phi2 + popA_hap_counts[i][j] + 1.0/no_haps[i]);
          t4h1B += gammln(phi1 + popB_hap_counts[i][j] + 1.0/no_haps[i]);
          t4h2B += gammln(phi2 + popB_hap_counts[i][j] + 1.0/no_haps[i]);
        }

        p1A = exp(T1A - t2A - TC4A + t4h1A);
        p1B = exp(T1B - t2B - TC4B + t4h1B);
        p2A = exp(T1A - t2A - TC4A + t4h2A);
        p2B = exp(T1B - t2B - TC4B + t4h2B);
        probV[k/2] = identity2_hap(hap1, hap2) * (p1A*p2B + p2A*p1B) +
                     (1 - identity2_hap(hap1, hap2)) * (p1A*p2B);
      }
      norm_factor = max_element(probV, numPairs);
      for (int a = 0; a < numPairs; a++)
        probV[a] = exp(log(probV[a]) - log(2.0) - log(norm_factor));
      probL = log(kahanSum(probV, numPairs)) + log(2.0) + log(norm_factor);
      logL += probL;
      continue;
    }

    /* Score and sort haplotype pairs */
    double total_score = 0.0;
    for (int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k += 2) {
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      double score = score_diplotype_combined(hap1, hap2,
                                              popA_hap_counts[i], popB_hap_counts[i],
                                              haplist[i], no_haps[i],
                                              noSamplesPopA, noSamplesPopB);
      pairs[k/2].pair_index = k;
      pairs[k/2].score = score;
      total_score += score;
    }

    /* Sort by score descending */
    qsort(pairs, numPairs, sizeof(struct scored_pair), compare_scored_pairs);

    /* Compute likelihoods for pairs until cumulative score exceeds threshold */
    double cumulative_score = 0.0;
    int pairsUsed = 0;

    for (int p = 0; p < numPairs; p++) {
      int k = pairs[p].pair_index;
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      /* First compute t2A, t2B (baseline sums) */
      t2A = 0.0;
      t2B = 0.0;
      for (int j = 0; j < no_haps[i]; j++) {
        t2A += gammln(popA_hap_counts[i][j] + 1.0/no_haps[i]);
        t2B += gammln(popB_hap_counts[i][j] + 1.0/no_haps[i]);
      }

      /* Then compute t4h1/t4h2 (sums including haplotype contributions) */
      t4h1A = 0.0; t4h2A = 0.0; t4h1B = 0.0; t4h2B = 0.0;
      for (int j = 0; j < no_haps[i]; j++) {
        int phi1 = identity1_hap(haplist[i][j], hap1);
        int phi2 = identity1_hap(haplist[i][j], hap2);
        t4h1A += gammln(phi1 + popA_hap_counts[i][j] + 1.0/no_haps[i]);
        t4h2A += gammln(phi2 + popA_hap_counts[i][j] + 1.0/no_haps[i]);
        t4h1B += gammln(phi1 + popB_hap_counts[i][j] + 1.0/no_haps[i]);
        t4h2B += gammln(phi2 + popB_hap_counts[i][j] + 1.0/no_haps[i]);
      }

      p1A = exp(T1A - t2A - TC4A + t4h1A);
      p1B = exp(T1B - t2B - TC4B + t4h1B);
      p2A = exp(T1A - t2A - TC4A + t4h2A);
      p2B = exp(T1B - t2B - TC4B + t4h2B);
      probV[pairsUsed] = identity2_hap(hap1, hap2) * (p1A*p2B + p2A*p1B) +
                         (1 - identity2_hap(hap1, hap2)) * (p1A*p2B);
      pairsUsed++;

      cumulative_score += pairs[p].score;
      if (cumulative_score / total_score >= threshold) break;
    }

    norm_factor = max_element(probV, pairsUsed);
    for (int a = 0; a < pairsUsed; a++)
      probV[a] = exp(log(probV[a]) - log(2.0) - log(norm_factor));
    probL = log(kahanSum(probV, pairsUsed)) + log(2.0) + log(norm_factor);
    logL += probL;
  }

  free(probV);
  free(pairs);
  return logL;
}

/*
 * Compute U(h) using sparse enumeration with k=1 (max 1 recombination).
 * Only considers ancestry configurations with at most 1 switch:
 * - z = 0 (all from A)
 * - z = 2^L - 1 (all from B)
 * - z with single switch at each position
 * This is O(L) configurations instead of O(2^L).
 */
static double compute_U_k1(unsigned int hap, int noloci, double recomb_rate,
                           unsigned long* positions,
                           int** popA_hap_counts, int** popB_hap_counts,
                           unsigned int** haplist, int* no_haps,
                           int noSamplesPopA, int noSamplesPopB, int chrom)
{
  /* We'll compute log terms first for numerical stability */
  int max_configs = 2 + 2 * (noloci - 1);  /* 2L configurations max */
  double* log_terms = (double*)malloc(max_configs * sizeof(double));
  int n_configs = 0;

  unsigned int all_A = 0;
  unsigned int all_B = (1U << noloci) - 1;

  /* z = all 0s (all from A, no recombination) */
  double Q_z = prQ(all_A, noloci, recomb_rate, positions);
  double log_Q = (Q_z > 0) ? log(Q_z) : -1e300;
  double log_U_z = log_U_given_z(hap, all_A, noloci, popA_hap_counts, popB_hap_counts,
                                 haplist, no_haps, noSamplesPopA, noSamplesPopB, chrom);
  log_terms[n_configs++] = log_Q + log_U_z;

  /* z = all 1s (all from B, no recombination) */
  Q_z = prQ(all_B, noloci, recomb_rate, positions);
  log_Q = (Q_z > 0) ? log(Q_z) : -1e300;
  log_U_z = log_U_given_z(hap, all_B, noloci, popA_hap_counts, popB_hap_counts,
                          haplist, no_haps, noSamplesPopA, noSamplesPopB, chrom);
  log_terms[n_configs++] = log_Q + log_U_z;

  /* z with exactly 1 switch */
  for (int switch_pos = 1; switch_pos < noloci; switch_pos++) {
    /* Start with A (0s), switch to B (1s) at position switch_pos */
    unsigned int z = all_B ^ ((1U << switch_pos) - 1);
    Q_z = prQ(z, noloci, recomb_rate, positions);
    log_Q = (Q_z > 0) ? log(Q_z) : -1e300;
    log_U_z = log_U_given_z(hap, z, noloci, popA_hap_counts, popB_hap_counts,
                            haplist, no_haps, noSamplesPopA, noSamplesPopB, chrom);
    log_terms[n_configs++] = log_Q + log_U_z;

    /* Start with B (1s), switch to A (0s) at position switch_pos */
    z = (1U << switch_pos) - 1;
    Q_z = prQ(z, noloci, recomb_rate, positions);
    log_Q = (Q_z > 0) ? log(Q_z) : -1e300;
    log_U_z = log_U_given_z(hap, z, noloci, popA_hap_counts, popB_hap_counts,
                            haplist, no_haps, noSamplesPopA, noSamplesPopB, chrom);
    log_terms[n_configs++] = log_Q + log_U_z;
  }

  /* Sum using log-sum-exp trick */
  double max_log = log_terms[0];
  for (int i = 1; i < n_configs; i++) {
    if (log_terms[i] > max_log) max_log = log_terms[i];
  }

  double sum = 0.0;
  for (int i = 0; i < n_configs; i++) {
    sum += exp(log_terms[i] - max_log);
  }

  free(log_terms);
  return exp(max_log + log(sum));
}

/*
 * Score a haplotype pair for backcross model (Mb or Me).
 * For Mb: one haplotype from pure A, other is recombinant
 * Score = max(P_A(h1) * U_k1(h2), P_A(h2) * U_k1(h1))
 */
static double score_diplotype_backcross(unsigned int hap1, unsigned int hap2,
                                        int* pure_pop_hap_counts,
                                        int** popA_hap_counts, int** popB_hap_counts,
                                        unsigned int* haplist, int no_haps,
                                        int noSamplesPure, int noSamplesPopA, int noSamplesPopB,
                                        int noloci, double recomb_rate,
                                        unsigned long* positions, int chrom,
                                        unsigned int** haplist_all, int* no_haps_all)
{
  int idx1 = find_hap_index(hap1, haplist, no_haps);
  int idx2 = find_hap_index(hap2, haplist, no_haps);

  double denom = 2.0 * noSamplesPure + no_haps;
  double p1_pure = (idx1 >= 0 ? pure_pop_hap_counts[idx1] + 1.0 : 1.0) / denom;
  double p2_pure = (idx2 >= 0 ? pure_pop_hap_counts[idx2] + 1.0 : 1.0) / denom;

  /* Compute U_k1 for each haplotype */
  double U1 = compute_U_k1(hap1, noloci, recomb_rate, positions,
                           popA_hap_counts, popB_hap_counts,
                           haplist_all, no_haps_all,
                           noSamplesPopA, noSamplesPopB, chrom);
  double U2 = compute_U_k1(hap2, noloci, recomb_rate, positions,
                           popA_hap_counts, popB_hap_counts,
                           haplist_all, no_haps_all,
                           noSamplesPopA, noSamplesPopB, chrom);

  /* For backcross: one haplotype is pure, other is recombinant */
  double score1 = p1_pure * U2;  /* h1 is pure, h2 is recombinant */
  double score2 = p2_pure * U1;  /* h2 is pure, h1 is recombinant */

  return (score1 > score2) ? score1 : score2;
}

/*
 * Pruned likelihood for backcross models (Mb, Me).
 * Only sums over haplotype pairs up to cumulative probability threshold.
 */
double lik_b_e_pruned(struct likelihood_cache* cache, int indivIndex,
                      struct indiv** hybrid_indiv,
                      int** popB_hap_counts, int** popA_hap_counts,
                      unsigned int** haplist, int* no_haps,
                      int noSamplesPopB, int noSamplesPopA, int noChr,
                      enum modelType model, double threshold)
{
  unsigned int hap1, hap2;
  double logL = 0.0;
  double *probV = (double *)malloc(MAXHAPS * sizeof(double));
  struct scored_pair *pairs = (struct scored_pair *)malloc(MAXHAPS * sizeof(struct scored_pair));

  /* Determine which population is the "pure" one */
  int** pure_hap_counts = (model == MODEL_B) ? popA_hap_counts : popB_hap_counts;
  int noSamplesPure = (model == MODEL_B) ? noSamplesPopA : noSamplesPopB;

  for (int i = 0; i < noChr; i++) {
    int numDiplotypes = hybrid_indiv[i][indivIndex].numHaps / 2;
    int noloci = cache->no_loci[i];
    double recomb_rate = cache->recomb_rate;
    unsigned long* positions = cache->marker_positions[i];

    /* Score all haplotype pairs */
    double total_score = 0.0;
    for (int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k += 2) {
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      double score = score_diplotype_backcross(hap1, hap2,
                                               pure_hap_counts[i],
                                               popA_hap_counts, popB_hap_counts,
                                               haplist[i], no_haps[i],
                                               noSamplesPure, noSamplesPopA, noSamplesPopB,
                                               noloci, recomb_rate, positions, i,
                                               haplist, no_haps);
      pairs[k/2].pair_index = k;
      pairs[k/2].score = score;
      total_score += score;
    }

    /* Sort pairs by score (descending) */
    qsort(pairs, numDiplotypes, sizeof(struct scored_pair), compare_scored_pairs);

    /* Sum over pairs until we reach the threshold */
    double cumulative_score = 0.0;
    int pairsUsed = 0;

    for (int p = 0; p < numDiplotypes; p++) {
      int k = pairs[p].pair_index;
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      /* Get cached U values (uses actual k value from cache) */
      double U1 = compute_U_cached(cache, hap1, i, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB);
      double U2 = compute_U_cached(cache, hap2, i, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB);

      /* Get cached log P values for pure haplotype */
      int idx1 = find_hap_index(hap1, haplist[i], no_haps[i]);
      int idx2 = find_hap_index(hap2, haplist[i], no_haps[i]);
      double p1_pure, p2_pure;
      if (model == MODEL_B) {
        p1_pure = (idx1 >= 0) ? exp(cache->logP_A[i][idx1]) : 0.0;
        p2_pure = (idx2 >= 0) ? exp(cache->logP_A[i][idx2]) : 0.0;
      } else {
        p1_pure = (idx1 >= 0) ? exp(cache->logP_B[i][idx1]) : 0.0;
        p2_pure = (idx2 >= 0) ? exp(cache->logP_B[i][idx2]) : 0.0;
      }

      if (identity2_hap(hap1, hap2)) {
        probV[pairsUsed] = p1_pure * U2 + p2_pure * U1;
      } else {
        probV[pairsUsed] = p1_pure * U2;
      }
      pairsUsed++;

      cumulative_score += pairs[p].score;
      if (cumulative_score / total_score >= threshold) break;
    }

    double norm_factor = max_element(probV, pairsUsed);
    double probL;
    if (norm_factor > 0) {
      double log_norm = log(norm_factor);
      for (int a = 0; a < pairsUsed; a++) {
        probV[a] = exp(log(probV[a]) - cache->log2 - log_norm);
      }
      probL = log(kahanSum(probV, pairsUsed)) + cache->log2 + log_norm;
    } else {
      probL = -1e300;
    }
    logL += probL;
  }

  free(probV);
  free(pairs);
  return logL;
}

/*
 * Pruned likelihood for Mf model (F2 - both haplotypes recombinant).
 * Only sums over haplotype pairs up to cumulative probability threshold.
 */
double lik_f_pruned(struct likelihood_cache* cache, int indivIndex,
                    struct indiv** hybrid_indiv,
                    int** popB_hap_counts, int** popA_hap_counts,
                    unsigned int** haplist, int* no_haps,
                    int noSamplesPopB, int noSamplesPopA, int noChr,
                    double threshold)
{
  unsigned int hap1, hap2;
  double logL = 0.0;
  double *probV = (double *)malloc(MAXHAPS * sizeof(double));
  struct scored_pair *pairs = (struct scored_pair *)malloc(MAXHAPS * sizeof(struct scored_pair));

  for (int i = 0; i < noChr; i++) {
    int numDiplotypes = hybrid_indiv[i][indivIndex].numHaps / 2;
    int noloci = cache->no_loci[i];
    double recomb_rate = cache->recomb_rate;
    unsigned long* positions = cache->marker_positions[i];

    /* Score all haplotype pairs using U_k1 for both haplotypes */
    double total_score = 0.0;
    for (int k = 0; k < hybrid_indiv[i][indivIndex].numHaps; k += 2) {
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      double U1 = compute_U_k1(hap1, noloci, recomb_rate, positions,
                               popA_hap_counts, popB_hap_counts,
                               haplist, no_haps,
                               noSamplesPopA, noSamplesPopB, i);
      double U2 = compute_U_k1(hap2, noloci, recomb_rate, positions,
                               popA_hap_counts, popB_hap_counts,
                               haplist, no_haps,
                               noSamplesPopA, noSamplesPopB, i);

      double score = U1 * U2;
      pairs[k/2].pair_index = k;
      pairs[k/2].score = score;
      total_score += score;
    }

    /* Sort pairs by score (descending) */
    qsort(pairs, numDiplotypes, sizeof(struct scored_pair), compare_scored_pairs);

    /* Sum over pairs until we reach the threshold */
    double cumulative_score = 0.0;
    int pairsUsed = 0;

    for (int p = 0; p < numDiplotypes; p++) {
      int k = pairs[p].pair_index;
      hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
      hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];

      /* Get cached U values (uses actual k value from cache) */
      double U1 = compute_U_cached(cache, hap1, i, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB);
      double U2 = compute_U_cached(cache, hap2, i, popA_hap_counts, popB_hap_counts,
                                   haplist, no_haps, noSamplesPopA, noSamplesPopB);

      if (identity2_hap(hap1, hap2)) {
        probV[pairsUsed] = 2.0 * U1 * U2;
      } else {
        probV[pairsUsed] = U1 * U2;
      }
      pairsUsed++;

      cumulative_score += pairs[p].score;
      if (cumulative_score / total_score >= threshold) break;
    }

    double norm_factor = max_element(probV, pairsUsed);
    double probL;
    if (norm_factor > 0) {
      double log_norm = log(norm_factor);
      for (int a = 0; a < pairsUsed; a++) {
        probV[a] = exp(log(probV[a]) - cache->log2 - log_norm);
      }
      probL = log(kahanSum(probV, pairsUsed)) + cache->log2 + log_norm;
    } else {
      probL = -1e300;
    }
    logL += probL;
  }

  free(probV);
  free(pairs);
  return logL;
}

/*============================================================================
 * PHASED HYBRID LIKELIHOOD FUNCTIONS
 *
 * For phased hybrid data, we know exactly which haplotypes the individual has.
 * This allows direct likelihood computation without enumerating haplotype pairs.
 *============================================================================*/

/*
 * Direct likelihood for Ma/Md with phased hybrids.
 * Computes joint probability P(g1, g2 | pop) using Dirichlet-multinomial.
 * This matches the original lik_a_d calculation for a single compatible pair.
 */
double lik_a_d_phased(int indivIndex, struct indiv** hybrid_indiv,
                      int** pop_hap_counts, unsigned int** haplist, int* no_haps,
                      int noSamplesPop, int noChr)
{
  double logL = 0.0;
  double theta = 1.0 + 2.0 * noSamplesPop;
  double t1 = gammln(theta);
  double t4 = gammln(theta + 2.0);  /* normalization for drawing 2 haplotypes */

  for (int c = 0; c < noChr; c++) {
    unsigned int g1 = hybrid_indiv[c][indivIndex].genotype1;
    unsigned int g2 = hybrid_indiv[c][indivIndex].genotype2;
    int H = no_haps[c];
    double alpha = 1.0 / H;

    /* t3 = log(2) if heterozygous (hap1 != hap2), 0 if homozygous */
    double t3 = identity2_hap(g1, g2) * log(2.0);

    /* t5 = sum over haplotypes of [gammln(phi + count + alpha) - gammln(count + alpha)]
     * where phi = number of times this haplotype appears in {g1, g2} */
    double t5 = 0.0;
    for (int j = 0; j < H; j++) {
      int phi = identity1_hap(haplist[c][j], g1) + identity1_hap(haplist[c][j], g2);
      t5 += gammln(phi + pop_hap_counts[c][j] + alpha) -
            gammln(pop_hap_counts[c][j] + alpha);
    }

    logL += t1 + t3 - t4 + t5;
  }

  return logL;
}

/*
 * Direct likelihood for Mc (F1) with phased hybrids.
 * L = P(g1|A)*P(g2|B) + P(g1|B)*P(g2|A)
 */
double lik_c_phased(int indivIndex, struct indiv** hybrid_indiv,
                    int** popB_hap_counts, int** popA_hap_counts,
                    unsigned int** haplist, int* no_haps,
                    int noSamplesPopB, int noSamplesPopA, int noChr)
{
  double logL = 0.0;

  for (int c = 0; c < noChr; c++) {
    unsigned int g1 = hybrid_indiv[c][indivIndex].genotype1;
    unsigned int g2 = hybrid_indiv[c][indivIndex].genotype2;

    /* Compute log P for each haplotype in each population */
    double logP1_A = log_prob_single_hap(g1, popA_hap_counts, haplist, no_haps, noSamplesPopA, c);
    double logP2_A = log_prob_single_hap(g2, popA_hap_counts, haplist, no_haps, noSamplesPopA, c);
    double logP1_B = log_prob_single_hap(g1, popB_hap_counts, haplist, no_haps, noSamplesPopB, c);
    double logP2_B = log_prob_single_hap(g2, popB_hap_counts, haplist, no_haps, noSamplesPopB, c);

    /* L = P(g1|A)*P(g2|B) + P(g1|B)*P(g2|A) for het, just P(g|A)*P(g|B) for hom */
    double log_term1 = logP1_A + logP2_B;
    double logL_chr;

    if (identity2_hap(g1, g2)) {
      /* Heterozygous: sum both orderings using log-sum-exp */
      double log_term2 = logP1_B + logP2_A;
      double max_log = (log_term1 > log_term2) ? log_term1 : log_term2;
      logL_chr = max_log + log(exp(log_term1 - max_log) + exp(log_term2 - max_log));
    } else {
      /* Homozygous: just one term */
      logL_chr = log_term1;
    }

    logL += logL_chr;
  }

  return logL;
}

/*
 * Direct likelihood for Mb/Me (backcross) with phased hybrids.
 * L = P(g1|pure)*U(g2) + P(g2|pure)*U(g1)
 * For Mb: pure population is A
 * For Me: pure population is B
 */
double lik_b_e_phased(struct likelihood_cache* cache, int indivIndex,
                      struct indiv** hybrid_indiv,
                      int** popB_hap_counts, int** popA_hap_counts,
                      unsigned int** haplist, int* no_haps,
                      int noSamplesPopB, int noSamplesPopA, int noChr,
                      enum modelType model)
{
  double logL = 0.0;
  int** pure_counts = (model == MODEL_B) ? popA_hap_counts : popB_hap_counts;
  int noSamplesPure = (model == MODEL_B) ? noSamplesPopA : noSamplesPopB;

  for (int c = 0; c < noChr; c++) {
    unsigned int g1 = hybrid_indiv[c][indivIndex].genotype1;
    unsigned int g2 = hybrid_indiv[c][indivIndex].genotype2;

    /* Compute P(hap|pure) for each haplotype */
    double logP1_pure = log_prob_single_hap(g1, pure_counts, haplist, no_haps, noSamplesPure, c);
    double logP2_pure = log_prob_single_hap(g2, pure_counts, haplist, no_haps, noSamplesPure, c);

    /* Compute U(hap) for each haplotype */
    double U1 = compute_U_cached(cache, g1, c, popA_hap_counts, popB_hap_counts,
                                 haplist, no_haps, noSamplesPopA, noSamplesPopB);
    double U2 = compute_U_cached(cache, g2, c, popA_hap_counts, popB_hap_counts,
                                 haplist, no_haps, noSamplesPopA, noSamplesPopB);

    /* L = P(g1|pure)*U(g2) + P(g2|pure)*U(g1) for het, or just P*U for hom */
    double term1 = exp(logP1_pure) * U2;
    double L_chr;
    if (identity2_hap(g1, g2)) {
      /* Heterozygous: sum over both assignments */
      double term2 = exp(logP2_pure) * U1;
      L_chr = term1 + term2;
    } else {
      /* Homozygous: only one term */
      L_chr = term1;
    }
    logL += log(L_chr);
  }

  return logL;
}

/*
 * Direct likelihood for Mf (F2) with phased hybrids.
 * L = U(g1) * U(g2) * (2 if g1!=g2, 1 if g1==g2)
 */
double lik_f_phased(struct likelihood_cache* cache, int indivIndex,
                    struct indiv** hybrid_indiv,
                    int** popB_hap_counts, int** popA_hap_counts,
                    unsigned int** haplist, int* no_haps,
                    int noSamplesPopB, int noSamplesPopA, int noChr)
{
  double logL = 0.0;

  for (int c = 0; c < noChr; c++) {
    unsigned int g1 = hybrid_indiv[c][indivIndex].genotype1;
    unsigned int g2 = hybrid_indiv[c][indivIndex].genotype2;

    /* Compute U(hap) for each haplotype */
    double U1 = compute_U_cached(cache, g1, c, popA_hap_counts, popB_hap_counts,
                                 haplist, no_haps, noSamplesPopA, noSamplesPopB);
    double U2 = compute_U_cached(cache, g2, c, popA_hap_counts, popB_hap_counts,
                                 haplist, no_haps, noSamplesPopA, noSamplesPopB);

    /* L = U(g1) * U(g2) * factor */
    /* identity2_hap returns 1 if different (het), 0 if same (hom) */
    double L_chr;
    if (identity2_hap(g1, g2)) {
      L_chr = 2.0 * U1 * U2;  /* heterozygous: two orderings */
    } else {
      L_chr = U1 * U2;  /* homozygous: only one way */
    }

    logL += log(L_chr);
  }

  return logL;
}
