#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "mongrail.h"

int verbose = 0;

/*
 * Compute posterior probabilities from log-likelihoods using log-sum-exp trick.
 * Assumes uniform priors over the 6 genealogical classes.
 *
 * Input: logL array of 6 log-likelihoods
 * Output: post array of 6 posterior probabilities (normalized to sum to 1)
 */
static void compute_posteriors(const double *logL, double *post, int n)
{
  /* Find maximum log-likelihood for numerical stability */
  double max_ll = logL[0];
  for (int i = 1; i < n; i++) {
    if (logL[i] > max_ll) max_ll = logL[i];
  }

  /* Compute sum of exp(logL[i] - max_ll) */
  double sum_exp = 0.0;
  for (int i = 0; i < n; i++) {
    sum_exp += exp(logL[i] - max_ll);
  }

  /* Compute posteriors: exp(logL[i] - max_ll) / sum_exp */
  for (int i = 0; i < n; i++) {
    post[i] = exp(logL[i] - max_ll) / sum_exp;
  }
}

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s [options] popA popB hybrids\n", progname);
  fprintf(stderr, "\nClassify hybrid individuals using Bayesian genealogical models.\n");
  fprintf(stderr, "\nArguments:\n");
  fprintf(stderr, "  popA      Genotype file for reference population A\n");
  fprintf(stderr, "  popB      Genotype file for reference population B\n");
  fprintf(stderr, "  hybrids   Genotype file for putative hybrid individuals\n");
  fprintf(stderr, "\nOptions:\n");
  fprintf(stderr, "  -c        Read VCF format (requires phased, biallelic SNPs)\n");
  fprintf(stderr, "  -p THRESH Prune haplotype pairs by posterior probability (0-1, e.g., 0.99)\n");
  fprintf(stderr, "  -r RATE   Recombination rate per bp (default: 1e-8)\n");
  fprintf(stderr, "  -k K      Max recombinations for sparse enumeration (required if >16 loci)\n");
  fprintf(stderr, "  -v        Verbose output (show input file summary)\n");
  fprintf(stderr, "  -h        Show this help message\n");
  fprintf(stderr, "  -V        Show version information\n");
  fprintf(stderr, "\nModels:\n");
  fprintf(stderr, "  Ma  Pure population A      Md  Pure population B\n");
  fprintf(stderr, "  Mb  Backcross to A         Me  Backcross to B\n");
  fprintf(stderr, "  Mc  F1 hybrid              Mf  F2 hybrid\n");
}

static void version(void)
{
  fprintf(stderr, "mongrail %s\n", VERSION);
}

int main(int argc, char *argv[])
{
  char popAfileNm[MAXFILENMSZ];
  char popBfileNm[MAXFILENMSZ];
  char hybridfileNm[MAXFILENMSZ];

  int noChrom = 0;
  int noSamplesPopA = 0;
  int noSamplesPopB = 0;
  int noSamplesPophybrid = 0;
  unsigned int** popA_haplotypes = NULL;
  unsigned int** popB_haplotypes = NULL;
  unsigned int** hybrid_haplotypes = NULL;
  unsigned int** hybrid_missing_masks = NULL;  /* tracks missing loci per hybrid */
  int** popA_hap_counts;
  int** popB_hap_counts;
  unsigned int** haplist;
  int* nohaps;
  struct indiv** hybrid_indiv = NULL;
  unsigned long** marker_positions = NULL;
  char chr_names[MAXCHRNUM][MAXNAMESZ];
  int no_loci[MAXCHRNUM] = {0};
  double recomb_rate_per_bp = 1e-8;
  int max_recomb = -1;  /* -1 means not specified (use exhaustive) */
  int vcf_mode = 0;
  int hybrid_phased = 0;  /* set by data loading if hybrids have phased genotypes */
  double prune_threshold = -1.0;  /* -1 means no pruning */

  int opt;
  while ((opt = getopt(argc, argv, "cp:r:k:vhV")) != -1) {
    switch (opt) {
    case 'c':
      vcf_mode = 1;
      break;
    case 'p':
      prune_threshold = atof(optarg);
      if (prune_threshold <= 0 || prune_threshold > 1.0) {
	fprintf(stderr, "mongrail: error: prune threshold must be between 0 and 1\n");
	exit(1);
      }
      break;
    case 'r':
      recomb_rate_per_bp = atof(optarg);
      if (recomb_rate_per_bp <= 0) {
	fprintf(stderr, "mongrail: error: recombination rate must be positive\n");
	exit(1);
      }
      break;
    case 'k':
      max_recomb = atoi(optarg);
      if (max_recomb < 1 || max_recomb > 31) {
	fprintf(stderr, "mongrail: error: max recombinations must be between 1 and 31\n");
	exit(1);
      }
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
      usage(argv[0]);
      exit(0);
    case 'V':
      version();
      exit(0);
    default:
      usage(argv[0]);
      exit(1);
    }
  }

  if (argc - optind != 3) {
    usage(argv[0]);
    exit(1);
  }

  strncpy(popAfileNm, argv[optind], MAXFILENMSZ-1);
  popAfileNm[MAXFILENMSZ-1] = '\0';
  strncpy(popBfileNm, argv[optind+1], MAXFILENMSZ-1);
  popBfileNm[MAXFILENMSZ-1] = '\0';
  strncpy(hybridfileNm, argv[optind+2], MAXFILENMSZ-1);
  hybridfileNm[MAXFILENMSZ-1] = '\0';

  /* read data files */
  if (vcf_mode) {
    if (readVCFFiles(popAfileNm, popBfileNm, hybridfileNm,
		     &noSamplesPopA, &noSamplesPopB, &noSamplesPophybrid,
		     &noChrom, chr_names, no_loci, &marker_positions,
		     &popA_haplotypes, &popB_haplotypes, &hybrid_haplotypes,
		     &hybrid_missing_masks) < 0) {
      exit(1);
    }
    hybrid_phased = 1;  /* VCF mode requires phased data */
  } else {
    readDataFiles(popAfileNm, popBfileNm, hybridfileNm,
		  &noSamplesPopA, &noSamplesPopB, &noSamplesPophybrid,
		  &noChrom, chr_names, no_loci, &marker_positions,
		  &popA_haplotypes, &popB_haplotypes, &hybrid_haplotypes,
		  &hybrid_missing_masks, &hybrid_phased);
  }

  /* Warn if pruning specified with phased data (pruning is irrelevant) */
  if (hybrid_phased && prune_threshold > 0) {
    fprintf(stderr, "mongrail: warning: -p option ignored for phased data "
	    "(direct computation used)\n");
  }

  /* Check loci counts and validate -k option */
  int max_loci = 0;
  for (int i = 0; i < noChrom; i++) {
    if (no_loci[i] > max_loci) max_loci = no_loci[i];
  }

  if (max_loci > 32) {
    fprintf(stderr, "mongrail: error: chromosome has %d loci (max 32 supported)\n", max_loci);
    exit(1);
  }

  if (max_recomb < 0) {
    /* -k not specified: use exhaustive enumeration */
    if (max_loci > 16) {
      fprintf(stderr, "mongrail: error: chromosome has %d loci (max 16 for exhaustive enumeration)\n", max_loci);
      fprintf(stderr, "  Use -k option to enable sparse enumeration for larger datasets\n");
      exit(1);
    }
    max_recomb = 32;  /* effectively exhaustive for <=16 loci */
  }

  /* create hybrid individual data structures */
  hybrid_indiv = (struct indiv **)malloc(noChrom * sizeof(struct indiv *));
  for (int i = 0; i < noChrom; i++) {
    hybrid_indiv[i] = (struct indiv *)malloc(noSamplesPophybrid * sizeof(struct indiv));
    for (int j = 0; j < noSamplesPophybrid; j++)
      hybrid_indiv[i][j].compHaps = calloc(MAXHAPS, sizeof(unsigned int));
  }

  for (int i = 0; i < noChrom; i++) {
    for (int j = 0; j < noSamplesPophybrid; j++) {
      hybrid_indiv[i][j].genotype1 = hybrid_haplotypes[i][j];
      hybrid_indiv[i][j].genotype2 = hybrid_haplotypes[i][j + noSamplesPophybrid];
      hybrid_indiv[i][j].missing_mask = hybrid_missing_masks[i][j];

      /* If there's missing data, expand over all possible haplotypes */
      if (hybrid_indiv[i][j].missing_mask != 0) {
	expand_missing_haplotypes(&hybrid_indiv[i][j], no_loci[i]);
      } else if (hybrid_phased) {
	/* Phased data - haplotypes are known exactly */
	hybrid_indiv[i][j].numHaps = 2;
	hybrid_indiv[i][j].compHaps[0] = hybrid_indiv[i][j].genotype1;
	hybrid_indiv[i][j].compHaps[1] = hybrid_indiv[i][j].genotype2;
      } else {
	/* Unphased data - enumerate compatible haplotypes */
	hybrid_indiv[i][j].numHaps = count_haplotypes(hybrid_indiv[i][j].genotype1,
						      hybrid_indiv[i][j].genotype2,
						      no_loci[i]);
	if (hybrid_indiv[i][j].numHaps > MAXHAPS) {
	  fprintf(stderr, "mongrail: error: individual %d has %d compatible haplotypes "
		  "(exceeds MAXHAPS=%d)\n", j, hybrid_indiv[i][j].numHaps, MAXHAPS);
	  fprintf(stderr, "  This happens with too many heterozygous loci. "
		  "Consider using phased data (-c with VCF).\n");
	  exit(1);
	}
	compatible_haps(hybrid_indiv[i][j].compHaps,
			hybrid_indiv[i][j].genotype1,
			hybrid_indiv[i][j].genotype2);
	sortDiplotypes(&hybrid_indiv[i][j]);
      }
    }
  }

  /* build haplotype list */
  haplist = (unsigned int **)malloc(noChrom * sizeof(unsigned int *));
  for (int i = 0; i < noChrom; i++) {
    haplist[i] = (unsigned int *)malloc(MAXHAPS * sizeof(unsigned int));
    if (haplist[i] == NULL) {
      perror("mongrail");
      return EXIT_FAILURE;
    }
  }
  nohaps = (int *)calloc(noChrom, sizeof(int));

  for (int z = 0; z < noChrom; z++) {
    for (int j = 0; j < 2 * noSamplesPopA; j++)
      add_hap(popA_haplotypes[z][j], haplist, nohaps, z);
    for (int j = 0; j < 2 * noSamplesPopB; j++)
      add_hap(popB_haplotypes[z][j], haplist, nohaps, z);
  }

  for (int z = 0; z < noChrom; z++) {
    for (int j = 0; j < noSamplesPophybrid; j++) {
      for (int h = 0; h < hybrid_indiv[z][j].numHaps; h += 2) {
	add_hap(hybrid_indiv[z][j].compHaps[h], haplist, nohaps, z);
	add_hap(hybrid_indiv[z][j].compHaps[h+1], haplist, nohaps, z);
      }
    }
  }

  /* allocate and initialize haplotype counts (indexed by haplist position) */
  popA_hap_counts = (int **)malloc(noChrom * sizeof(int *));
  popB_hap_counts = (int **)malloc(noChrom * sizeof(int *));
  for (int i = 0; i < noChrom; i++) {
    popA_hap_counts[i] = (int *)calloc(MAXHAPS, sizeof(int));
    popB_hap_counts[i] = (int *)calloc(MAXHAPS, sizeof(int));
  }

  get_hap_counts(popA_haplotypes, popA_hap_counts, haplist, nohaps, noChrom, noSamplesPopA);
  get_hap_counts(popB_haplotypes, popB_hap_counts, haplist, nohaps, noChrom, noSamplesPopB);

  /* initialize cache for fast likelihood computation */
  struct likelihood_cache *cache = cache_init(noChrom, no_loci, recomb_rate_per_bp,
					      marker_positions,
					      popA_hap_counts, popB_hap_counts,
					      haplist, nohaps,
					      noSamplesPopA, noSamplesPopB,
					      max_recomb);

  /* output header */
  printf("# mongrail %s\n", VERSION);
  printf("# recomb_rate=%.2e  chroms=%d  n_A=%d  n_B=%d  n_hyb=%d\n",
	 recomb_rate_per_bp, noChrom, noSamplesPopA, noSamplesPopB, noSamplesPophybrid);
  printf("#\n");
  printf("# Log-likelihoods followed by posterior probabilities (assuming uniform priors)\n");
  printf("#\n");
  printf("# indiv      Ma          Mb          Mc          Md          Me          Mf       P(Ma)   P(Mb)   P(Mc)   P(Md)   P(Me)   P(Mf)\n");

  /* compute and output log-likelihoods and posterior probabilities */
  for (int i = 0; i < noSamplesPophybrid; i++) {
    double lik_d, lik_cc, lik_a;
    double lik_b, lik_e, lik_ff;

    if (hybrid_phased) {
      /* Phased data: use direct likelihood computation (no pair enumeration) */
      lik_a = lik_a_d_phased(i, hybrid_indiv, popB_hap_counts,
			     haplist, nohaps, noSamplesPopB, noChrom);
      lik_d = lik_a_d_phased(i, hybrid_indiv, popA_hap_counts,
			     haplist, nohaps, noSamplesPopA, noChrom);
      lik_cc = lik_c_phased(i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			    haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom);
      lik_b = lik_b_e_phased(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			     haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom, MODEL_B);
      lik_e = lik_b_e_phased(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			     haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom, MODEL_E);
      lik_ff = lik_f_phased(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			    haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom);
    } else if (prune_threshold > 0) {
      /* Unphased with pruning threshold: use pre-ranked functions */
      lik_d = lik_a_d_preranked(i, hybrid_indiv, popA_hap_counts,
			       haplist, nohaps,
			       noSamplesPopA, noSamplesPopB, noChrom,
			       MODEL_D, prune_threshold);
      lik_cc = lik_c_preranked(i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			      haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom,
			      prune_threshold);
      lik_a = lik_a_d_preranked(i, hybrid_indiv, popB_hap_counts,
			       haplist, nohaps,
			       noSamplesPopA, noSamplesPopB, noChrom,
			       MODEL_A, prune_threshold);
      lik_b = lik_b_e_preranked(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			       haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom,
			       MODEL_B, prune_threshold);
      lik_e = lik_b_e_preranked(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			       haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom,
			       MODEL_E, prune_threshold);
      lik_ff = lik_f_preranked(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			      haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom,
			      prune_threshold);
    } else {
      /* Unphased without pruning: use original functions */
      lik_d = lik_a_d(i, hybrid_indiv, popA_hap_counts, haplist, nohaps,
		      noSamplesPopA, noSamplesPopB, noChrom, MODEL_D);
      lik_cc = lik_c(i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
		     haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom);
      lik_a = lik_a_d(i, hybrid_indiv, popB_hap_counts, haplist, nohaps,
		      noSamplesPopA, noSamplesPopB, noChrom, MODEL_A);
      lik_b = lik_b_e_cached(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			     haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom, MODEL_B);
      lik_e = lik_b_e_cached(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			     haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom, MODEL_E);
      lik_ff = lik_f_cached(cache, i, hybrid_indiv, popB_hap_counts, popA_hap_counts,
			    haplist, nohaps, noSamplesPopB, noSamplesPopA, noChrom);
    }

    /* Compute posterior probabilities using log-sum-exp trick */
    double logL[6] = {lik_a, lik_b, lik_cc, lik_d, lik_e, lik_ff};
    double post[6];
    compute_posteriors(logL, post, 6);

    printf("%7d  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	   i, lik_a, lik_b, lik_cc, lik_d, lik_e, lik_ff,
	   post[0], post[1], post[2], post[3], post[4], post[5]);
  }

  cache_free(cache);
  return 0;
}
