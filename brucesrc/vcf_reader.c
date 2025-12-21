/*
 * vcf_reader.c - VCF file reading for MONGRAIL using htslib
 *
 * Reads phased, biallelic SNP data from VCF files for three populations.
 * Enforces a limit of MAXLOCI (12) SNPs per chromosome.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "mongrail.h"

extern int verbose;

/* Check if a VCF record is a valid biallelic SNP */
static int is_valid_snp(bcf1_t *rec, bcf_hdr_t *hdr)
{
  /* Must be biallelic */
  if (rec->n_allele != 2) {
    return 0;
  }

  /* Both alleles must be single nucleotides (SNP, not indel) */
  if (strlen(rec->d.allele[0]) != 1 || strlen(rec->d.allele[1]) != 1) {
    return 0;
  }

  return 1;
}

/* Check if all non-missing genotypes are phased for a record.
 * Returns 1 if valid (all non-missing are phased), 0 if invalid (unphased found).
 * Sets missing_samples bitmask to indicate which samples have missing data.
 */
static int check_phased_with_missing(bcf1_t *rec, bcf_hdr_t *hdr, int *gt_arr, int ngt,
				     unsigned int *missing_samples)
{
  int nsamples = bcf_hdr_nsamples(hdr);
  *missing_samples = 0;

  for (int i = 0; i < nsamples; i++) {
    if (ngt < 2 * (i + 1)) return 0;

    int allele1 = gt_arr[2*i];
    int allele2 = gt_arr[2*i + 1];

    /* Check for missing genotypes - mark but don't reject */
    if (bcf_gt_is_missing(allele1) || bcf_gt_is_missing(allele2)) {
      *missing_samples |= (1U << i);
      continue;  /* Skip phasing check for missing data */
    }

    /* Check if phased - second allele should have phase bit set */
    if (!bcf_gt_is_phased(allele2)) {
      return 0;  /* Unphased non-missing genotype is an error */
    }
  }

  return 1;
}

/* Read a single VCF file and extract haplotypes.
 * missing_masks[chrom][sample] is a bitmask of missing loci for each sample.
 */
static int read_single_vcf(const char *filename, int *nsamples, int *nchrom,
			   char chr_names[MAXCHRNUM][MAXNAMESZ],
			   int no_loci[MAXCHRNUM],
			   unsigned long **positions,
			   unsigned int **haplotypes,
			   unsigned int **missing_masks,
			   int *chrom_count)
{
  htsFile *fp = hts_open(filename, "r");
  if (!fp) {
    fprintf(stderr, "mongrail: error: cannot open VCF file '%s'\n", filename);
    return -1;
  }

  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  if (!hdr) {
    fprintf(stderr, "mongrail: error: cannot read VCF header from '%s'\n", filename);
    hts_close(fp);
    return -1;
  }

  *nsamples = bcf_hdr_nsamples(hdr);
  if (*nsamples == 0) {
    fprintf(stderr, "mongrail: error: no samples in VCF file '%s'\n", filename);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    return -1;
  }

  if (*nsamples > MAXINDIV) {
    fprintf(stderr, "mongrail: error: too many samples (%d > %d) in '%s'\n",
	    *nsamples, MAXINDIV, filename);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    return -1;
  }

  bcf1_t *rec = bcf_init();
  int *gt_arr = NULL;
  int ngt_arr = 0;

  char current_chrom[MAXNAMESZ] = "";
  int current_chrom_idx = -1;
  int total_snps = 0;

  /* Initialize chromosome tracking */
  *chrom_count = 0;

  while (bcf_read(fp, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_ALL);

    /* Get chromosome name */
    const char *chrom = bcf_hdr_id2name(hdr, rec->rid);

    /* Check if we're on a new chromosome */
    if (strcmp(chrom, current_chrom) != 0) {
      /* Find or add chromosome */
      int found = -1;
      for (int i = 0; i < *chrom_count; i++) {
	if (strcmp(chr_names[i], chrom) == 0) {
	  found = i;
	  break;
	}
      }

      if (found < 0) {
	if (*chrom_count >= MAXCHRNUM) {
	  fprintf(stderr, "mongrail: error: too many chromosomes (> %d)\n", MAXCHRNUM);
	  goto error;
	}
	found = *chrom_count;
	strncpy(chr_names[found], chrom, MAXNAMESZ - 1);
	chr_names[found][MAXNAMESZ - 1] = '\0';
	no_loci[found] = 0;
	(*chrom_count)++;
      }

      current_chrom_idx = found;
      strncpy(current_chrom, chrom, MAXNAMESZ - 1);
      current_chrom[MAXNAMESZ - 1] = '\0';
    }

    /* Validate: must be biallelic SNP */
    if (!is_valid_snp(rec, hdr)) {
      fprintf(stderr, "mongrail: error: variant at %s:%ld is not a biallelic SNP\n",
	      chrom, (long)(rec->pos + 1));
      goto error;
    }

    /* Check SNP limit per chromosome */
    if (no_loci[current_chrom_idx] >= MAXLOCI) {
      fprintf(stderr, "mongrail: error: too many SNPs on chromosome %s (limit: %d)\n",
	      chrom, MAXLOCI);
      goto error;
    }

    /* Get genotypes */
    int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
    if (ngt <= 0) {
      fprintf(stderr, "mongrail: error: cannot read genotypes at %s:%ld\n",
	      chrom, (long)(rec->pos + 1));
      goto error;
    }

    /* Validate: must be phased (missing genotypes are allowed) */
    unsigned int missing_samples = 0;
    if (!check_phased_with_missing(rec, hdr, gt_arr, ngt, &missing_samples)) {
      fprintf(stderr, "mongrail: error: unphased genotype at %s:%ld\n",
	      chrom, (long)(rec->pos + 1));
      goto error;
    }

    /* Store position */
    int locus_idx = no_loci[current_chrom_idx];
    positions[current_chrom_idx][locus_idx] = rec->pos + 1;  /* 1-based */

    /* Extract haplotypes for each sample */
    for (int s = 0; s < *nsamples; s++) {
      /* Check if this sample has missing data at this locus */
      if (missing_samples & (1U << s)) {
	/* Mark this locus as missing for this sample */
	missing_masks[current_chrom_idx][s] |= (1U << locus_idx);
	/* Leave haplotype bits as 0 - will be expanded later */
	continue;
      }

      int allele1 = bcf_gt_allele(gt_arr[2*s]);
      int allele2 = bcf_gt_allele(gt_arr[2*s + 1]);

      /* Set bit in haplotype if allele is 1 */
      if (allele1 == 1) {
	haplotypes[current_chrom_idx][s] |= (1U << locus_idx);
      }
      if (allele2 == 1) {
	haplotypes[current_chrom_idx][s + *nsamples] |= (1U << locus_idx);
      }
    }

    no_loci[current_chrom_idx]++;
    total_snps++;
  }

  if (total_snps == 0) {
    fprintf(stderr, "mongrail: error: no valid SNPs found in '%s'\n", filename);
    goto error;
  }

  free(gt_arr);
  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  hts_close(fp);

  *nchrom = *chrom_count;
  return 0;

error:
  free(gt_arr);
  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  hts_close(fp);
  return -1;
}

/*
 * Read three VCF files (popA, popB, hybrids) and validate consistency.
 * All files must have the same chromosomes and SNP positions.
 * Missing genotypes in hybrids are tracked via hybrid_missing_masks.
 */
int readVCFFiles(char popAfileNm[], char popBfileNm[], char hybridfileNm[],
		 int* noSamplesPopA, int* noSamplesPopB, int* noSamplesPophybrid,
		 int* noChrom, char chr_names[MAXCHRNUM][MAXNAMESZ], int no_loci[MAXCHRNUM],
		 unsigned long*** marker_positions, unsigned int*** popA_haplotypes,
		 unsigned int*** popB_haplotypes, unsigned int*** hybrid_haplotypes,
		 unsigned int*** hybrid_missing_masks)
{
  /* Temporary storage for each file */
  char chr_names_A[MAXCHRNUM][MAXNAMESZ];
  char chr_names_B[MAXCHRNUM][MAXNAMESZ];
  char chr_names_H[MAXCHRNUM][MAXNAMESZ];

  int no_loci_A[MAXCHRNUM] = {0};
  int no_loci_B[MAXCHRNUM] = {0};
  int no_loci_H[MAXCHRNUM] = {0};

  int nchrom_A, nchrom_B, nchrom_H;

  /* Allocate temporary position arrays */
  unsigned long *positions_A[MAXCHRNUM];
  unsigned long *positions_B[MAXCHRNUM];
  unsigned long *positions_H[MAXCHRNUM];

  for (int i = 0; i < MAXCHRNUM; i++) {
    positions_A[i] = calloc(MAXLOCI, sizeof(unsigned long));
    positions_B[i] = calloc(MAXLOCI, sizeof(unsigned long));
    positions_H[i] = calloc(MAXLOCI, sizeof(unsigned long));
  }

  /* Allocate temporary haplotype arrays */
  unsigned int *haps_A[MAXCHRNUM];
  unsigned int *haps_B[MAXCHRNUM];
  unsigned int *haps_H[MAXCHRNUM];

  for (int i = 0; i < MAXCHRNUM; i++) {
    haps_A[i] = calloc(2 * MAXINDIV, sizeof(unsigned int));
    haps_B[i] = calloc(2 * MAXINDIV, sizeof(unsigned int));
    haps_H[i] = calloc(2 * MAXINDIV, sizeof(unsigned int));
  }

  /* Allocate temporary missing mask arrays (only for hybrids) */
  unsigned int *missing_A[MAXCHRNUM];
  unsigned int *missing_B[MAXCHRNUM];
  unsigned int *missing_H[MAXCHRNUM];

  for (int i = 0; i < MAXCHRNUM; i++) {
    missing_A[i] = calloc(MAXINDIV, sizeof(unsigned int));
    missing_B[i] = calloc(MAXINDIV, sizeof(unsigned int));
    missing_H[i] = calloc(MAXINDIV, sizeof(unsigned int));
  }

  /* Read population A */
  if (verbose) fprintf(stderr, "Reading VCF: %s\n", popAfileNm);
  if (read_single_vcf(popAfileNm, noSamplesPopA, &nchrom_A, chr_names_A,
		      no_loci_A, positions_A, haps_A, missing_A, &nchrom_A) < 0) {
    return -1;
  }

  /* Read population B */
  if (verbose) fprintf(stderr, "Reading VCF: %s\n", popBfileNm);
  if (read_single_vcf(popBfileNm, noSamplesPopB, &nchrom_B, chr_names_B,
		      no_loci_B, positions_B, haps_B, missing_B, &nchrom_B) < 0) {
    return -1;
  }

  /* Read hybrids */
  if (verbose) fprintf(stderr, "Reading VCF: %s\n", hybridfileNm);
  if (read_single_vcf(hybridfileNm, noSamplesPophybrid, &nchrom_H, chr_names_H,
		      no_loci_H, positions_H, haps_H, missing_H, &nchrom_H) < 0) {
    return -1;
  }

  /* Validate: same number of chromosomes */
  if (nchrom_A != nchrom_B || nchrom_A != nchrom_H) {
    fprintf(stderr, "mongrail: error: VCF files have different numbers of chromosomes "
	    "(A=%d, B=%d, hybrids=%d)\n", nchrom_A, nchrom_B, nchrom_H);
    return -1;
  }

  *noChrom = nchrom_A;

  /* Validate: same chromosomes and positions */
  for (int c = 0; c < *noChrom; c++) {
    if (strcmp(chr_names_A[c], chr_names_B[c]) != 0 ||
	strcmp(chr_names_A[c], chr_names_H[c]) != 0) {
      fprintf(stderr, "mongrail: error: chromosome mismatch at index %d: A=%s, B=%s, H=%s\n",
	      c, chr_names_A[c], chr_names_B[c], chr_names_H[c]);
      return -1;
    }

    if (no_loci_A[c] != no_loci_B[c] || no_loci_A[c] != no_loci_H[c]) {
      fprintf(stderr, "mongrail: error: different number of SNPs on %s: A=%d, B=%d, H=%d\n",
	      chr_names_A[c], no_loci_A[c], no_loci_B[c], no_loci_H[c]);
      return -1;
    }

    for (int l = 0; l < no_loci_A[c]; l++) {
      if (positions_A[c][l] != positions_B[c][l] ||
	  positions_A[c][l] != positions_H[c][l]) {
	fprintf(stderr, "mongrail: error: position mismatch on %s at locus %d\n",
		chr_names_A[c], l);
	return -1;
      }
    }
  }

  /* Copy to output arrays */
  for (int c = 0; c < *noChrom; c++) {
    strncpy(chr_names[c], chr_names_A[c], MAXNAMESZ - 1);
    chr_names[c][MAXNAMESZ - 1] = '\0';
    no_loci[c] = no_loci_A[c];
  }

  /* Allocate output arrays */
  *marker_positions = malloc(*noChrom * sizeof(unsigned long*));
  *popA_haplotypes = malloc(*noChrom * sizeof(unsigned int*));
  *popB_haplotypes = malloc(*noChrom * sizeof(unsigned int*));
  *hybrid_haplotypes = malloc(*noChrom * sizeof(unsigned int*));
  *hybrid_missing_masks = malloc(*noChrom * sizeof(unsigned int*));

  for (int c = 0; c < *noChrom; c++) {
    (*marker_positions)[c] = malloc(no_loci[c] * sizeof(unsigned long));
    (*popA_haplotypes)[c] = malloc(2 * (*noSamplesPopA) * sizeof(unsigned int));
    (*popB_haplotypes)[c] = malloc(2 * (*noSamplesPopB) * sizeof(unsigned int));
    (*hybrid_haplotypes)[c] = malloc(2 * (*noSamplesPophybrid) * sizeof(unsigned int));
    (*hybrid_missing_masks)[c] = malloc((*noSamplesPophybrid) * sizeof(unsigned int));

    /* Copy positions */
    for (int l = 0; l < no_loci[c]; l++) {
      (*marker_positions)[c][l] = positions_A[c][l];
    }

    /* Copy haplotypes */
    for (int s = 0; s < 2 * (*noSamplesPopA); s++) {
      (*popA_haplotypes)[c][s] = haps_A[c][s];
    }
    for (int s = 0; s < 2 * (*noSamplesPopB); s++) {
      (*popB_haplotypes)[c][s] = haps_B[c][s];
    }
    for (int s = 0; s < 2 * (*noSamplesPophybrid); s++) {
      (*hybrid_haplotypes)[c][s] = haps_H[c][s];
    }
    /* Copy missing masks for hybrids */
    for (int s = 0; s < *noSamplesPophybrid; s++) {
      (*hybrid_missing_masks)[c][s] = missing_H[c][s];
    }
  }

  /* Free temporary arrays */
  for (int i = 0; i < MAXCHRNUM; i++) {
    free(positions_A[i]);
    free(positions_B[i]);
    free(positions_H[i]);
    free(haps_A[i]);
    free(haps_B[i]);
    free(haps_H[i]);
    free(missing_A[i]);
    free(missing_B[i]);
    free(missing_H[i]);
  }

  /* Print summary if verbose */
  if (verbose) {
    fprintf(stderr, "\nVCF Summary:\n");
    fprintf(stderr, "  Chromosomes: %d\n", *noChrom);
    fprintf(stderr, "  Pop A samples: %d\n", *noSamplesPopA);
    fprintf(stderr, "  Pop B samples: %d\n", *noSamplesPopB);
    fprintf(stderr, "  Hybrid samples: %d\n", *noSamplesPophybrid);
    fprintf(stderr, "  SNPs per chromosome:\n");
    for (int c = 0; c < *noChrom; c++) {
      fprintf(stderr, "    %s: %d\n", chr_names[c], no_loci[c]);
    }
    fprintf(stderr, "\n");
  }

  return 0;
}
