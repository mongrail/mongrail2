/*
 * vcf_reader.c - VCF file reading for MONGRAIL using htslib
 *
 * Reads phased, biallelic SNP data from VCF files for three populations.
 * Supports files with different loci - missing loci are treated as missing data.
 * Enforces a limit of MAXLOCI SNPs per chromosome.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "mongrail.h"

extern int verbose;

/*
 * Structure to hold a locus (chromosome + position)
 */
struct vcf_locus {
  char chr[MAXNAMESZ];
  unsigned long pos;
};

/*
 * Compare two loci for sorting
 */
static int compare_vcf_loci(const void *a, const void *b)
{
  const struct vcf_locus *la = (const struct vcf_locus *)a;
  const struct vcf_locus *lb = (const struct vcf_locus *)b;
  int chr_cmp = strcmp(la->chr, lb->chr);
  if (chr_cmp != 0) return chr_cmp;
  if (la->pos < lb->pos) return -1;
  if (la->pos > lb->pos) return 1;
  return 0;
}

/*
 * First pass: collect loci (chr:pos) from a VCF file.
 * Returns number of loci and populates the loci array.
 */
static int collect_vcf_loci(const char *filename, struct vcf_locus *loci, int max_loci)
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

  bcf1_t *rec = bcf_init();
  int count = 0;

  while (bcf_read(fp, hdr, rec) == 0) {
    if (count >= max_loci) {
      fprintf(stderr, "mongrail: error: too many loci in '%s'\n", filename);
      bcf_destroy(rec);
      bcf_hdr_destroy(hdr);
      hts_close(fp);
      return -1;
    }

    const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
    strncpy(loci[count].chr, chrom, MAXNAMESZ - 1);
    loci[count].chr[MAXNAMESZ - 1] = '\0';
    loci[count].pos = rec->pos + 1;  /* 1-based */
    count++;
  }

  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  hts_close(fp);

  return count;
}

/*
 * Merge loci from three files into a sorted, deduplicated list.
 */
static int merge_vcf_loci(struct vcf_locus *a, int n_a,
			  struct vcf_locus *b, int n_b,
			  struct vcf_locus *h, int n_h,
			  struct vcf_locus *merged)
{
  int total = 0;
  for (int i = 0; i < n_a; i++) merged[total++] = a[i];
  for (int i = 0; i < n_b; i++) merged[total++] = b[i];
  for (int i = 0; i < n_h; i++) merged[total++] = h[i];

  qsort(merged, total, sizeof(struct vcf_locus), compare_vcf_loci);

  if (total == 0) return 0;
  int unique = 1;
  for (int i = 1; i < total; i++) {
    if (strcmp(merged[i].chr, merged[unique-1].chr) != 0 ||
	merged[i].pos != merged[unique-1].pos) {
      merged[unique++] = merged[i];
    }
  }

  return unique;
}

/*
 * Build chromosome list from merged loci
 */
static int build_chrom_list_from_vcf_loci(struct vcf_locus *merged, int n_merged,
					  char chr_names[MAXCHRNUM][MAXNAMESZ],
					  int no_loci[MAXCHRNUM])
{
  if (n_merged == 0) return 0;

  int n_chrom = 0;
  char current_chr[MAXNAMESZ] = "";

  for (int i = 0; i < n_merged; i++) {
    if (strcmp(merged[i].chr, current_chr) != 0) {
      if (n_chrom >= MAXCHRNUM) {
	fprintf(stderr, "mongrail: error: too many chromosomes (limit: %d)\n", MAXCHRNUM);
	exit(1);
      }
      strncpy(chr_names[n_chrom], merged[i].chr, MAXNAMESZ - 1);
      chr_names[n_chrom][MAXNAMESZ - 1] = '\0';
      strncpy(current_chr, merged[i].chr, MAXNAMESZ - 1);
      current_chr[MAXNAMESZ - 1] = '\0';
      no_loci[n_chrom] = 1;
      n_chrom++;
    } else {
      no_loci[n_chrom - 1]++;
      if (no_loci[n_chrom - 1] > MAXLOCI) {
	fprintf(stderr, "mongrail: error: too many SNPs (%d) on chromosome '%s' (limit: %d)\n",
		no_loci[n_chrom - 1], current_chr, MAXLOCI);
	exit(1);
      }
    }
  }

  return n_chrom;
}

/*
 * Find merged locus index for a given chr:pos
 */
static int find_merged_locus(struct vcf_locus *merged, int n_merged,
			     const char *chr, unsigned long pos)
{
  int lo = 0, hi = n_merged - 1;
  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    int cmp = strcmp(merged[mid].chr, chr);
    if (cmp == 0) {
      if (merged[mid].pos == pos) return mid;
      if (merged[mid].pos < pos) lo = mid + 1;
      else hi = mid - 1;
    } else if (cmp < 0) {
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }
  return -1;
}

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
 * Read a VCF file and fill haplotypes at merged locus positions.
 * Loci not in this file will be marked as missing.
 */
static int read_vcf_with_merged_loci(const char *filename,
				     struct vcf_locus *merged, int n_merged,
				     int *nsamples,
				     unsigned int **haplotypes,
				     unsigned int **missing_masks)
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

  /* First, mark all loci as missing for all samples */
  for (int i = 0; i < n_merged; i++) {
    for (int s = 0; s < *nsamples; s++) {
      /* Find which chromosome this merged locus belongs to */
      int chrom_idx = 0;
      int locus_in_chrom = i;
      char current_chr[MAXNAMESZ] = "";
      for (int j = 0; j <= i; j++) {
	if (strcmp(merged[j].chr, current_chr) != 0) {
	  if (j > 0) chrom_idx++;
	  strncpy(current_chr, merged[j].chr, MAXNAMESZ - 1);
	  locus_in_chrom = 0;
	} else if (j < i) {
	  locus_in_chrom++;
	}
      }
      missing_masks[chrom_idx][s] |= (1U << locus_in_chrom);
    }
  }

  /* Read VCF and update haplotypes for loci that exist */
  while (bcf_read(fp, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_ALL);

    const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
    unsigned long pos = rec->pos + 1;

    /* Find this locus in merged list */
    int merged_idx = find_merged_locus(merged, n_merged, chrom, pos);
    if (merged_idx < 0) {
      /* This shouldn't happen if collect_vcf_loci worked correctly */
      continue;
    }

    /* Find chromosome index and locus index within chromosome */
    int chrom_idx = 0;
    int locus_in_chrom = 0;
    char current_chr[MAXNAMESZ] = "";
    for (int j = 0; j <= merged_idx; j++) {
      if (strcmp(merged[j].chr, current_chr) != 0) {
	if (j > 0) chrom_idx++;
	strncpy(current_chr, merged[j].chr, MAXNAMESZ - 1);
	locus_in_chrom = 0;
      } else {
	locus_in_chrom++;
      }
    }

    /* Validate SNP */
    if (!is_valid_snp(rec, hdr)) {
      fprintf(stderr, "mongrail: error: variant at %s:%ld is not a biallelic SNP\n",
	      chrom, (long)pos);
      goto error;
    }

    /* Get genotypes */
    int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
    if (ngt <= 0) {
      fprintf(stderr, "mongrail: error: cannot read genotypes at %s:%ld\n",
	      chrom, (long)pos);
      goto error;
    }

    /* Process each sample */
    for (int s = 0; s < *nsamples; s++) {
      if (ngt < 2 * (s + 1)) continue;

      int allele1 = gt_arr[2*s];
      int allele2 = gt_arr[2*s + 1];

      /* Check for missing genotypes */
      if (bcf_gt_is_missing(allele1) || bcf_gt_is_missing(allele2)) {
	/* Leave as missing (already marked above) */
	continue;
      }

      /* Check phasing */
      if (!bcf_gt_is_phased(allele2)) {
	fprintf(stderr, "mongrail: error: unphased genotype at %s:%ld\n",
		chrom, (long)pos);
	goto error;
      }

      /* Clear missing bit for this locus/sample */
      missing_masks[chrom_idx][s] &= ~(1U << locus_in_chrom);

      /* Set haplotype bits */
      if (bcf_gt_allele(allele1) == 1) {
	haplotypes[chrom_idx][s] |= (1U << locus_in_chrom);
      }
      if (bcf_gt_allele(allele2) == 1) {
	haplotypes[chrom_idx][s + *nsamples] |= (1U << locus_in_chrom);
      }
    }
  }

  free(gt_arr);
  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  hts_close(fp);
  return 0;

error:
  free(gt_arr);
  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  hts_close(fp);
  return -1;
}

/*
 * Read three VCF files (popA, popB, hybrids).
 * Files may have different loci - missing loci are treated as missing data.
 */
int readVCFFiles(char popAfileNm[], char popBfileNm[], char hybridfileNm[],
		 int* noSamplesPopA, int* noSamplesPopB, int* noSamplesPophybrid,
		 int* noChrom, char chr_names[MAXCHRNUM][MAXNAMESZ], int no_loci[MAXCHRNUM],
		 unsigned long*** marker_positions, unsigned int*** popA_haplotypes,
		 unsigned int*** popB_haplotypes, unsigned int*** hybrid_haplotypes,
		 unsigned int*** hybrid_missing_masks)
{
  /* First pass: collect loci from all three files */
  int max_loci = MAXCHRNUM * MAXLOCI;
  struct vcf_locus *loci_A = malloc(max_loci * sizeof(struct vcf_locus));
  struct vcf_locus *loci_B = malloc(max_loci * sizeof(struct vcf_locus));
  struct vcf_locus *loci_H = malloc(max_loci * sizeof(struct vcf_locus));

  if (verbose) fprintf(stderr, "Collecting loci from VCF files...\n");

  int n_loci_A = collect_vcf_loci(popAfileNm, loci_A, max_loci);
  if (n_loci_A < 0) return -1;

  int n_loci_B = collect_vcf_loci(popBfileNm, loci_B, max_loci);
  if (n_loci_B < 0) return -1;

  int n_loci_H = collect_vcf_loci(hybridfileNm, loci_H, max_loci);
  if (n_loci_H < 0) return -1;

  /* Merge loci from all three files */
  struct vcf_locus *merged = malloc((n_loci_A + n_loci_B + n_loci_H) * sizeof(struct vcf_locus));
  int n_merged = merge_vcf_loci(loci_A, n_loci_A, loci_B, n_loci_B,
				loci_H, n_loci_H, merged);

  if (n_merged == 0) {
    fprintf(stderr, "mongrail: error: no loci found in VCF files\n");
    return -1;
  }

  /* Report if files have different loci */
  if (verbose && (n_loci_A != n_merged || n_loci_B != n_merged || n_loci_H != n_merged)) {
    fprintf(stderr, "Note: VCF files have different loci (A=%d, B=%d, H=%d, merged=%d)\n",
	    n_loci_A, n_loci_B, n_loci_H, n_merged);
    fprintf(stderr, "      Missing loci will be treated as missing data\n");
  }

  /* Build chromosome list from merged loci */
  memset(no_loci, 0, MAXCHRNUM * sizeof(int));
  *noChrom = build_chrom_list_from_vcf_loci(merged, n_merged, chr_names, no_loci);

  /* Allocate marker positions */
  *marker_positions = malloc(*noChrom * sizeof(unsigned long*));
  int locus_idx = 0;
  for (int c = 0; c < *noChrom; c++) {
    (*marker_positions)[c] = malloc(no_loci[c] * sizeof(unsigned long));
    for (int l = 0; l < no_loci[c]; l++) {
      (*marker_positions)[c][l] = merged[locus_idx++].pos;
    }
  }

  /* Allocate haplotype and missing mask arrays */
  *popA_haplotypes = malloc(*noChrom * sizeof(unsigned int*));
  *popB_haplotypes = malloc(*noChrom * sizeof(unsigned int*));
  *hybrid_haplotypes = malloc(*noChrom * sizeof(unsigned int*));
  *hybrid_missing_masks = malloc(*noChrom * sizeof(unsigned int*));

  unsigned int **popA_missing = malloc(*noChrom * sizeof(unsigned int*));
  unsigned int **popB_missing = malloc(*noChrom * sizeof(unsigned int*));

  for (int c = 0; c < *noChrom; c++) {
    (*popA_haplotypes)[c] = calloc(2 * MAXINDIV, sizeof(unsigned int));
    (*popB_haplotypes)[c] = calloc(2 * MAXINDIV, sizeof(unsigned int));
    (*hybrid_haplotypes)[c] = calloc(2 * MAXINDIV, sizeof(unsigned int));
    (*hybrid_missing_masks)[c] = calloc(MAXINDIV, sizeof(unsigned int));
    popA_missing[c] = calloc(MAXINDIV, sizeof(unsigned int));
    popB_missing[c] = calloc(MAXINDIV, sizeof(unsigned int));
  }

  /* Second pass: read genotypes from each file */
  if (verbose) fprintf(stderr, "Reading genotypes from VCF files...\n");

  if (read_vcf_with_merged_loci(popAfileNm, merged, n_merged,
				noSamplesPopA, *popA_haplotypes, popA_missing) < 0) {
    return -1;
  }

  if (read_vcf_with_merged_loci(popBfileNm, merged, n_merged,
				noSamplesPopB, *popB_haplotypes, popB_missing) < 0) {
    return -1;
  }

  if (read_vcf_with_merged_loci(hybridfileNm, merged, n_merged,
				noSamplesPophybrid, *hybrid_haplotypes, *hybrid_missing_masks) < 0) {
    return -1;
  }

  /* Verbose output */
  if (verbose) {
    fprintf(stderr, "\nVCF Summary:\n");
    fprintf(stderr, "  Pop A:    %s (%d loci, %d samples)\n", popAfileNm, n_loci_A, *noSamplesPopA);
    fprintf(stderr, "  Pop B:    %s (%d loci, %d samples)\n", popBfileNm, n_loci_B, *noSamplesPopB);
    fprintf(stderr, "  Hybrids:  %s (%d loci, %d samples)\n", hybridfileNm, n_loci_H, *noSamplesPophybrid);
    fprintf(stderr, "  Merged:   %d unique loci across %d chromosomes\n\n", n_merged, *noChrom);
  }

  /* Free temporary arrays */
  free(loci_A);
  free(loci_B);
  free(loci_H);
  free(merged);
  for (int c = 0; c < *noChrom; c++) {
    free(popA_missing[c]);
    free(popB_missing[c]);
  }
  free(popA_missing);
  free(popB_missing);

  return 0;
}
