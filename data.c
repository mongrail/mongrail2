#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mongrail.h"

extern int verbose;

/* Check if file appears to be VCF format */
static int is_vcf_file(const char *filename)
{
  FILE *fp = fopen(filename, "r");
  if (!fp) return 0;

  char line[256];
  if (fgets(line, sizeof(line), fp) != NULL) {
    fclose(fp);
    /* VCF files start with ## or #CHROM */
    if (strncmp(line, "##", 2) == 0 || strncmp(line, "#CHROM", 6) == 0)
      return 1;
  } else {
    fclose(fp);
  }
  return 0;
}

/*
 * Structure to hold a single locus (chromosome + position)
 */
struct locus {
  char chr[MAXNAMESZ];
  unsigned long pos;
};

/*
 * Compare two loci for sorting: first by chromosome name, then by position
 */
static int compare_loci(const void *a, const void *b)
{
  const struct locus *la = (const struct locus *)a;
  const struct locus *lb = (const struct locus *)b;
  int chr_cmp = strcmp(la->chr, lb->chr);
  if (chr_cmp != 0) return chr_cmp;
  if (la->pos < lb->pos) return -1;
  if (la->pos > lb->pos) return 1;
  return 0;
}

/*
 * Extract loci (chr:pos) from raw file data into locus array.
 * Returns number of loci extracted.
 */
static int extract_loci(char** raw_data, int nlines, struct locus* loci)
{
  char buffer[MAXLINESZ];
  const char delim[] = ":";

  for (int i = 0; i < nlines; i++) {
    strncpy(buffer, raw_data[i], MAXLINESZ);
    buffer[MAXLINESZ-1] = '\0';

    char *chr = strtok(buffer, delim);
    char *pos_str = strtok(NULL, delim);

    if (chr == NULL || pos_str == NULL) {
      fprintf(stderr, "mongrail: error: malformed line %d (expected chr:pos format)\n", i+1);
      exit(1);
    }

    strncpy(loci[i].chr, chr, MAXNAMESZ);
    loci[i].chr[MAXNAMESZ-1] = '\0';
    loci[i].pos = strtoul(pos_str, NULL, 10);
  }
  return nlines;
}

/*
 * Merge loci from three files into a sorted, deduplicated list.
 * Returns the number of unique loci.
 */
static int merge_loci(struct locus* loci_a, int n_a,
		      struct locus* loci_b, int n_b,
		      struct locus* loci_h, int n_h,
		      struct locus* merged)
{
  /* Copy all loci into merged array */
  int total = 0;
  for (int i = 0; i < n_a; i++) merged[total++] = loci_a[i];
  for (int i = 0; i < n_b; i++) merged[total++] = loci_b[i];
  for (int i = 0; i < n_h; i++) merged[total++] = loci_h[i];

  /* Sort by chromosome then position */
  qsort(merged, total, sizeof(struct locus), compare_loci);

  /* Remove duplicates */
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
 * Find the index of a locus in raw_data, or return -1 if not found.
 * Uses binary search since loci are sorted.
 */
static int find_locus_in_file(struct locus* file_loci, int n_loci,
			      const char* chr, unsigned long pos)
{
  int lo = 0, hi = n_loci - 1;
  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    int cmp = strcmp(file_loci[mid].chr, chr);
    if (cmp == 0) {
      if (file_loci[mid].pos == pos) return mid;
      if (file_loci[mid].pos < pos) lo = mid + 1;
      else hi = mid - 1;
    } else if (cmp < 0) {
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }
  return -1;  /* Not found */
}

/*
 * Set missing genotypes for all individuals at a locus.
 * Used when a locus is absent from a file.
 */
static void set_missing_genotypes(struct genotype* genos, int n_indiv)
{
  for (int i = 0; i < n_indiv; i++) {
    genos[i].g1 = MISSING_ALLELE;
    genos[i].g2 = MISSING_ALLELE;
    genos[i].phase = '|';
  }
}

/*
 * Build chromosome list from merged loci.
 * Returns number of chromosomes and populates chr_names and no_loci arrays.
 */
static int build_chrom_list_from_merged(struct locus* merged, int n_merged,
					char chr_names[MAXCHRNUM][MAXNAMESZ],
					int no_loci[MAXCHRNUM])
{
  if (n_merged == 0) return 0;

  int n_chrom = 0;
  char current_chr[MAXNAMESZ] = "";

  for (int i = 0; i < n_merged; i++) {
    if (strcmp(merged[i].chr, current_chr) != 0) {
      /* New chromosome */
      if (n_chrom >= MAXCHRNUM) {
	fprintf(stderr, "mongrail: error: too many chromosomes (limit: %d)\n", MAXCHRNUM);
	exit(1);
      }
      strncpy(chr_names[n_chrom], merged[i].chr, MAXNAMESZ);
      chr_names[n_chrom][MAXNAMESZ-1] = '\0';
      strncpy(current_chr, merged[i].chr, MAXNAMESZ);
      current_chr[MAXNAMESZ-1] = '\0';
      no_loci[n_chrom] = 1;
      n_chrom++;
    } else {
      no_loci[n_chrom-1]++;
      if (no_loci[n_chrom-1] > MAXLOCI) {
	fprintf(stderr, "mongrail: error: too many SNPs (%d) on chromosome '%s' (limit: %d)\n",
		no_loci[n_chrom-1], current_chr, MAXLOCI);
	exit(1);
      }
    }
  }

  return n_chrom;
}

void readDataFiles(char popAfileNm[], char popBfileNm[], char hybridfileNm[],
		   int* noSamplesPopA, int* noSamplesPopB, int* noSamplesPophybrid,
		   int* noChrom, char chr_names[MAXCHRNUM][MAXNAMESZ], int no_loci[MAXCHRNUM],
		   unsigned long*** marker_positions, unsigned int*** popA_haplotypes,
		   unsigned int*** popB_haplotypes, unsigned int*** hybrid_haplotypes,
		   unsigned int*** hybrid_missing_masks)
{
  /* Check if any input file is VCF format */
  if (is_vcf_file(popAfileNm) || is_vcf_file(popBfileNm) || is_vcf_file(hybridfileNm)) {
    fprintf(stderr, "mongrail: error: input file appears to be VCF format\n");
    fprintf(stderr, "Use the -c option to read VCF files: mongrail2 -c popA.vcf popB.vcf hybrids.vcf\n");
    exit(1);
  }

  /* Count lines in each file (may be different) */
  int n_lines_A = countfilelines(popAfileNm);
  int n_lines_B = countfilelines(popBfileNm);
  int n_lines_H = countfilelines(hybridfileNm);

  /* Read all three files */
  char** raw_data_popA = (char**)malloc(n_lines_A * sizeof(char*));
  char** raw_data_popB = (char**)malloc(n_lines_B * sizeof(char*));
  char** raw_data_hybrids = (char**)malloc(n_lines_H * sizeof(char*));

  file_to_array(raw_data_popA, popAfileNm);
  file_to_array(raw_data_popB, popBfileNm);
  file_to_array(raw_data_hybrids, hybridfileNm);

  /* Extract loci from each file */
  struct locus* loci_A = malloc(n_lines_A * sizeof(struct locus));
  struct locus* loci_B = malloc(n_lines_B * sizeof(struct locus));
  struct locus* loci_H = malloc(n_lines_H * sizeof(struct locus));

  extract_loci(raw_data_popA, n_lines_A, loci_A);
  extract_loci(raw_data_popB, n_lines_B, loci_B);
  extract_loci(raw_data_hybrids, n_lines_H, loci_H);

  /* Sort each file's loci for binary search later */
  qsort(loci_A, n_lines_A, sizeof(struct locus), compare_loci);
  qsort(loci_B, n_lines_B, sizeof(struct locus), compare_loci);
  qsort(loci_H, n_lines_H, sizeof(struct locus), compare_loci);

  /* Merge loci from all three files */
  int max_merged = n_lines_A + n_lines_B + n_lines_H;
  struct locus* merged_loci = malloc(max_merged * sizeof(struct locus));
  int n_merged = merge_loci(loci_A, n_lines_A, loci_B, n_lines_B,
			    loci_H, n_lines_H, merged_loci);

  if (n_merged == 0) {
    fprintf(stderr, "mongrail: error: no loci found in input files\n");
    exit(1);
  }

  /* Build chromosome list from merged loci */
  memset(no_loci, 0, MAXCHRNUM * sizeof(int));
  *noChrom = build_chrom_list_from_merged(merged_loci, n_merged, chr_names, no_loci);

  /* Report if files have different loci */
  if (verbose && (n_lines_A != n_merged || n_lines_B != n_merged || n_lines_H != n_merged)) {
    fprintf(stderr, "Note: Files have different loci (A=%d, B=%d, H=%d, merged=%d)\n",
	    n_lines_A, n_lines_B, n_lines_H, n_merged);
    fprintf(stderr, "      Missing loci will be treated as missing data\n");
  }

  /* Allocate marker positions from merged loci */
  *marker_positions = (unsigned long**)malloc(*noChrom * sizeof(unsigned long*));
  int locus_idx = 0;
  for (int c = 0; c < *noChrom; c++) {
    (*marker_positions)[c] = (unsigned long*)malloc(no_loci[c] * sizeof(unsigned long));
    for (int l = 0; l < no_loci[c]; l++) {
      (*marker_positions)[c][l] = merged_loci[locus_idx++].pos;
    }
  }

  /* Get sample counts from first locus in each file */
  struct genotype temp_genos[MAXINDIV];
  *noSamplesPopA = get_genotypes(raw_data_popA[0], temp_genos);
  *noSamplesPopB = get_genotypes(raw_data_popB[0], temp_genos);
  *noSamplesPophybrid = get_genotypes(raw_data_hybrids[0], temp_genos);

  if (*noSamplesPopA == 0) {
    fprintf(stderr, "mongrail: error: no individuals found in population A file\n");
    exit(1);
  }
  if (*noSamplesPopB == 0) {
    fprintf(stderr, "mongrail: error: no individuals found in population B file\n");
    exit(1);
  }
  if (*noSamplesPophybrid == 0) {
    fprintf(stderr, "mongrail: error: no individuals found in hybrids file\n");
    exit(1);
  }

  /* Allocate genotype data structures */
  struct genotype*** popA_genotypes = (struct genotype***)malloc(*noChrom * sizeof(struct genotype**));
  struct genotype*** popB_genotypes = (struct genotype***)malloc(*noChrom * sizeof(struct genotype**));
  struct genotype*** pophybrid_genotypes = (struct genotype***)malloc(*noChrom * sizeof(struct genotype**));

  for (int c = 0; c < *noChrom; c++) {
    popA_genotypes[c] = (struct genotype**)malloc(no_loci[c] * sizeof(struct genotype*));
    popB_genotypes[c] = (struct genotype**)malloc(no_loci[c] * sizeof(struct genotype*));
    pophybrid_genotypes[c] = (struct genotype**)malloc(no_loci[c] * sizeof(struct genotype*));
    for (int l = 0; l < no_loci[c]; l++) {
      popA_genotypes[c][l] = (struct genotype*)malloc(MAXINDIV * sizeof(struct genotype));
      popB_genotypes[c][l] = (struct genotype*)malloc(MAXINDIV * sizeof(struct genotype));
      pophybrid_genotypes[c][l] = (struct genotype*)malloc(MAXINDIV * sizeof(struct genotype));
    }
  }

  /* Process each merged locus for each population */
  locus_idx = 0;
  for (int c = 0; c < *noChrom; c++) {
    for (int l = 0; l < no_loci[c]; l++) {
      struct locus* mloc = &merged_loci[locus_idx];

      /* Population A */
      int file_idx_A = find_locus_in_file(loci_A, n_lines_A, mloc->chr, mloc->pos);
      if (file_idx_A >= 0) {
	get_genotypes(raw_data_popA[file_idx_A], popA_genotypes[c][l]);
      } else {
	set_missing_genotypes(popA_genotypes[c][l], *noSamplesPopA);
      }

      /* Population B */
      int file_idx_B = find_locus_in_file(loci_B, n_lines_B, mloc->chr, mloc->pos);
      if (file_idx_B >= 0) {
	get_genotypes(raw_data_popB[file_idx_B], popB_genotypes[c][l]);
      } else {
	set_missing_genotypes(popB_genotypes[c][l], *noSamplesPopB);
      }

      /* Hybrids */
      int file_idx_H = find_locus_in_file(loci_H, n_lines_H, mloc->chr, mloc->pos);
      if (file_idx_H >= 0) {
	get_genotypes(raw_data_hybrids[file_idx_H], pophybrid_genotypes[c][l]);
      } else {
	set_missing_genotypes(pophybrid_genotypes[c][l], *noSamplesPophybrid);
      }

      locus_idx++;
    }
  }

  /* Allocate haplotype data structures */
  *popA_haplotypes = (unsigned int**)malloc(*noChrom * sizeof(unsigned int*));
  *popB_haplotypes = (unsigned int**)malloc(*noChrom * sizeof(unsigned int*));
  *hybrid_haplotypes = (unsigned int**)malloc(*noChrom * sizeof(unsigned int*));
  *hybrid_missing_masks = (unsigned int**)malloc(*noChrom * sizeof(unsigned int*));

  /* Also need missing masks for popA and popB now */
  unsigned int** popA_missing_masks = (unsigned int**)malloc(*noChrom * sizeof(unsigned int*));
  unsigned int** popB_missing_masks = (unsigned int**)malloc(*noChrom * sizeof(unsigned int*));

  for (int c = 0; c < *noChrom; c++) {
    (*popA_haplotypes)[c] = (unsigned int*)malloc(*noSamplesPopA * 2 * sizeof(unsigned int));
    (*popB_haplotypes)[c] = (unsigned int*)malloc(*noSamplesPopB * 2 * sizeof(unsigned int));
    (*hybrid_haplotypes)[c] = (unsigned int*)malloc(*noSamplesPophybrid * 2 * sizeof(unsigned int));
    (*hybrid_missing_masks)[c] = (unsigned int*)calloc(*noSamplesPophybrid, sizeof(unsigned int));
    popA_missing_masks[c] = (unsigned int*)calloc(*noSamplesPopA, sizeof(unsigned int));
    popB_missing_masks[c] = (unsigned int*)calloc(*noSamplesPopB, sizeof(unsigned int));
  }

  /* Get haplotypes with missing data tracking for all populations */
  get_haplotypes_with_missing(*popA_haplotypes, popA_missing_masks,
			      popA_genotypes, *noChrom, *noSamplesPopA, no_loci);
  get_haplotypes_with_missing(*popB_haplotypes, popB_missing_masks,
			      popB_genotypes, *noChrom, *noSamplesPopB, no_loci);
  get_haplotypes_with_missing(*hybrid_haplotypes, *hybrid_missing_masks,
			      pophybrid_genotypes, *noChrom, *noSamplesPophybrid, no_loci);

  /* Verbose output */
  if (verbose) {
    fprintf(stderr, "Input files:\n");
    fprintf(stderr, "  Pop A:    %s (%d loci)\n", popAfileNm, n_lines_A);
    fprintf(stderr, "  Pop B:    %s (%d loci)\n", popBfileNm, n_lines_B);
    fprintf(stderr, "  Hybrids:  %s (%d loci)\n", hybridfileNm, n_lines_H);
    fprintf(stderr, "  Merged:   %d unique loci across %d chromosomes\n", n_merged, *noChrom);
    fprintf(stderr, "\nSamples: n_A=%d, n_B=%d, n_hyb=%d\n\n",
	    *noSamplesPopA, *noSamplesPopB, *noSamplesPophybrid);
  }

  /* Free temporary storage */
  for (int i = 0; i < n_lines_A; i++) free(raw_data_popA[i]);
  for (int i = 0; i < n_lines_B; i++) free(raw_data_popB[i]);
  for (int i = 0; i < n_lines_H; i++) free(raw_data_hybrids[i]);
  free(raw_data_popA);
  free(raw_data_popB);
  free(raw_data_hybrids);
  free(loci_A);
  free(loci_B);
  free(loci_H);
  free(merged_loci);

  for (int c = 0; c < *noChrom; c++) {
    for (int l = 0; l < no_loci[c]; l++) {
      free(popA_genotypes[c][l]);
      free(popB_genotypes[c][l]);
      free(pophybrid_genotypes[c][l]);
    }
    free(popA_genotypes[c]);
    free(popB_genotypes[c]);
    free(pophybrid_genotypes[c]);
    free(popA_missing_masks[c]);
    free(popB_missing_masks[c]);
  }
  free(popA_genotypes);
  free(popB_genotypes);
  free(pophybrid_genotypes);
  free(popA_missing_masks);
  free(popB_missing_masks);
} 
  
/* find if string is in array of strings */
int findstring(char strarray[MAXCHRNUM][MAXNAMESZ], char *strtarget, int len)
{
  for(int i=0; i<len; i++)
    if(strcmp(strarray[i],strtarget)==0)
      return(1);
  return(0);
}

/* count number of lines in file */
int countfilelines(char filename[])
{
  FILE *file;
  int ch;  /* fgetc returns int, not char, to properly detect EOF */
  int line_count=0;
  file = fopen(filename, "r");
  if (file == NULL) {
    fprintf(stderr, "mongrail: error: cannot open file '%s'\n", filename);
    exit(1);
  }
  while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            line_count++;
        }
    }
    fclose(file);

    if (line_count == 0) {
      fprintf(stderr, "mongrail: error: file '%s' is empty or has no newlines\n", filename);
      exit(1);
    }

    return line_count;
}

/* read file into array of strings */
void file_to_array (char *target_array[MAXLINESZ], char filename[])
{
  char buffer[MAXLINESZ];
  FILE *fp;
  int linepos = 0;
  if ((fp = fopen(filename, "r")) == NULL)
    {
      printf("Can't open %s\n", filename);
      exit(1);
    }
  while (fgets(buffer, MAXLINESZ, fp) != NULL) {
    size_t len = strlen(buffer);
    if (len > 0) {
      buffer[len - 1] = '\0';
    }
    target_array[linepos] = strdup(buffer);
    /* add memory allocation error check here */
    linepos++;
  }
  fclose(fp);
}

/* create chromosome list, count number of loci per chromosome and number of chromosomes */
int mk_chrom_list (char* file_array[MAXLINESZ], char names[MAXCHRNUM][MAXNAMESZ], int* no_loci, int nolines)
{
  const char delim[] = ":";
  char *token;
  int chrpos=0;
  char buffer[MAXLINESZ];
  for(int i=0; i<nolines; i++)
    {
      strncpy(buffer,file_array[i],MAXLINESZ);
      token = strtok(buffer,delim);
      if (token == NULL) {
	fprintf(stderr, "mongrail: error: malformed line %d (missing chromosome name)\n", i+1);
	exit(1);
      }
      if(chrpos > 0)
	{
	  if(findstring(names,token,chrpos)==0)
	    {
	      if (chrpos >= MAXCHRNUM) {
		fprintf(stderr, "mongrail: error: too many chromosomes (limit: %d)\n", MAXCHRNUM);
		exit(1);
	      }
	      strncpy(names[chrpos],token,MAXNAMESZ);
	      names[chrpos][MAXNAMESZ-1] = '\0';
	      no_loci[chrpos] = no_loci[chrpos] + 1;
	      chrpos++;
	    }
	  else
	    {
	      no_loci[chrpos-1] = no_loci[chrpos-1] + 1;
	      if (no_loci[chrpos-1] > MAXLOCI) {
		fprintf(stderr, "mongrail: error: too many SNPs (%d) on chromosome '%s' (limit: %d)\n",
			no_loci[chrpos-1], names[chrpos-1], MAXLOCI);
		exit(1);
	      }
	    }
	}
      else
	{
	  strncpy(names[chrpos],token,MAXNAMESZ);
	  names[chrpos][MAXNAMESZ-1] = '\0';
	  no_loci[chrpos] = no_loci[chrpos] + 1;;
	  chrpos++;
	}
    }
  return chrpos;
}

/* parse genotype string into genotype data structure
 * Detects missing genotypes: ./. or .|. or -/- or -|- patterns
 * Missing alleles are stored as MISSING_ALLELE (-1)
 */
int get_genotypes(char* gstring, struct genotype* chr_locus)
{
  char *token;
  const char delim[] = ": \t";
  int numindivs = 0;
  struct genotype geno[MAXINDIV];
  int g1, g2;
  char phase;

  char *str_copy = strdup(gstring);
  if (str_copy == NULL) {
    fprintf(stderr, "mongrail: error: memory allocation failed\n");
    exit(1);
  }
  token = strtok(str_copy,delim); // Tokenize based on spaces and :
  token = strtok(NULL,delim);
  token = strtok(NULL,delim);
  while (token != NULL) {
    /* Validate genotype format: must be at least 3 characters like "0|1" */
    size_t len = strlen(token);
    if (len < 3) {
      fprintf(stderr, "mongrail: error: malformed genotype '%s' (expected format like 0|1)\n", token);
      free(str_copy);
      exit(1);
    }

    /* Check bounds before adding */
    if (numindivs >= MAXINDIV) {
      fprintf(stderr, "mongrail: error: too many individuals (limit: %d)\n", MAXINDIV);
      free(str_copy);
      exit(1);
    }

    /* Check for missing data: ./. or .|. or -/- or -|- */
    if (token[0] == '.' || token[0] == '-') {
      g1 = MISSING_ALLELE;
    } else {
      g1 = token[0] - '0';
    }
    phase = token[1];
    if (token[2] == '.' || token[2] == '-') {
      g2 = MISSING_ALLELE;
    } else {
      g2 = token[2] - '0';
    }
    geno[numindivs].g1 = g1;
    geno[numindivs].g2 = g2;
    geno[numindivs].phase = phase;
    numindivs++;
    token = strtok(NULL,delim);
  }
  for(int i=0; i<numindivs; i++)
    {
      chr_locus[i].g1 = geno[i].g1;
      chr_locus[i].g2 = geno[i].g2;
      chr_locus[i].phase = geno[i].phase;
    }
  free(str_copy); // Free the copy of the string
  return numindivs;
}

/* get sampled population haplotypes as unsigned ints */
/* number of haplotypes = 2 * number of individuals sampled */
void get_haplotypes(unsigned int** haplotypes, struct genotype*** gdata, int nChr, int noInd, int* nloci)
{
  int i,j,k;
  /* Buffer size allows up to MAXLOCI loci plus null terminator */
  char hapstring_g1[MAXLOCI + 1];
  char hapstring_g2[MAXLOCI + 1];
  for(i=0; i<nChr; i++)
    {
      if(nloci[i] > MAXLOCI) {
	fprintf(stderr, "Error: number of loci (%d) exceeds maximum (%d) for chromosome %d\n", nloci[i], MAXLOCI, i);
	exit(1);
      }
      for(j=0; j<noInd; j++)
	{
	  for(k=0; k<nloci[i]; k++)
	    {
	      if(gdata[i][k][j].g1==0)
		hapstring_g1[k] = '0';
	      else
		hapstring_g1[k] = '1';
	      if(gdata[i][k][j].g2==0)
		hapstring_g2[k] = '0';
	      else
		hapstring_g2[k] = '1';
	    }
	  hapstring_g1[nloci[i]]='\0';
	  hapstring_g2[nloci[i]]='\0';
	  /* array places all g1 haplotypes first, then g2 haplotypes */
	  haplotypes[i][j]=binaryToDecimal(hapstring_g1);
	  haplotypes[i][j+noInd]=binaryToDecimal(hapstring_g2);
	}
    }
}

/* get sampled population haplotypes as unsigned ints, also tracking missing loci
 * missing_masks[chrom][indiv] is a bitmask where bit k = 1 means locus k is missing
 * For missing loci, haplotype bits are set to 0 (will be expanded later)
 */
void get_haplotypes_with_missing(unsigned int** haplotypes, unsigned int** missing_masks,
				 struct genotype*** gdata, int nChr, int noInd, int* nloci)
{
  int i,j,k;
  char hapstring_g1[MAXLOCI + 1];
  char hapstring_g2[MAXLOCI + 1];

  for(i=0; i<nChr; i++)
    {
      if(nloci[i] > MAXLOCI) {
	fprintf(stderr, "Error: number of loci (%d) exceeds maximum (%d) for chromosome %d\n", nloci[i], MAXLOCI, i);
	exit(1);
      }
      for(j=0; j<noInd; j++)
	{
	  unsigned int missing = 0;
	  for(k=0; k<nloci[i]; k++)
	    {
	      /* Check for missing data (MISSING_ALLELE = -1) */
	      if(gdata[i][k][j].g1 == MISSING_ALLELE || gdata[i][k][j].g2 == MISSING_ALLELE)
		{
		  missing |= (1U << k);
		  hapstring_g1[k] = '0';  /* placeholder, will be expanded */
		  hapstring_g2[k] = '0';
		}
	      else
		{
		  if(gdata[i][k][j].g1==0)
		    hapstring_g1[k] = '0';
		  else
		    hapstring_g1[k] = '1';
		  if(gdata[i][k][j].g2==0)
		    hapstring_g2[k] = '0';
		  else
		    hapstring_g2[k] = '1';
		}
	    }
	  hapstring_g1[nloci[i]]='\0';
	  hapstring_g2[nloci[i]]='\0';
	  haplotypes[i][j]=binaryToDecimal(hapstring_g1);
	  haplotypes[i][j+noInd]=binaryToDecimal(hapstring_g2);
	  missing_masks[i][j] = missing;
	}
    }
}

/* get marker positions for each chromosome */
void get_positions(int noChr,int* no_loci, char** raw_data, unsigned long** positions)
{
  int linepos=0;
  char* token;
  const char delim[] = ":";
  char buffer[MAXLINESZ];
  for(int i=0; i<noChr; i++)
    for(int j=0; j<no_loci[i]; j++)
      {
	if(i==0)
	  linepos=j;
	else
	  {
	    linepos = 0;
	    for(int k=0; k<i; k++)
	      linepos += no_loci[k];
	    linepos+=j;
	  }
	strncpy(buffer,raw_data[linepos],MAXLINESZ);
	token = strtok(buffer,delim);
	if (token == NULL) {
	  fprintf(stderr, "mongrail: error: malformed line %d (missing chromosome name)\n", linepos+1);
	  exit(1);
	}
	token = strtok(NULL,delim);
	if (token == NULL) {
	  fprintf(stderr, "mongrail: error: malformed line %d (missing position)\n", linepos+1);
	  exit(1);
	}
	positions[i][j] = strtol(token, NULL, 10);
      }
}

/* summarize pop samples as haplotype counts in array of length MAXHAPS with base_10 haplotype as array index */
void get_hap_counts(unsigned int** haplotypes, int** hap_counts, int noChr, int no_indiv)
{
  for(int i=0; i<noChr; i++)
    for(int j=0; j<2*no_indiv; j++)
      hap_counts[i][haplotypes[i][j]] = hap_counts[i][haplotypes[i][j]] + 1;
}

/* get number of polymorphic sites and number of haplotypes for hybrid indivs */
unsigned int count_haplotypes(unsigned int hap1, unsigned int hap2, int noloci)
{
  unsigned int hap1hap2;
  unsigned int no1bits=0;
  hap1hap2 = hap1^hap2; /* 1 if 1^0 or 0^1, 0 otherwise */
  for (int i = 0; i < 32 && i < noloci; ++i)
    if (hap1hap2 >> i & 0x1) no1bits++; /* if ith bit of hap1hap2 = 1 is polymorphic */
  if (no1bits == 0) return 2; /* if no polymorphic sites, only 2 haplotypes */
  return pow(2,no1bits); /* no haplotypes = 2^no1bits */
}

/* find all unique haplotypes compatible with genotypes */
void compatible_haps(unsigned int *hapvec, unsigned int genotype1, unsigned int genotype2)
{
  long int g1g2 = genotype1 ^ genotype2;
  int hapNo=1;
  unsigned int mask=1;
  for(int j=0; j<MAXLOCI; j++)
    {
      if((mask << j) & g1g2)
	{
	  hapNo = hapNo*2;
	  for(int k=(hapNo/2);k<hapNo;k++)
	    hapvec[k]=hapvec[k-(hapNo/2)];
	  for(int i=0;i<(hapNo/2);i++)
	    {
	      int allele=0; int pos=0;
	      while(pos<hapNo && allele < 2)
		{
		  if(hapvec[i] == hapvec[pos])
		    {
		      if(allele == 1)
			  hapvec[pos] = (mask << j) ^ hapvec[pos];
		      allele++;
		    }
		  pos++;
		}
	    }
	}
      else
	{
	  if((mask << j) & genotype1)
	    for(int i=0;i<hapNo;i++)
	      hapvec[i] = (mask << j) ^ hapvec[i];
	}
    }
}

/*
 * Expand haplotypes when there is missing data.
 * For k missing loci, we enumerate all 4^k possible completions
 * (2^k for hap1 x 2^k for hap2) and collect all unique diplotypes.
 *
 * The algorithm:
 * 1. Identify missing loci from missing_mask
 * 2. For each possible completion of missing values in both haplotypes
 * 3. Generate compatible haplotypes for that completion
 * 4. Add unique diplotypes to the result
 */
void expand_missing_haplotypes(struct indiv* ind, int noloci)
{
  unsigned int missing_mask = ind->missing_mask;
  unsigned int g1_base = ind->genotype1;
  unsigned int g2_base = ind->genotype2;

  /* Count missing loci */
  int num_missing = 0;
  for (int i = 0; i < noloci; i++) {
    if (missing_mask & (1U << i)) {
      num_missing++;
    }
  }

  /* Each missing locus can be 0 or 1 for each haplotype
   * So we have 4^num_missing = 2^(2*num_missing) total completions */
  int num_completions = 1 << (2 * num_missing);

  /* Temporary storage for diplotypes */
  unsigned int temp_haps[MAXHAPS * 4];  /* Plenty of space */
  int temp_count = 0;

  /* For each completion of missing data */
  for (int c = 0; c < num_completions; c++) {
    unsigned int g1 = g1_base;
    unsigned int g2 = g2_base;

    /* Apply this completion to both genotypes */
    int bit_idx = 0;
    for (int i = 0; i < noloci; i++) {
      if (missing_mask & (1U << i)) {
	/* bit_idx: 2*k and 2*k+1 control hap1 and hap2 at missing locus */
	int hap1_val = (c >> (2 * bit_idx)) & 1;
	int hap2_val = (c >> (2 * bit_idx + 1)) & 1;

	if (hap1_val) g1 |= (1U << i);
	else g1 &= ~(1U << i);

	if (hap2_val) g2 |= (1U << i);
	else g2 &= ~(1U << i);

	bit_idx++;
      }
    }

    /* Get compatible haplotypes for this completion */
    unsigned int comp_haps[MAXHAPS];
    memset(comp_haps, 0, sizeof(comp_haps));
    compatible_haps(comp_haps, g1, g2);

    int num_haps = count_haplotypes(g1, g2, noloci);

    /* Create temp indiv to sort diplotypes */
    struct indiv temp_ind;
    temp_ind.genotype1 = g1;
    temp_ind.genotype2 = g2;
    temp_ind.compHaps = comp_haps;
    temp_ind.numHaps = num_haps;
    sortDiplotypes(&temp_ind);

    /* Add unique diplotypes to our collection */
    for (int h = 0; h < num_haps; h += 2) {
      unsigned int h1 = temp_ind.compHaps[h];
      unsigned int h2 = temp_ind.compHaps[h + 1];

      /* Check if this diplotype already exists */
      int found = 0;
      for (int t = 0; t < temp_count; t += 2) {
	if ((temp_haps[t] == h1 && temp_haps[t+1] == h2) ||
	    (temp_haps[t] == h2 && temp_haps[t+1] == h1)) {
	  found = 1;
	  break;
	}
      }

      if (!found && temp_count < MAXHAPS * 4 - 2) {
	temp_haps[temp_count++] = h1;
	temp_haps[temp_count++] = h2;
      }
    }
  }

  /* Copy results to indiv structure */
  for (int i = 0; i < temp_count; i++) {
    ind->compHaps[i] = temp_haps[i];
  }
  ind->numHaps = temp_count;

  /* Store representative genotypes (first completion, for reference) */
  ind->genotype1 = g1_base;
  ind->genotype2 = g2_base;
}

/* arrange haplotypes into successive pairs (hap1,hap2) of diplotypes: even = hap1, odd = hap2 */
void sortDiplotypes(struct indiv* sampleInd)
{
  int j;
  int i=0;
  if(sampleInd->numHaps > 2)
    {
      while(i < sampleInd->numHaps)
	{
	  j=1;
	  while(j < sampleInd->numHaps)
	    {
	      if((sampleInd->compHaps[i] ^ sampleInd->compHaps[j]) ==
		 (sampleInd->genotype1 ^ sampleInd->genotype2))
		{
		  unsigned int temp = sampleInd->compHaps[i+1];
		  sampleInd->compHaps[i+1] = sampleInd->compHaps[j];
		  sampleInd->compHaps[j] = temp;
		  i = i+2;
		  break;
		}
	      else
		j++;
	    }
	}
    }
  else
    {
    sampleInd->compHaps[0] = sampleInd->genotype1;
    sampleInd->compHaps[1] = sampleInd->genotype2;
    }
}


