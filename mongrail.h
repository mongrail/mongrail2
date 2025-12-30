 #define VERSION "2.1"

#define MAXINDIV 1000
#define MAXLINESZ 100000
#define MAXNAMESZ 1000
#define MAXCHRNUM 1000
#define MAXLOCI 32        /* max SNPs per chromosome (limited by 32-bit haplotype representation) */
#define MAXFILENMSZ 1000
#define MAXHAPS 65536     /* max haplotypes to track per chromosome */
#define DEFAULT_MAX_RECOMB 4  /* default max recombinations for sparse enumeration */

/* Cache structure for precomputed values */
struct likelihood_cache {
  int initialized;
  int noChr;
  int* no_loci;
  int max_recomb;           /* max recombinations for sparse enumeration */
  double recomb_rate;       /* recombination rate per bp */
  unsigned long** marker_positions;  /* marker positions for Q computation */
  /* U(hap) cache per chromosome: U_cache[chrom][hap] (-1 = not computed) */
  double** U_cache;
  /* log P(hap|popA) cache: logP_A[chrom][hap] */
  double** logP_A;
  /* log P(hap|popB) cache: logP_B[chrom][hap] */
  double** logP_B;
  /* precomputed constants */
  double log2;
};

#define THETA_A 1.0 + 2.0*noSamplesPopA
#define THETA_B 1.0 + 2.0*noSamplesPopB
#define T1A gammln(THETA_A)
#define T1B gammln(THETA_B)
#define T4A gammln(THETA_A + 2.0)
#define T4B gammln(THETA_B + 2.0)
#define TC4A gammln(THETA_A + 1.0)
#define TC4B gammln(THETA_B + 1.0)
#define L2 log(2.0)

/* Special value for missing alleles */
#define MISSING_ALLELE -1

struct genotype {
  int g1, g2;  /* use -1 (MISSING_ALLELE) for missing data */
  char phase;
};

struct indiv
{
  unsigned int genotype1;
  unsigned int genotype2;
  unsigned int* compHaps;
  int numHaps;
  unsigned int missing_mask;  /* bitmask: bit i = 1 means locus i is missing */
};

enum modelType { MODEL_A, MODEL_B, MODEL_C, MODEL_D, MODEL_E, MODEL_F };

int findstring (char s[MAXCHRNUM][MAXNAMESZ], char*, int);
unsigned int binaryToDecimal (const char*);
int countfilelines (char*);
void file_to_array (char* t[MAXLINESZ], char*);
int mk_chrom_list (char* f[MAXLINESZ], char n[MAXCHRNUM][MAXNAMESZ], int*, int);
int get_genotypes (char*, struct genotype*);
void get_haplotypes (unsigned int**, struct genotype***, int, int, int*);
void get_positions(int,int*, char**, unsigned long**);
void pr_haplotypes(int, char c[MAXCHRNUM][MAXNAMESZ], unsigned int**, unsigned int**, int, int);
void pr_datastats(int, char c[MAXCHRNUM][MAXNAMESZ], unsigned long**, int**, int**, int**, int n[MAXCHRNUM]);
void pr_chr_SNP_stats(int, char c[MAXCHRNUM][MAXNAMESZ],int no[MAXCHRNUM], unsigned long**);
void err_line_n(char p[MAXFILENMSZ], char q[MAXFILENMSZ], char h[MAXFILENMSZ]);
void pr_summary(char p[MAXFILENMSZ], char q[MAXFILENMSZ], char h[MAXFILENMSZ], int, int n[MAXCHRNUM],
		char c[MAXCHRNUM][MAXNAMESZ], int**, int**, int**,unsigned long**);
float gammln(float);
void get_hap_counts(unsigned int**, int**, int, int);
unsigned int count_haplotypes(unsigned int, unsigned int, int);
void compatible_haps(unsigned int*, unsigned int, unsigned int);
void sortDiplotypes(struct indiv*);
void pr_hapcounts(int**, char c[MAXCHRNUM][MAXNAMESZ], int);
void pr_compatible_haps_hybrids(int, char c[MAXCHRNUM][MAXNAMESZ], struct indiv**, int, int n[MAXCHRNUM]);
double lik_a_d(int, struct indiv**, int**, unsigned int**, int*, int, int, int, enum modelType);
double lik_c(int, struct indiv**, int**, int**, unsigned int**, int*, int,  int, int);
double lik_b_e(int, struct indiv**, int**, int**, unsigned int**, int*, int, int, int,
	       int*, unsigned long**, double, enum modelType);
double lik_f(int, struct indiv**, int**, int**, unsigned int**, int*, int, int, int,
	     int*, unsigned long**, double);
/* Helper functions for Models B and E */
int get_distinct_subhaps(unsigned int, int**, unsigned int**, int*, int, unsigned int*, int*);
double log_prob_single_hap(unsigned int, int**, unsigned int**, int*, int, int);
double log_prob_subhap(unsigned int, unsigned int, int**, unsigned int**, int*, int, int);
double log_U_given_z(unsigned int, unsigned int, int, int**, int**, unsigned int**, int*, int, int, int);
double compute_U(unsigned int, int, double, unsigned long*, int**, int**, unsigned int**, int*, int, int, int);
void add_hap(unsigned int, unsigned int**, int*, int);
int identity2_hap(unsigned int, unsigned int);
int identity1_hap(unsigned int, unsigned int);
int match_sub_hap(unsigned int, unsigned int, unsigned int);
int get_marginal_hap_counts(unsigned int, int, unsigned int, int**, unsigned int**, int*, int);
double prPopSwitch(double, unsigned long);
double prQ(unsigned int, int, double, unsigned long*);
float gammln(float);
unsigned int binaryToDecimal(const char*);
void prn_binary(unsigned int, int);
unsigned int intpow(unsigned int, unsigned int);
unsigned int get_anc_complement(unsigned int,int);
/* Missing data handling */
void get_haplotypes_with_missing(unsigned int**, unsigned int**, struct genotype***, int, int, int*);
void expand_missing_haplotypes(struct indiv*, int);
double kahanSum(double*, int);
double max_element(double*, int);
void init_likelihood_globals(int, int);
/* Cache functions */
struct likelihood_cache* cache_init(int noChr, int* no_loci, double recomb_rate,
				    unsigned long** marker_positions,
				    int** popA_hap_counts, int** popB_hap_counts,
				    unsigned int** haplist, int* no_haps,
				    int noSamplesPopA, int noSamplesPopB,
				    int max_recomb);
/* Sparse enumeration helpers */
int count_switches(unsigned int z, int noloci);
void cache_free(struct likelihood_cache* cache);
double compute_U_cached(struct likelihood_cache* cache, unsigned int hap, int chrom,
			int** popA_hap_counts, int** popB_hap_counts,
			unsigned int** haplist, int* no_haps,
			int noSamplesPopA, int noSamplesPopB);
double lik_b_e_cached(struct likelihood_cache* cache, int indivIndex,
		      struct indiv** hybrid_indiv,
		      int** popB_hap_counts, int** popA_hap_counts,
		      unsigned int** haplist, int* no_haps,
		      int noSamplesPopB, int noSamplesPopA, int noChr,
		      enum modelType model);
double lik_f_cached(struct likelihood_cache* cache, int indivIndex,
		    struct indiv** hybrid_indiv,
		    int** popB_hap_counts, int** popA_hap_counts,
		    unsigned int** haplist, int* no_haps,
		    int noSamplesPopB, int noSamplesPopA, int noChr);
void readDataFiles(char popAfileNm[], char popBfileNm[], char hybridfileNm[],
		   int* noSamplesPopA, int* noSamplesPopB, int* noSamplesPophybrid,
		   int* noChrom, char chr_names[MAXCHRNUM][MAXNAMESZ], int no_loci[MAXCHRNUM],
		   unsigned long*** marker_positions, unsigned int*** popA_haplotypes,
		   unsigned int*** popB_haplotypes, unsigned int*** hybrid_haplotypes,
		   unsigned int*** hybrid_missing_masks);

/* VCF file reading (requires htslib) */
int readVCFFiles(char popAfileNm[], char popBfileNm[], char hybridfileNm[],
		 int* noSamplesPopA, int* noSamplesPopB, int* noSamplesPophybrid,
		 int* noChrom, char chr_names[MAXCHRNUM][MAXNAMESZ], int no_loci[MAXCHRNUM],
		 unsigned long*** marker_positions, unsigned int*** popA_haplotypes,
		 unsigned int*** popB_haplotypes, unsigned int*** hybrid_haplotypes,
		 unsigned int*** hybrid_missing_masks);

