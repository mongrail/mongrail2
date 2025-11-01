 #define VERSION "2.0"

#define MAXINDIV 1000
#define MAXLINESZ 100000
#define MAXNAMESZ 1000
#define MAXCHRNUM 1000
#define MAXLOCI 10
#define MAXFILENMSZ 1000
#define MAXHAPS 2048

#define THETA_A 1.0 + 2.0*noSamplesPopA
#define THETA_B 1.0 + 2.0*noSamplesPopB
#define T1A gammln(THETA_A)
#define T1B gammln(THETA_B)
#define T4A gammln(THETA_A + 2.0)
#define T4B gammln(THETA_B + 2.0)
#define TC4A gammln(THETA_A + 1.0)
#define TC4B gammln(THETA_B + 1.0)
#define L2 log(2.0)

struct genotype {
  int g1, g2;
  char phase;
};

struct indiv
{
  unsigned int genotype1;
  unsigned int genotype2;
  unsigned int* compHaps;
  int numHaps;
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
void pr_progname();
void pr_chr_SNP_stats(int, char c[MAXCHRNUM][MAXNAMESZ],int no[MAXCHRNUM], unsigned long**);
void err_line_n(char p[MAXFILENMSZ], char q[MAXFILENMSZ], char h[MAXFILENMSZ]);
void pr_summary(char p[MAXFILENMSZ], char q[MAXFILENMSZ], char h[MAXFILENMSZ], int, int n[MAXCHRNUM],
		char c[MAXCHRNUM][MAXNAMESZ], int**, int**, int**,unsigned long**);
float gammln(float);
void get_hap_counts(unsigned int**, int**, int, int);
unsigned int count_haplotypes(unsigned int, unsigned int, int);
void compatible_haps(unsigned int*, unsigned int, unsigned int);
void sortDiplotypes(struct indiv);
void pr_hapcounts(int**, char c[MAXCHRNUM][MAXNAMESZ], int);
void pr_compatible_haps_hybrids(int, char c[MAXCHRNUM][MAXNAMESZ], struct indiv**, int, int n[MAXCHRNUM]);
double lik_a_d(int, struct indiv**, int**, unsigned int**, int*, int, int, int, enum modelType);
double lik_c(int, struct indiv**, int**, int**, unsigned int**, int*, int,  int, int);
double lik_b_e(int, struct indiv**, int**, int**, unsigned int**, int*, int, int, int, enum modelType);
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
double kahanSum(double*, int);
double max_element(double*, int);
void init_likelihood_globals(int, int);
void readDataFiles(char popAfileNm[], char popBfileNm[], char hybridfileNm[],
		   int* noSamplesPopA, int* noSamplesPopB, int* noSamplesPophybrid,
		   int* noChrom, char chr_names[MAXCHRNUM][MAXNAMESZ], int no_loci[MAXCHRNUM],
		   unsigned long** marker_positions, unsigned int*** popA_haplotypes,
		   unsigned int*** popB_haplotypes, unsigned int*** hybrid_haplotypes);

