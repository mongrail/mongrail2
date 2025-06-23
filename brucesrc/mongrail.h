#define MAXINDIV 1000
#define MAXLINESZ 100000
#define MAXNAMESZ 1000
#define MAXCHRNUM 1000
#define MAXLOCI 12
#define MAXFILENMSZ 1000
#define MAXHAPS 2048

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

int findstring (char strarray[MAXCHRNUM][MAXNAMESZ], char *strtarget, int len);
unsigned int binaryToDecimal (const char *binaryStr);
int countfilelines (char filename[]);
void file_to_array (char *target_array[MAXLINESZ], char filename[]);
int mk_chrom_list (char* file_array[MAXLINESZ], char names[MAXCHRNUM][MAXNAMESZ], int* no_loci, int nolines);
int get_genotypes (char* gstring, struct genotype* chr_locus);
void get_haplotypes (unsigned int** haplotypes, struct genotype*** gdata, int nChr, int noInd, int* nloci);
void get_positions(int noChr,int* no_loci, char** raw_data, unsigned long** positions);
void pr_haplotypes(int noChr, char chr_nm[MAXCHRNUM][MAXNAMESZ], unsigned int** popA_haps, unsigned int** popB_haps, int nsPopA, int nsPopB);
void pr_datastats(int noChr, char chr_nm[MAXCHRNUM][MAXNAMESZ], unsigned long** positions, int** popA_noIndivs, int** popB_noIndivs, int** pophybrid_noIndivs, int no_loci[MAXCHRNUM]);
void pr_progname(char* version);
void pr_chr_SNP_stats(int noChr, char chr_nm[MAXCHRNUM][MAXNAMESZ],int no_loci[MAXCHRNUM], unsigned long** positions);
void err_line_n(char pAfn[MAXFILENMSZ], char pBfn[MAXFILENMSZ], char hybfn[MAXFILENMSZ]);
void pr_summary(char popAfileNm[MAXFILENMSZ], char popBfileNm[MAXFILENMSZ], char hybridfileNm[MAXFILENMSZ], int noChr, int no_loci[MAXCHRNUM], char chr_nm[MAXCHRNUM][MAXNAMESZ], int** popA_noIndivs, int** popB_noIndivs, int** pophybrid_noIndivs,unsigned long** positions);
float gammln(float xx);
void get_hap_counts(unsigned int** haplotypes, int** hap_counts, int noChr, int no_indiv);
unsigned int count_haplotypes(unsigned int hap1, unsigned int hap2, int noloci);
void compatible_haps(unsigned int *hapvec, unsigned int genotype1, unsigned int genotype2);
void sortDiplotypes(struct indiv sampleInd);
void pr_hapcounts(int** hap_counts, char chr_nm[MAXCHRNUM][MAXNAMESZ], int noChr);
double lik_a_d(int indivIndex, struct indiv** hybrid_indiv, int** popY_hap_counts, unsigned int** haplist, int* no_haps, int noSamplesPopY, int noChr);
void add_hap(unsigned int hap, unsigned int** haplist, int* no_haps, int chrom);
int find_hap(unsigned int hap, unsigned int** haplist, int* no_haps, int chrom);
int identity2_hap(unsigned int hap1, unsigned int hap2);
int identity1_hap(unsigned int hap1, unsigned int hap2);
