#define MAXINDIV 1000
#define MAXLINESZ 100000
#define MAXNAMESZ 1000
#define MAXCHRNUM 1000
#define MAXLOCI 20
#define MAXFILENMSZ 1000
#define MAXBINARY 2048

struct genotype {
  int g1, g2;
  char phase;
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
void pr_hapcounts(int** hap_counts, char chr_nm[MAXCHRNUM][MAXNAMESZ], int noChr);
