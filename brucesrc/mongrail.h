#define MAXINDIV 1000
#define MAXLINESZ 100000
#define MAXNAMESZ 1000
#define MAXCHRNUM 1000
#define MAXLOCI 20

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
