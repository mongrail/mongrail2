
struct indiv
{
  unsigned int genotype1;
  unsigned int genotype2;
  unsigned int* compHaps;
  int numHaps;
};

double lik_a(unsigned int hap1, unsigned int hap2, int** popB_hap_counts, unsigned int** haplist, int* no_haps);
void add_hap(unsigned int hap, unsigned int** haplist, int* no_haps, int chrom);
