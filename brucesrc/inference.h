#include "mongrail.h"

double lik_a(struct indiv** hybrid_indiv, int** popB_hap_counts, unsigned int** haplist, int* no_haps);
void add_hap(unsigned int hap, unsigned int** haplist, int* no_haps, int chrom);
int find_hap(unsigned int hap, unsigned int** haplist, int* no_haps, int chrom);
