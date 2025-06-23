#include "mongrail.h"

double lik_a(unsigned int hap1, unsigned int hap2, int** popB_hap_counts, unsigned int** haplist, int* no_haps)
{
  return 1.0;
}

void add_hap(unsigned int hap, unsigned int** haplist, int* no_haps, int chrom)
{
  int hap_found = 0;
  if(no_haps[chrom] == 0)
    {
      haplist[chrom][0] = hap;
      no_haps[chrom] = 1;
    }
  else
    {
      for(int i=0; i < no_haps[chrom]; i++) 
	{
	  if(haplist[chrom][i] == hap)
	    {
	      hap_found = 1;
	      break;
	    }
	}
      if(hap_found == 0)
	{
	  haplist[chrom][no_haps[chrom]] = hap;
	  no_haps[chrom] = no_haps[chrom] + 1;
	}
    }
}

