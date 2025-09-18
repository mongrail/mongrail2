#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "bit.h"


/* get number of polymorphic sites and number of haplotypes for hybrid indivs */
unsigned int count_haplotypes(unsigned int hap1, unsigned int hap2, int noloci)
{
  unsigned int hap1hap2;
  unsigned int no1bits=0;
  hap1hap2 = hap1^hap2; /* 1 if 1^0 or 0^1, 0 otherwise */
  for (int i = 0; i < 32 && i < noloci; ++i)
    if (hap1hap2 >> i & 0x1) no1bits++; /* if ith bit of hap1hap2 = 1 is polymorphic */
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




int main(int argc, char *argv[])
{
   struct indiv sampleInd;
   sampleInd.compHaps = calloc(MAXHAPS,sizeof(unsigned int));
   sampleInd.genotype1 = 3;
   sampleInd.genotype2 = 0;
   compatible_haps(sampleInd.compHaps, sampleInd.genotype1, sampleInd.genotype2);
   sampleInd.numHaps = count_haplotypes(1, 0, 2);
   printf("hcount: %d\n",sampleInd.numHaps);
   for(int i=0; i<sampleInd.numHaps; i++)
     printf("hap: %d is %u\n",i,sampleInd.compHaps[i]);
}
