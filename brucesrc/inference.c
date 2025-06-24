#include "mongrail.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double lik_a_d(int indivIndex, struct indiv** hybrid_indiv, int** popY_hap_counts, unsigned int** haplist, int* no_haps, int noSamplesPopY, int noChr)
{
  unsigned int hap1, hap2;
  int nhaps=0;
  unsigned int* hlist;
  int new_hap1 = 0;
  int new_hap2 = 0;
  int phi = 0;
  double logL = 0.0;
  double probL = 0.0;
  double t1, t2, t3, t4, t5;
  double thetaB = 1.0 + 2.0*noSamplesPopY;
  double l2 = log(2.0);
  t1 = gammln(thetaB);
  t4 = gammln(thetaB + 2.0);
  hlist = (unsigned int *)malloc(MAXHAPS*sizeof(unsigned int)); 
  for(int i=0; i<noChr; i++)
    {
      probL = 0.0;
      for(int k=0; k<hybrid_indiv[i][indivIndex].numHaps; k=k+2)
	{
	  hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
	  hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];
	  /* create local copy of haplist */
	  for(int h=0; h<no_haps[i]; h++)
	    hlist[h] = haplist[i][h]; 
	  nhaps = no_haps[i];
	  /* add hybrid haplotypes */
	  add_hap_lcopy(hap1,hlist,&nhaps);
	  add_hap_lcopy(hap2,hlist,&nhaps);
	  /* likelihood calculations */
	  t3 = identity2_hap(hap1, hap2)*l2; 
	  t2=0.0;
	  t5=0.0;
	  for(int j=0; j<nhaps; j++)
		t2 += gammln(popY_hap_counts[i][hlist[j]]+1.0/nhaps);
	      for(int j=0; j<nhaps; j++)
		{
		  phi = identity1_hap(hlist[j], hap1) + identity1_hap(hlist[j], hap2);
		  t5 += gammln(phi + popY_hap_counts[i][hlist[j]] + 1.0/nhaps);
		}
	      probL += exp(t1 - t2 + t3 - t4 + t5);
	}
      logL += log(probL);
    }
  free(hlist);
  return logL;
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

void add_hap_lcopy(unsigned int hap, unsigned int* hlist, int* nhaps)
{
  int hap_found = 0;
  if(*nhaps == 0)
    {
      hlist[0] = hap;
      *nhaps = 1;
    }
  else
    {
      for(int i=0; i < *nhaps; i++) 
	{
	  if(hlist[i] == hap)
	    {
	      hap_found = 1;
	      break;
	    }
	}
      if(hap_found == 0)
	{
	  hlist[*nhaps] = hap;
	  *nhaps = *nhaps + 1;
	}
    }
}

int identity2_hap(unsigned int hap1, unsigned int hap2)
{
  if(hap1==hap2)
    return 0;
  else
    return 1;
}

int identity1_hap(unsigned int hap1, unsigned int hap2)
{
  if(hap1==hap2)
    return 1;
  else
    return 0;
}
