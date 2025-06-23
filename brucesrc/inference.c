#include "mongrail.h"
#include <math.h>
#include <stdio.h>

double lik_a_d(int indivIndex, struct indiv** hybrid_indiv, int** popY_hap_counts, unsigned int** haplist, int* no_haps, int noSamplesPopY, int noChr)
{
  unsigned int hap1, hap2;
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
  for(int i=0; i<noChr; i++)
    {
      probL = 0.0;
      for(int k=0; k<hybrid_indiv[i][indivIndex].numHaps; k=k+2)
	{
	  hap1 = hybrid_indiv[i][indivIndex].compHaps[k];
	  hap2 = hybrid_indiv[i][indivIndex].compHaps[k+1];
	  new_hap1 = find_hap(hap1, haplist, no_haps, i);
	  new_hap2 = find_hap(hap2, haplist, no_haps, i);
	  t3 = identity2_hap(hap1, hap2)*l2; 
	  if(new_hap1&&new_hap2)
	    {
	      t2=0.0;
	      t5=0.0;
	      for(int j=0; j<no_haps[i]; j++)
		t2 += gammln(popY_hap_counts[i][haplist[i][j]]+1.0/no_haps[i]);
	      for(int j=0; j<no_haps[i]; j++)
		{
		  phi = identity1_hap(haplist[i][j], hap1) + identity1_hap(haplist[i][j], hap2);
		  t5 += gammln(phi + popY_hap_counts[i][haplist[i][j]] + 1.0/no_haps[i]);
		}
	      probL += exp(t1 - t2 + t3 - t4 + t5);
	    }
	  else
	    if((new_hap1)||(new_hap2))
	      {
		t2=0.0;
		t5=0.0;
		for(int j=0; j<no_haps[i]; j++)
		  t2 += gammln(popY_hap_counts[i][haplist[i][j]]+1.0/(no_haps[i]+1));
		t2 += 1.0/(no_haps[i]+1);
		for(int j=0; j<no_haps[i]; j++)
		  {
		    phi = identity1_hap(haplist[i][j], hap1) + identity1_hap(haplist[i][j], hap2);
		    t5 += gammln(phi + popY_hap_counts[i][haplist[i][j]] + 1.0/(no_haps[i]+1));
		  }
		t5 += gammln(1.0 + 1.0/(no_haps[i]+1));
	      probL += exp(t1 - t2 + t3 - t4 + t5);
	      }
	  else
	    {
	      t2=0.0;
	      t5=0.0;
	      for(int j=0; j<no_haps[i]; j++)
		t2 += gammln(popY_hap_counts[i][haplist[i][j]]+1.0/(no_haps[i]+2));
	      t2 += 2.0/(no_haps[i]+2);
	      for(int j=0; j<no_haps[i]; j++)
		{
		  phi = identity1_hap(haplist[i][j], hap1) + identity1_hap(haplist[i][j], hap2);
		  t5 += gammln(phi + popY_hap_counts[i][haplist[i][j]] + 1.0/(no_haps[i]+2));
		}
		t5 += 2.0*gammln(1.0 + 1.0/(no_haps[i]+2));
	      probL += exp(t1 - t2 + t3 - t4 + t5);
	    }
	}	
      logL += log(probL);
      //      printf("logL: %f\n",logL);
    }
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

int find_hap(unsigned int hap, unsigned int** haplist, int* no_haps, int chrom)
{
  int hap_found = 0;
  for(int i=0; i < no_haps[chrom]; i++) 
    {
      if(haplist[chrom][i] == hap)
	{
	  hap_found = 1;
	  break;
	}
    }
  return hap_found;
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
