#include "mongrail.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double lik_a_d(int indivIndex, struct indiv** hybrid_indiv, int** popY_hap_counts, unsigned int** haplist, int* no_haps, int noSamplesPopY, int noChr)
{
  unsigned int hap1, hap2;
  int nhaps=0;
  unsigned int* hlist;
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
	  /* debugging begins */
 	  /* printf("hybrid indiv: %d chrom: %d Haplotype list: ",indivIndex,i); */
	  /* for(int ha=0; ha<nhaps; ha++) */
	  /*   printf(" %d ", hlist[ha]); */
	  /* printf("\n"); */
	  /* debugging ends */
	  /* likelihood calculations */
	  t3 = identity2_hap(hap1, hap2)*l2; 
	  t2=0.0;
	  t5=0.0;
	  for(int j=0; j<nhaps; j++)
	    {
	      phi = identity1_hap(hlist[j], hap1) + identity1_hap(hlist[j], hap2);
	      t5 += gammln(phi + popY_hap_counts[i][hlist[j]] + 1.0/nhaps) - gammln(popY_hap_counts[i][hlist[j]]+1.0/nhaps);
	    }
	  probL += exp(t1+ t3 - t4 + t5);
	}
      logL += log(probL);
    }
  free(hlist);
  return logL;
}

double lik_c(int indivIndex, struct indiv** hybrid_indiv, int** popB_hap_counts, int** popA_hap_counts,
	     unsigned int** haplist, int* no_haps, int noSamplesPopB,  int noSamplesPopA, int noChr)
{
  unsigned int hap1, hap2;
  int nhaps=0;
  unsigned int* hlist;
  int phi1 = 0;
  int phi2 = 0;
  double logL = 0.0;
  double probL = 0.0;
  double thetaA = 1.0 + 2.0*noSamplesPopA;
  double thetaB = 1.0 + 2.0*noSamplesPopB;
  double lgt1A = gammln(thetaA);
  double lgt2A = gammln(1.0 + thetaA);
  double lgt1B = gammln(thetaB);
  double lgt2B = gammln(1.0 + thetaB);
  double t2A, t2B;
  double t4h1A, t4h2A, t4h1B, t4h2B;
  double p1A, p1B, p2A, p2B;
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
	  t2A = 0.0;
	  t2B = 0.0;
	  t4h1A = 0.0;
	  t4h2A = 0.0;
	  t4h1B = 0.0;
	  t4h2B = 0.0;
	  p1A = 0.0;
	  p1B = 0.0;
	  p2A = 0.0;
	  p2B = 0.0;
	  
	  for(int j=0; j<nhaps; j++)
	    {
	      t2A += gammln(popA_hap_counts[i][hlist[j]]+1.0/nhaps);
	      t2B += gammln(popB_hap_counts[i][hlist[j]]+1.0/nhaps);
	    }
	  for(int j=0; j<nhaps; j++)
	    {
	      phi1 = identity1_hap(hlist[j], hap1);
	      phi2 = identity1_hap(hlist[j], hap2);
	      t4h1A += gammln(phi1 + popA_hap_counts[i][hlist[j]] + 1.0/nhaps);
	      t4h2A += gammln(phi2 + popA_hap_counts[i][hlist[j]] + 1.0/nhaps);
	      t4h1B += gammln(phi1 + popB_hap_counts[i][hlist[j]] + 1.0/nhaps);
	      t4h2B += gammln(phi2 + popB_hap_counts[i][hlist[j]] + 1.0/nhaps);
	    }
	  p1A = exp(lgt1A - t2A - lgt2A + t4h1A);
	  p1B = exp(lgt1B - t2B - lgt2B + t4h1B);
	  p2A = exp(lgt1A - t2A - lgt2A + t4h2A);
	  p2B = exp(lgt1B - t2B - lgt2B + t4h2B);
	  probL += identity1_hap(hap1, hap2)*(p1A*p2B+p2A*p1B)+(1 - identity1_hap(hap1, hap2))*(p1A*p2B); 
	}
      logL += log(probL);
    }
  free(hlist);
  return logL;
}

double lik_b_e(int indivIndex, struct indiv** hybrid_indiv, int** popB_hap_counts, int** popA_hap_counts,
	       unsigned int** haplist, int* no_haps, int noSamplesPopB,  int noSamplesPopA, int noChr, char model)
{
  unsigned int hap1, hap2;
  int nhaps=0;
  unsigned int* hlist;
  int phi1 = 0;
  int phi2 = 0;
  double logL = 0.0;
  double probL = 0.0;
  double thetaA = 1.0 + 2.0*noSamplesPopA;
  double thetaB = 1.0 + 2.0*noSamplesPopB;
  double lgt1A = gammln(thetaA);
  double lgt2A = gammln(1.0 + thetaA);
  double lgt1B = gammln(thetaB);
  double lgt2B = gammln(1.0 + thetaB);
  double t2A, t2B;
  double t4h1A, t4h2A, t4h1B, t4h2B;
  double p1A, p1B, p2A, p2B;
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
	  t2A = 0.0;
	  t2B = 0.0;
	  t4h1A = 0.0;
	  t4h2A = 0.0;
	  t4h1B = 0.0;
	  t4h2B = 0.0;
	  p1A = 0.0;
	  p1B = 0.0;
	  p2A = 0.0;
	  p2B = 0.0;
	  
	  for(int j=0; j<nhaps; j++)
	    {
	      t2A += gammln(popA_hap_counts[i][hlist[j]]+1.0/nhaps);
	      t2B += gammln(popB_hap_counts[i][hlist[j]]+1.0/nhaps);
	    }
	  for(int j=0; j<nhaps; j++)
	    {
	      phi1 = identity1_hap(hlist[j], hap1);
	      phi2 = identity1_hap(hlist[j], hap2);
	      t4h1A += gammln(phi1 + popA_hap_counts[i][hlist[j]] + 1.0/nhaps);
	      t4h2A += gammln(phi2 + popA_hap_counts[i][hlist[j]] + 1.0/nhaps);
	      t4h1B += gammln(phi1 + popB_hap_counts[i][hlist[j]] + 1.0/nhaps);
	      t4h2B += gammln(phi2 + popB_hap_counts[i][hlist[j]] + 1.0/nhaps);
	    }
	  if(model=='a')
	    {
	      p1A = exp(lgt1A - t2A - lgt2A + t4h1A);
	      p2A = exp(lgt1A - t2A - lgt2A + t4h2A);
	    }
	  else
	    {
	      p1B = exp(lgt1B - t2B - lgt2B + t4h1B);
	      p2B = exp(lgt1B - t2B - lgt2B + t4h2B);
	    }
	  probL += identity1_hap(hap1, hap2)*(p1A*p2B+p2A*p1B)+(1 - identity1_hap(hap1, hap2))*(p1A*p2B); 
	}
      logL += log(probL);
    }
  return 0.0;
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

int match_sub_hap(unsigned int hyb_hap, unsigned int pop_hap, unsigned int ancvec)
{
  if((hyb_hap & ancvec)==(pop_hap & ancvec))
    return 1;
  else
    return 0;
}

int get_marginal_hap_counts(unsigned int hyb_hap, int population, unsigned int ancvec, int** popY_hap_counts, unsigned int** haplist, int* no_haps, int chrom)
{
  int mhcount=0;
  for(int i=0; i<no_haps[chrom]; i++)
    if(match_sub_hap(hyb_hap,haplist[chrom][i],ancvec))
      mhcount += popY_hap_counts[chrom][haplist[chrom][i]];
  return mhcount;
}

unsigned int get_anc_complement(unsigned int ancvec,int noloci)
{
  unsigned int cav = ~ancvec; 
  return cav - intpow(2,32) + intpow(2,noloci);;
}


double prPopSwitch(double r, unsigned long d)
{
  return 1.0 - cosh(r*d)*exp(-r*d); 
}

double prQ(unsigned int anc, int noloci, double rec_rate_bp, unsigned long* position)
{
  double t1 = 0.0;
  double t2 = 0.0;
  double Q = 0.0;
  double p = 0.0;
  double pprod = 1.0;
  unsigned int mask1 = intpow(2,noloci) - 2;
  p = prPopSwitch(rec_rate_bp,position[0]);
  if((anc & mask1) == 0)
    {
      t1 = 0.5*(1.0 - p);
      t2 = 0.5*p;
    }
  else
    {
      t1 = 0.5*p;
      t2 = 0.5*(1.0 - p);
    }
  int curr_anc = 0;
  int prev_anc = ((anc >> 0)&1);
  for(int i=1; i<noloci; i++)
    {
      curr_anc = ((anc >> i)&1);
      p = prPopSwitch(rec_rate_bp,position[i]-position[i-1]); 
      if(curr_anc != prev_anc)
	pprod *= p;
      else
	pprod *= 1.0 - p;
      prev_anc = ((anc >> i)&1);
    }
  Q = pprod*(t1 + t2);
  return Q;
}

