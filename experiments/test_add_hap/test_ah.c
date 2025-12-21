#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
// #include "mongrail.h"

#define MAXHAPS 16

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

int main(int argc, char *argv[])
{
  unsigned int hap;
  int chrom = 0;
  int noChrom = 2;
  unsigned int** haplist;
  int* nohaps;
  haplist = (unsigned int **)malloc(noChrom*sizeof(unsigned int *));
  for(int i=0; i<noChrom; i++)
    haplist[i] = (unsigned int *)malloc(MAXHAPS*sizeof(unsigned int));
  nohaps = (int *)malloc(noChrom*sizeof(int));
  nohaps[0]=2;
  nohaps[1]=3;
  haplist[0][0] = 1;
  haplist[0][1] = 3;
  haplist[1][0] = 1;
  haplist[1][1] = 2;
  haplist[1][2] = 3;
  add_hap(3, haplist, nohaps, 0);
  add_hap(4, haplist, nohaps, 1);
  for(int i=0; i < noChrom; i++)
    for(int j=0; j<nohaps[i]; j++)
      printf("Chrom: %d NoHaps: %d Hap: %u\n",i,nohaps[i],haplist[i][j]);
		  
  printf("testing add_hap\n");

}

