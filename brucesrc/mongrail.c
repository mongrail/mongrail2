#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mongrail.h"

 

int main(int argc, char *argv[])
{
  char version[]="2.0";
  int verbose=2;
  int linepos=0;
  int no_file_lines = 0;
  int noChrom = 0;
  int noSamplesPopA = 0;
  int noSamplesPopB = 0;
  int noSamplesPophybrid = 0;
  char** raw_data_popA = NULL; 
  char** raw_data_popB = NULL;
  char** raw_data_hybrids = NULL;
  struct genotype*** popA_genotypes = NULL;
  struct genotype*** popB_genotypes = NULL;
  struct genotype*** pophybrid_genotypes = NULL;
  unsigned int** popA_haplotypes = NULL;
  unsigned int** popB_haplotypes = NULL;
  int** popA_noIndivs = NULL;
  int** popB_noIndivs = NULL;
  int** pophybrid_noIndivs = NULL;
  unsigned long** marker_positions = NULL;
  char chr_names[MAXCHRNUM][MAXNAMESZ];
  int no_loci[MAXCHRNUM] = {0};
 
  if (argc != 4)
    {
      printf("Usage: %s popA popB hybrids\n", argv[0]);
      exit(1);
    }

  pr_progname(version);
  
  /* read raw GT data into arrays of lines as strings */

  no_file_lines = countfilelines(argv[1]);
  raw_data_popA = (char**)malloc(no_file_lines*sizeof(char *));
  raw_data_popB = (char**)malloc(no_file_lines*sizeof(char *));
  raw_data_hybrids = (char**)malloc(no_file_lines*sizeof(char *));
  file_to_array(raw_data_popA,argv[1]);
  file_to_array(raw_data_popB,argv[2]);
  file_to_array(raw_data_hybrids,argv[3]);

  /* get chromosomes names, number of chromsomes (noChrom) and number of markers for each chromosome (no_loci) */  

  noChrom = mk_chrom_list (raw_data_popA,chr_names,no_loci,no_file_lines);

  /* get marker positions for each chromosome */

  marker_positions = (unsigned long**)malloc(noChrom*sizeof(unsigned long *));
  for(int i=0; i<noChrom; i++)
    marker_positions[i] = (unsigned long*)malloc(no_loci[i]*sizeof(unsigned long));

  get_positions(noChrom, no_loci, raw_data_popA, marker_positions);

  /* get individual genotypes for each population sample and counts of individuals */

  popA_genotypes = (struct genotype***)malloc(noChrom*sizeof(struct genotype *));
  for(int i=0; i<noChrom; i++)
    {
      popA_genotypes[i] = (struct genotype**)malloc(no_loci[i]*sizeof(struct genotype *));
      for(int j=0; j<no_loci[i]; j++)
	popA_genotypes[i][j] = (struct genotype*)malloc(MAXINDIV*sizeof(struct genotype ));
	
    }

  popB_genotypes = (struct genotype***)malloc(noChrom*sizeof(struct genotype *));
  for(int i=0; i<noChrom; i++)
    {
      popB_genotypes[i] = (struct genotype**)malloc(no_loci[i]*sizeof(struct genotype *));
      for(int j=0; j<no_loci[i]; j++)
	popB_genotypes[i][j] = (struct genotype*)malloc(MAXINDIV*sizeof(struct genotype ));
	
    }

    pophybrid_genotypes = (struct genotype***)malloc(noChrom*sizeof(struct genotype *));
  for(int i=0; i<noChrom; i++)
    {
      pophybrid_genotypes[i] = (struct genotype**)malloc(no_loci[i]*sizeof(struct genotype *));
      for(int j=0; j<no_loci[i]; j++)
	pophybrid_genotypes[i][j] = (struct genotype*)malloc(MAXINDIV*sizeof(struct genotype ));
	
    }

  popA_noIndivs = (int**)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    popA_noIndivs[i] = (int*)malloc(no_loci[i]*sizeof(int));

  popB_noIndivs = (int**)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    popB_noIndivs[i] = (int*)malloc(no_loci[i]*sizeof(int));

   pophybrid_noIndivs = (int**)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    pophybrid_noIndivs[i] = (int*)malloc(no_loci[i]*sizeof(int));

  linepos=0;
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<no_loci[i]; j++)
    {
      if(i==0)
	linepos=j;
      else
	{
	  linepos = 0;
	  for(int k=0; k<i; k++)
	    linepos += no_loci[k];
	  linepos+=j;
	}
      popA_noIndivs[i][j] = get_genotypes(raw_data_popA[linepos],popA_genotypes[i][j]);
    }
  noSamplesPopA = popA_noIndivs[0][0];

 linepos=0;
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<no_loci[i]; j++)
    {
      if(i==0)
	linepos=j;
      else
	{
	  linepos = 0;
	  for(int k=0; k<i; k++)
	    linepos += no_loci[k];
	  linepos+=j;
	}
      popB_noIndivs[i][j] = get_genotypes(raw_data_popB[linepos],popB_genotypes[i][j]);
    }
  noSamplesPopB = popB_noIndivs[0][0];

linepos=0;
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<no_loci[i]; j++)
    {
      if(i==0)
	linepos=j;
      else
	{
	  linepos = 0;
	  for(int k=0; k<i; k++)
	    linepos += no_loci[k];
	  linepos+=j;
	}
      pophybrid_noIndivs[i][j] = get_genotypes(raw_data_hybrids[linepos],pophybrid_genotypes[i][j]);
    }
  noSamplesPophybrid = pophybrid_noIndivs[0][0];
 
  
  /* get population haplotypes as unsigned ints */
  
  popA_haplotypes = (unsigned int **)malloc(noChrom*sizeof(unsigned int *));
  for(int i=0; i<noChrom; i++)
    popA_haplotypes[i] = (unsigned int *)malloc(2*noSamplesPopA*sizeof(unsigned int));
  get_haplotypes(popA_haplotypes,popA_genotypes,noChrom,noSamplesPopA,no_loci);

  popB_haplotypes = (unsigned int **)malloc(noChrom*sizeof(unsigned int *));
  for(int i=0; i<noChrom; i++)
    popB_haplotypes[i] = (unsigned int *)malloc(2*noSamplesPopB*sizeof(unsigned int));
  get_haplotypes(popB_haplotypes,popB_genotypes,noChrom,noSamplesPopB,no_loci);

  /* printing information to screen */
  
  if(verbose==1)
    {
      pr_haplotypes(noChrom, chr_names, popA_haplotypes, popB_haplotypes, noSamplesPopA, noSamplesPopB);
    }
  if(verbose==2)
    {
      pr_datastats(noChrom, chr_names, marker_positions, popA_noIndivs, popB_noIndivs, pophybrid_noIndivs,no_loci);
    }
  if(verbose==3)
    {
      pr_chr_SNP_stats(noChrom, chr_names,no_loci, marker_positions);
      
    } 
      return 0;
}

void pr_progname(char* version)
{
  printf("MONGRAIL v%s\nhttps://github.com/mongrail/mongrail2\n\n",version);
}
