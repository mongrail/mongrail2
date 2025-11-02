#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mongrail.h"

int verbose=4;

int main(int argc, char *argv[])
{
  char popAfileNm[MAXFILENMSZ];
  char popBfileNm[MAXFILENMSZ];
  char hybridfileNm[MAXFILENMSZ];
  


  int noChrom = 0;
  int noSamplesPopA = 0;
  int noSamplesPopB = 0;
  int noSamplesPophybrid = 0;
  unsigned int** popA_haplotypes = NULL;
  unsigned int** popB_haplotypes = NULL;
  unsigned int** hybrid_haplotypes = NULL;
  int** popA_hap_counts;
  int** popB_hap_counts;
  unsigned int** haplist;
  int* nohaps;
  struct indiv** hybrid_indiv = NULL;
  unsigned long** marker_positions = NULL;
  char chr_names[MAXCHRNUM][MAXNAMESZ];
  int no_loci[MAXCHRNUM] = {0};
  double recomb_rate_per_bp = 0.000001;
  unsigned int pop_anc1;
  int pop_anc2;

  pr_progname();
  
  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s popA popB hybrids\n", argv[0]);
      exit(1);
    }

  /* read input file names from command line */
  strncpy(popAfileNm,argv[1],MAXFILENMSZ);
  strncpy(popBfileNm,argv[2],MAXFILENMSZ);
  strncpy(hybridfileNm,argv[3],MAXFILENMSZ);

  /* read data files and create haplotype data structures */
  readDataFiles(popAfileNm, popBfileNm, hybridfileNm, &noSamplesPopA, &noSamplesPopB, &noSamplesPophybrid,
		&noChrom, chr_names, no_loci, marker_positions, &popA_haplotypes, &popB_haplotypes, &hybrid_haplotypes);

  /* allocate memory for haplotype counts of populations A and B */
  popA_hap_counts = (int **)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    popA_hap_counts[i] = (int *)malloc(MAXHAPS*sizeof(int));
  popB_hap_counts = (int **)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    popB_hap_counts[i] = (int *)malloc(MAXHAPS*sizeof(int));
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<MAXHAPS; j++)
      {
	popA_hap_counts[i][j] = 0;
	popB_hap_counts[i][j] = 0;
      }

  /* summarize pop samples as haplotype counts in array of length MAXHAPS with base_10 haplotype as array index */ 
  get_hap_counts(popA_haplotypes, popA_hap_counts, noChrom, noSamplesPopA);
  get_hap_counts(popB_haplotypes, popB_hap_counts, noChrom, noSamplesPopB);

  /* create data structures for hybrid individuals */
  hybrid_indiv = (struct indiv** )malloc(noChrom*sizeof(struct indiv*));
  for(int i=0; i<noChrom; i++)
    hybrid_indiv[i] = (struct indiv*)malloc(noSamplesPophybrid*sizeof(struct indiv));
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<noSamplesPophybrid; j++)
      hybrid_indiv[i][j].compHaps = calloc(MAXHAPS,sizeof(unsigned int));

  /* assign genotypes to hybrid individuals and determine compatible haplotypes */
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<noSamplesPophybrid; j++)
      {
	hybrid_indiv[i][j].genotype1 = hybrid_haplotypes[i][j];
	hybrid_indiv[i][j].genotype2 = hybrid_haplotypes[i][j+noSamplesPophybrid];
      }
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<noSamplesPophybrid; j++)
      {
	hybrid_indiv[i][j].numHaps = count_haplotypes(hybrid_indiv[i][j].genotype1,hybrid_indiv[i][j].genotype2,no_loci[i]);
	compatible_haps(hybrid_indiv[i][j].compHaps, hybrid_indiv[i][j].genotype1, hybrid_indiv[i][j].genotype2);
	sortDiplotypes(hybrid_indiv[i][j]);
      }

  /* allocate memory for list (and total number) of unique haplotypes for each chromosome */
  haplist = (unsigned int **)malloc(noChrom*sizeof(unsigned int *));
  for(int i=0; i<noChrom; i++)
    {
      haplist[i] = (unsigned int *)malloc(MAXHAPS*sizeof(unsigned int));
      if (haplist[i] == NULL)
	{
            perror("Memory allocation failed");
            return EXIT_FAILURE;
	}
    }
  nohaps = (int *)calloc(4,noChrom*sizeof(int));
 
  /* get list of all possible unique haplotypes for each chromosome in pop A and B samples */
    for(int z=0; z<noChrom; z++)
    {
      for(int j=0; j < 2*noSamplesPopA; j++)
	{
	  add_hap(popA_haplotypes[z][j], haplist, nohaps, z);
	} 
      for(int j=0; j < 2*noSamplesPopB; j++)
	{
	  add_hap(popB_haplotypes[z][j], haplist, nohaps, z);
	} 
    }

    /* add all possible haplotypes in hybrids to list of unique haplotypes */
    for(int z=0; z<noChrom; z++)
      for(int j=0; j<noSamplesPophybrid; j++)
	for(int h=0; h<hybrid_indiv[z][j].numHaps; h=h+2)
	  {
	    add_hap(hybrid_indiv[z][j].compHaps[h], haplist, nohaps, z);
	    add_hap(hybrid_indiv[z][j].compHaps[h+1], haplist, nohaps, z);
	  }    
  

  /* debugging likelihoods */

  printf("Indiv     Md                Ma                Mc\n");    
   for(int i=0; i<noSamplesPophybrid; i++)
     {
       printf("%5d      %.5f",i,lik_a_d(i,hybrid_indiv,popA_hap_counts,haplist,nohaps,noSamplesPopA,noSamplesPopB,noChrom,MODEL_D));
       printf("         %.5f",lik_a_d(i,hybrid_indiv,popB_hap_counts,haplist,nohaps,noSamplesPopA,noSamplesPopB,noChrom,MODEL_A));
       printf("         %.5f\n",lik_c(i,hybrid_indiv,popB_hap_counts,popA_hap_counts,haplist,nohaps,noSamplesPopB,noSamplesPopA,noChrom));
     }

      return 0;
}

void pr_progname()
{
  printf("\nMONGRAIL v%s\nhttps://github.com/mongrail/mongrail2\n\n",VERSION);
}
