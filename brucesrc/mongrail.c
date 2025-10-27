#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mongrail.h"

int verbose=1;

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

  
  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s popA popB hybrids\n", argv[0]);
      exit(1);
    }

  strncpy(popAfileNm,argv[1],MAXFILENMSZ);
  strncpy(popBfileNm,argv[2],MAXFILENMSZ);
  strncpy(hybridfileNm,argv[3],MAXFILENMSZ);

  
  pr_progname();
  readDataFiles(popAfileNm, popBfileNm, hybridfileNm, &noSamplesPopA, &noSamplesPopB, &noSamplesPophybrid,
		&noChrom, chr_names, no_loci, marker_positions, popA_haplotypes, popB_haplotypes,hybrid_haplotypes);

  /* summarize pop samples as haplotype counts in array of length MAXHAPS with base_10 haplotype as array index */ 
  
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
  get_hap_counts(popA_haplotypes, popA_hap_counts, noChrom, noSamplesPopA);
  get_hap_counts(popB_haplotypes, popB_hap_counts, noChrom, noSamplesPopB);

  /* get hybrid individuals */
  hybrid_indiv = (struct indiv** )malloc(noChrom*sizeof(struct indiv*));
  for(int i=0; i<noChrom; i++)
    hybrid_indiv[i] = (struct indiv*)malloc(noSamplesPophybrid*sizeof(struct indiv));
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<noSamplesPophybrid; j++)
      hybrid_indiv[i][j].compHaps = calloc(MAXHAPS,sizeof(unsigned int));
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
	/* debugging begins */
	/* printf("hybrid indiv: %d chrom: %d compatible haploypes: ",j,i); */
	/* for(int nh=0; nh < hybrid_indiv[i][j].numHaps; nh++) */
	/*   printf(" %u ",hybrid_indiv[i][j].compHaps[nh]); */
	/* printf("\n"); */
	/* debugging ends */
	
      }

  /* get list of all possible unique haplotypes for each chromosome in pop A and B samples */

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

 
  /* debugging commented out bug is here! */
  /* seg fault processing chromosome 8 */
  /* seems to be related to size parameter of calloc */

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

    /* add all possible haplotypes in hybrids */
    for(int z=0; z<noChrom; z++)
      for(int j=0; j<noSamplesPophybrid; j++)
	for(int h=0; h<hybrid_indiv[z][j].numHaps; h=h+2)
	  {
	    add_hap(hybrid_indiv[z][j].compHaps[h], haplist, nohaps, z);
	    add_hap(hybrid_indiv[z][j].compHaps[h+1], haplist, nohaps, z);
	  }    
  

  /* get intermarker distances in units of bps */
  


  /*  debugging marginal haplotype counts */
  /* printf("\nmarginal counts: %d\n",get_marginal_hap_counts(2,1,0,popB_hap_counts, haplist,nohaps,2)); */
  /* printf("counts for hap 2 in popB at chrom 2: %d\n",popB_hap_counts[2][2]); */
  
  

  /* debugging likelihoods */

  printf("Indiv     Md                Ma                Mc\n");    
   for(int i=0; i<noSamplesPophybrid; i++)
     {
       printf("%5d      %.5f",i,lik_a_d(i,hybrid_indiv,popA_hap_counts,haplist,nohaps,noSamplesPopA,noSamplesPopB,noChrom,MODEL_D));
       printf("         %.5f",lik_a_d(i,hybrid_indiv,popB_hap_counts,haplist,nohaps,noSamplesPopA,noSamplesPopB,noChrom,MODEL_A));
       printf("         %.5f\n",lik_c(i,hybrid_indiv,popB_hap_counts,popA_hap_counts,haplist,nohaps,noSamplesPopB,noSamplesPopA,noChrom));
     }


   /*    printf("Indiv: %d Ma logL: %f Md logL: %f Mc logL: %f\n",i,lik_a_d(i,hybrid_indiv,popB_hap_counts,haplist,nohaps,noSamplesPopB,noChrom),lik_a_d(i,hybrid_indiv,popA_hap_counts,haplist,nohaps,noSamplesPopA,noChrom), */

  
      return 0;
}

void pr_progname()
{
  printf("\nMONGRAIL v%s\nhttps://github.com/mongrail/mongrail2\n\n",VERSION);
}
