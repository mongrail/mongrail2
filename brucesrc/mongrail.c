#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mongrail.h"
#include "algorithms.h"
 

int main(int argc, char *argv[])
{
  char version[]="2.0";
  char popAfileNm[MAXFILENMSZ];
  char popBfileNm[MAXFILENMSZ];
  char hybridfileNm[MAXFILENMSZ];
  int verbose=4;
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
  unsigned int** hybrid_haplotypes = NULL;
  int** popA_hap_counts;
  int** popB_hap_counts;
  unsigned int** haplist;
  int* nohaps;
  int** popA_noIndivs = NULL;
  int** popB_noIndivs = NULL;
  int** pophybrid_noIndivs = NULL;
  struct indiv** hybrid_indiv = NULL;
  unsigned long** marker_positions = NULL;
  char chr_names[MAXCHRNUM][MAXNAMESZ];
  int no_loci[MAXCHRNUM] = {0};
 
  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s popA popB hybrids\n", argv[0]);
      exit(1);
    }

  strncpy(popAfileNm,argv[1],MAXFILENMSZ);
  strncpy(popBfileNm,argv[2],MAXFILENMSZ);
  strncpy(hybridfileNm,argv[3],MAXFILENMSZ);
  
  pr_progname(version);
  
  /* read raw GT data into arrays of lines as strings */
  /* total number of lines in popA.GT */ 
  no_file_lines = countfilelines(popAfileNm);
  /* error check: # lines equals # loci. Must be equal in popA.GT, popB.GT and hybrids.GT */
  err_line_n(popAfileNm, popBfileNm, hybridfileNm);
  /* create arrays to contain lines read from each input file */
  raw_data_popA = (char**)malloc(no_file_lines*sizeof(char *));
  raw_data_popB = (char**)malloc(no_file_lines*sizeof(char *));
  raw_data_hybrids = (char**)malloc(no_file_lines*sizeof(char *));
  /* fill arrays with lines from each input file */
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

  hybrid_haplotypes = (unsigned int **)malloc(noChrom*sizeof(unsigned int *));
  for(int i=0; i<noChrom; i++)
    hybrid_haplotypes[i] = (unsigned int *)malloc(2*noSamplesPophybrid*sizeof(unsigned int));
  get_haplotypes(hybrid_haplotypes,pophybrid_genotypes,noChrom,noSamplesPophybrid,no_loci);
  
  popA_hap_counts = (int **)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    popA_hap_counts[i] = (int *)malloc(MAXHAPS*sizeof(int));
  for(int i=0; i<noChrom; i++)
    for(int j=0; j<MAXHAPS; j++)
      popA_hap_counts[i][j] = 0;
  get_hap_counts(popA_haplotypes, popA_hap_counts, noChrom, noSamplesPopA);
  
  popB_hap_counts = (int **)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    popB_hap_counts[i] = (int *)malloc(MAXHAPS*sizeof(int));
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
      }

  /* debugging */
  for(int i=0; i<noChrom; i++)
    {
      printf("Chr: %s\n",chr_names[i]);
      for(int j=0; j<noSamplesPophybrid; j++)
	{
	  printf("Indiv: %d Numhaps: %d compatible haps:\n",j,hybrid_indiv[i][j].numHaps);
	  for(int k=0; k<hybrid_indiv[i][j].numHaps; k=k+2)
	    {
	      prn_binary(hybrid_indiv[i][j].compHaps[k],no_loci[i]);
	      printf(" ");
	      prn_binary(hybrid_indiv[i][j].compHaps[k+1],no_loci[i]);
	      printf("\n");
	    }

	  
	  /* printf("Indiv: %d Numhaps: %d compatible haps:\n",j,hybrid_indiv[i][j].numHaps); */
	  /* for(int k=0; k<hybrid_indiv[i][j].numHaps; k++) */
	  /*   printf(" %u ",hybrid_indiv[i][j].compHaps[k]); */
	  /* printf("\n"); */
	}
    }
  
  
  /* get list of all possible unique haplotypes for each chromosome in pop A and B samples */

  haplist = (unsigned int **)malloc(noChrom*sizeof(unsigned int *));
  for(int i=0; i<noChrom; i++)
    haplist[i] = (unsigned int *)malloc(MAXHAPS*sizeof(unsigned int));
  nohaps = (int *)malloc(noChrom*sizeof(int));
  for(int i=0; i<noChrom; i++)
    {
      for(int j=0; j < 2*popA_noIndivs[i][0]; j++)
	{
	  add_hap(popA_haplotypes[i][j], haplist, nohaps, i);
	}
      for(int j=0; j < 2*popB_noIndivs[i][0]; j++)
	{
	  add_hap(popB_haplotypes[i][j], haplist, nohaps, i);
	}
    }

  for(int i=0; i<noSamplesPophybrid; i++)
    printf("Indiv: %d Ma logL: %f Md logL: %f\n",i,lik_a_d(i,hybrid_indiv,popB_hap_counts,haplist,nohaps,noSamplesPopB,noChrom),lik_a_d(i,hybrid_indiv,popA_hap_counts,haplist,nohaps,noSamplesPopA,noChrom));

  
  /* printing information to screen */

  pr_summary(popAfileNm, popBfileNm, hybridfileNm, noChrom, no_loci, chr_names, popA_noIndivs, popB_noIndivs, pophybrid_noIndivs, marker_positions);

  /* for(int i=0; i<noChrom; i++) */
  /*   printf("Chr: %s nohaps[%d]:%d\n",chr_names[i],i,nohaps[i]); */
  /* printf("\npopulation A:\n"); */
  /* pr_hapcounts(popA_hap_counts, chr_names, noChrom); */
  /* printf("\npopulation B:\n"); */
  /* pr_hapcounts(popB_hap_counts, chr_names, noChrom); */
  
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
  printf("\nMONGRAIL v%s\nhttps://github.com/mongrail/mongrail2\n\n",version);
}
