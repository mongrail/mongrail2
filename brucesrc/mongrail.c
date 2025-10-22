#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mongrail.h"
 

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

  
  pr_progname(version);
  
  /* read raw GT data into arrays of lines as strings */
  
  /* total number of lines in popA.GT */ 


  
  no_file_lines = countfilelines(popAfileNm);

  /* error checking: # lines equals # loci. Must be equal in popA.GT, popB.GT and hybrids.GT */


  
  err_line_n(popAfileNm, popBfileNm, hybridfileNm);

  /* create arrays to contain lines read from each input file */

  raw_data_popA = (char**)malloc(no_file_lines*sizeof(char *));
  raw_data_popB = (char**)malloc(no_file_lines*sizeof(char *));
  raw_data_hybrids = (char**)malloc(no_file_lines*sizeof(char *));

  /* fill arrays with lines from each input file */

  file_to_array(raw_data_popA,argv[1]);
  file_to_array(raw_data_popB,argv[2]);
  file_to_array(raw_data_hybrids,argv[3]);

  /* parsing input data */
  
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

  /* get number of individuals sampled in pop A and B and putative hybrids */
  /* collect individual genotypes for each chromosome and locus */
  
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

  
  /* get sampled population haplotypes as unsigned ints */
  /* number of haplotypes = 2 * number of individuals sampled */
  
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
      for(int j=0; j < 2*popA_noIndivs[z][0]; j++)
	{
	  add_hap(popA_haplotypes[z][j], haplist, nohaps, z);
	} 
      for(int j=0; j < 2*popB_noIndivs[z][0]; j++)
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
  
  
  /* printing information to screen */

  pr_summary(popAfileNm, popBfileNm, hybridfileNm, noChrom, no_loci, chr_names, popA_noIndivs, popB_noIndivs, pophybrid_noIndivs, marker_positions);

  /* debugging likelihoods */
  init_likelihood_globals(noSamplesPopA,noSamplesPopB);
  printf("Indiv     Md                Ma                Mc\n");    
   for(int i=0; i<noSamplesPophybrid; i++)
     {
       printf("%5d      %.5f",i,lik_a_d(i,hybrid_indiv,popA_hap_counts,haplist,nohaps,noSamplesPopA,noChrom,MODEL_D));
       printf("         %.5f",lik_a_d(i,hybrid_indiv,popB_hap_counts,haplist,nohaps,noSamplesPopB,noChrom,MODEL_A));
       printf("         %.5f\n",lik_c(i,hybrid_indiv,popB_hap_counts,popA_hap_counts,haplist,nohaps,noSamplesPopB,noSamplesPopA,noChrom));
     }


   /*    printf("Indiv: %d Ma logL: %f Md logL: %f Mc logL: %f\n",i,lik_a_d(i,hybrid_indiv,popB_hap_counts,haplist,nohaps,noSamplesPopB,noChrom),lik_a_d(i,hybrid_indiv,popA_hap_counts,haplist,nohaps,noSamplesPopA,noChrom), */

  
  /* optional information printing to screen for debugging */
  
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
