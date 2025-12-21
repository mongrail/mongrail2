#include <stdio.h>
#include <stdlib.h>
#include "mongrail.h"

void err_line_n(char pAfn[MAXFILENMSZ], char pBfn[MAXFILENMSZ], char hybfn[MAXFILENMSZ])
{
  int noflA = countfilelines(pAfn);
  if(noflA != countfilelines(pBfn))
    {
      fprintf(stderr,"Error: number of markers in files %s and %s must be equal.\n",pAfn,pBfn);
      exit(1);
    }
  if(noflA != countfilelines(hybfn))
    {
      fprintf(stderr,"Error: number of markers in files %s and %s must be equal.\n",pAfn,hybfn);
      exit(1);
    }
}


/* print output for debugging */

/* print haplotypes for debugging */ 

void pr_haplotypes(int noChr, char chr_nm[MAXCHRNUM][MAXNAMESZ], unsigned int** popA_haps, unsigned int** popB_haps, int nsPopA, int nsPopB)
{
  printf("\nPopulation\tChromosome\tHaplotypes\n");
  printf("-------------------------------------------\n");
  /* print haplotypes for population A */
  for(int i=0; i<noChr; i++)
    {
      printf("\nA\t%s ",chr_nm[i]);
      for(int j=0; j<2*nsPopA; j++)
	printf(" %u",popA_haps[i][j]);
    }
  /* print haplotypes for population B */
  for(int i=0; i<noChr; i++)
    {
      printf("\nB\t%s ",chr_nm[i]);
      for(int j=0; j<2*nsPopB; j++)
	printf(" %u",popB_haps[i][j]);
    }
}

void pr_hapcounts(int** hap_counts, char chr_nm[MAXCHRNUM][MAXNAMESZ], int noChr)
{
  for(int i=0; i<noChr; i++)
    {
      printf("\n%s: ",chr_nm[i]);
      for(int j=0; j<MAXHAPS; j++)
	if(hap_counts[i][j] > 0)
	  printf("%d: %d ",j,hap_counts[i][j]);
    }
}

void pr_datastats(int noChr, char chr_nm[MAXCHRNUM][MAXNAMESZ], unsigned long** positions, int** popA_noIndivs, int** popB_noIndivs, int** pophybrid_noIndivs, int no_loci[MAXCHRNUM])
{
  /* print chromosome names, marker positions and samples sizes to screen */
  /* population A */
  printf("Population A:\n");
  for(int i=0; i<noChr; i++)
    {
      printf("Chr:%s Position(NoInds)\n",chr_nm[i]);
      for(int j=0; j<no_loci[i]; j++)
	{
	  printf("%ld",positions[i][j]);
	  printf("(%d) ",popA_noIndivs[i][j]);
	}
      printf("\n");
    }
  /* population B */
  printf("Population B:\n");
  for(int i=0; i<noChr; i++)
    {
      printf("Chr:%s Position(NoInds)\n",chr_nm[i]);
      for(int j=0; j<no_loci[i]; j++)
	{
	  printf("%ld",positions[i][j]);
	  printf("(%d) ",popB_noIndivs[i][j]);
	}
      printf("\n");
    }
  /* hybrids */
  printf("Possible Hybrids:\n");
  for(int i=0; i<noChr; i++)
    {
      printf("Chr:%s Position(NoInds)\n",chr_nm[i]);
      for(int j=0; j<no_loci[i]; j++)
	{
	  printf("%ld",positions[i][j]);
	  printf("(%d) ",pophybrid_noIndivs[i][j]);
	}
      printf("\n");
    }

}

void pr_chr_SNP_stats(int noChr, char chr_nm[MAXCHRNUM][MAXNAMESZ], int no_loci[MAXCHRNUM], unsigned long** positions)
{
  fprintf(stderr, "SNP positions (bp):\n");
  for (int i = 0; i < noChr; i++) {
    fprintf(stderr, "  %s:", chr_nm[i]);
    for (int j = 0; j < no_loci[i]; j++)
      fprintf(stderr, " %ld", positions[i][j]);
    fprintf(stderr, "\n");
  }
}

void pr_summary(char popAfileNm[MAXFILENMSZ], char popBfileNm[MAXFILENMSZ], char hybridfileNm[MAXFILENMSZ], int noChr, int no_loci[MAXCHRNUM], char chr_nm[MAXCHRNUM][MAXNAMESZ], int** popA_noIndivs, int** popB_noIndivs, int** pophybrid_noIndivs,unsigned long** positions)
{
  fprintf(stderr, "Input files:\n");
  fprintf(stderr, "  Pop A:    %s\n", popAfileNm);
  fprintf(stderr, "  Pop B:    %s\n", popBfileNm);
  fprintf(stderr, "  Hybrids:  %s\n", hybridfileNm);
  fprintf(stderr, "\nData summary:\n");
  fprintf(stderr, "  %-16s %7s %8s %8s %8s\n", "Chromosome", "SNPs", "n_A", "n_B", "n_hyb");
  for (int i = 0; i < noChr; i++) {
    fprintf(stderr, "  %-16s %7d %8d %8d %8d\n",
	    chr_nm[i], no_loci[i], popA_noIndivs[i][0], popB_noIndivs[i][0], pophybrid_noIndivs[i][0]);
  }
  fprintf(stderr, "\n");
  pr_chr_SNP_stats(noChr, chr_nm, no_loci, positions);
  fprintf(stderr, "\n");
}

  /* debugging */
void pr_compatible_haps_hybrids(int noChrom, char chr_names[MAXCHRNUM][MAXNAMESZ], struct indiv** hybrid_indiv, int noSamplesPophybrid, int no_loci[MAXCHRNUM])
{
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
}
