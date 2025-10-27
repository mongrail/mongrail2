#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mongrail.h"

extern int verbose;
  
void readDataFiles(char popAfileNm[], char popBfileNm[], char hybridfileNm[],
		   int* noSamplesPopA, int* noSamplesPopB, int* noSamplesPophybrid,
		   int* noChrom, char chr_names[MAXCHRNUM][MAXNAMESZ], int no_loci[MAXCHRNUM],
		   unsigned long** marker_positions, unsigned int** popA_haplotypes,
		   unsigned int** popB_haplotypes,unsigned int** hybrid_haplotypes)
{

  char** raw_data_popA = NULL; 
  char** raw_data_popB = NULL;
  char** raw_data_hybrids = NULL;
  struct genotype*** popA_genotypes = NULL;
  struct genotype*** popB_genotypes = NULL;
  struct genotype*** pophybrid_genotypes = NULL;
  int** popA_noIndivs = NULL;
  int** popB_noIndivs = NULL;
  int** pophybrid_noIndivs = NULL;

  /* total number of lines in popA.GT */ 
  int no_file_lines = countfilelines(popAfileNm);

  /* error checking: # lines equals # loci. Must be equal in popA.GT, popB.GT and hybrids.GT */
  err_line_n(popAfileNm, popBfileNm, hybridfileNm);

  int linepos=0;

  /* read popA data */
  raw_data_popA = (char**)malloc(no_file_lines*sizeof(char *));
  file_to_array(raw_data_popA, popAfileNm);
 
  
  /* read popB data */
  raw_data_popB = (char**)malloc(no_file_lines*sizeof(char *));
  file_to_array(raw_data_popB, popBfileNm);
 

  /* read hybrid data */
  raw_data_hybrids = (char**)malloc(no_file_lines*sizeof(char *));
  file_to_array(raw_data_hybrids, hybridfileNm);
 

  /* make chromosome list and get number of loci per chromosome */
  *noChrom = mk_chrom_list(raw_data_popA, chr_names, no_loci, no_file_lines);

  /* get marker positions for each chromosome */
  marker_positions = (unsigned long**)malloc(*noChrom*sizeof(unsigned long *));
  for(int i=0; i<*noChrom; i++)
    marker_positions[i] = (unsigned long*)malloc(no_loci[i]*sizeof(unsigned long));
  
  get_positions(*noChrom, no_loci, raw_data_popA, marker_positions);

  /* get individual genotypes for each population sample and counts of individuals */
  
  popA_genotypes = (struct genotype***)malloc(*noChrom*sizeof(struct genotype *));
  for(int i=0; i<*noChrom; i++)
    {
      popA_genotypes[i] = (struct genotype**)malloc(no_loci[i]*sizeof(struct genotype *));
      for(int j=0; j<no_loci[i]; j++)
	popA_genotypes[i][j] = (struct genotype*)malloc(MAXINDIV*sizeof(struct genotype ));
	
    }

  popB_genotypes = (struct genotype***)malloc(*noChrom*sizeof(struct genotype *));
  for(int i=0; i<*noChrom; i++)
    {
      popB_genotypes[i] = (struct genotype**)malloc(no_loci[i]*sizeof(struct genotype *));
      for(int j=0; j<no_loci[i]; j++)
	popB_genotypes[i][j] = (struct genotype*)malloc(MAXINDIV*sizeof(struct genotype ));
	
    }
  
  pophybrid_genotypes = (struct genotype***)malloc(*noChrom*sizeof(struct genotype *));
  for(int i=0; i<*noChrom; i++)
    {
      pophybrid_genotypes[i] = (struct genotype**)malloc(no_loci[i]*sizeof(struct genotype *));
      for(int j=0; j<no_loci[i]; j++)
	pophybrid_genotypes[i][j] = (struct genotype*)malloc(MAXINDIV*sizeof(struct genotype ));
	
    }

  /* get number of individuals sampled in pop A and B and putative hybrids */
  /* collect individual genotypes for each chromosome and locus */
  
  popA_noIndivs = (int**)malloc(*noChrom*sizeof(int *));
  for(int i=0; i<*noChrom; i++)
    popA_noIndivs[i] = (int*)malloc(no_loci[i]*sizeof(int));

  popB_noIndivs = (int**)malloc(*noChrom*sizeof(int *));
  for(int i=0; i<*noChrom; i++)
    popB_noIndivs[i] = (int*)malloc(no_loci[i]*sizeof(int));

   pophybrid_noIndivs = (int**)malloc(*noChrom*sizeof(int *));
  for(int i=0; i<*noChrom; i++)
    pophybrid_noIndivs[i] = (int*)malloc(no_loci[i]*sizeof(int));

  linepos=0;
  for(int i=0; i<*noChrom; i++)
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
  *noSamplesPopA = popA_noIndivs[0][0];

 linepos=0;
  for(int i=0; i<*noChrom; i++)
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
  *noSamplesPopB = popB_noIndivs[0][0];

linepos=0;
  for(int i=0; i<*noChrom; i++)
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
  *noSamplesPophybrid = pophybrid_noIndivs[0][0];

  /* get sampled population haplotypes as unsigned ints */
  /* number of haplotypes = 2 * number of individuals sampled */
  
  popA_haplotypes = (unsigned int **)malloc(*noChrom*sizeof(unsigned int *));
  for(int i=0; i<*noChrom; i++)
    popA_haplotypes[i] = (unsigned int *)malloc(*noSamplesPopA*2*sizeof(unsigned int));
  get_haplotypes(popA_haplotypes,popA_genotypes,*noChrom,*noSamplesPopA,no_loci);

  popB_haplotypes = (unsigned int **)malloc(*noChrom*sizeof(unsigned int *));
  for(int i=0; i<*noChrom; i++)
    popB_haplotypes[i] = (unsigned int *)malloc(*noSamplesPopB*2*sizeof(unsigned int));
  get_haplotypes(popB_haplotypes,popB_genotypes,*noChrom,*noSamplesPopB,no_loci);

  hybrid_haplotypes = (unsigned int **)malloc(*noChrom*sizeof(unsigned int *));
  for(int i=0; i<*noChrom; i++)
    hybrid_haplotypes[i] = (unsigned int *)malloc(*noSamplesPophybrid*2*sizeof(unsigned int));
  get_haplotypes(hybrid_haplotypes,pophybrid_genotypes,*noChrom,*noSamplesPophybrid,no_loci);

  /* printing information to screen */

  pr_summary(popAfileNm, popBfileNm, hybridfileNm, *noChrom, no_loci, chr_names,
	     popA_noIndivs, popB_noIndivs, pophybrid_noIndivs, marker_positions);


  /* optional information printing to screen for debugging */
  
  if(verbose==1)
    {
      pr_haplotypes(*noChrom, chr_names, popA_haplotypes, popB_haplotypes, *noSamplesPopA, *noSamplesPopB);
    }
  if(verbose==2)
    {
      pr_datastats(*noChrom, chr_names, marker_positions, popA_noIndivs, popB_noIndivs, pophybrid_noIndivs,no_loci);
    }
  if(verbose==3)
    {
      pr_chr_SNP_stats(*noChrom, chr_names,no_loci, marker_positions);
      
    }

  /* free allocated memory */
  
  free(raw_data_popA);
  free(raw_data_popB);
  free(raw_data_hybrids);
 /*  free(popA_genotypes); */
/*   free(popB_genotypes); */
/*   free(pophybrid_genotypes); */
/*   free(popA_noIndivs); */
/*   free(popB_noIndivs); */
/*   free(pophybrid_noIndivs); */
 } 
  


int findstring(char strarray[MAXCHRNUM][MAXNAMESZ], char *strtarget, int len)
{
  for(int i=0; i<len; i++)
    if(strcmp(strarray[i],strtarget)==0)
      return(1);
  return(0);
}

int countfilelines(char filename[])
{
  FILE *file;
  char ch;
  int line_count=0;
  file = fopen(filename, "r");
  if (file == NULL) {
    printf("Error opening the file.\n");
    return -1;
  }
  while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            line_count++;
        }
    }
    fclose(file);
    return line_count;
}

void file_to_array (char *target_array[MAXLINESZ], char filename[])
{
  char buffer[MAXLINESZ];
  FILE *fp;
  int linepos = 0;
  if ((fp = fopen(filename, "r")) == NULL)
    {
      printf("Can't open %s\n", filename);
      exit(1);
    }
  while (fgets(buffer, MAXLINESZ, fp) != NULL) {
    size_t len = strlen(buffer);
    if (len > 0) {
      buffer[len - 1] = '\0';
    }
    target_array[linepos] = strdup(buffer);
    /* add memory allocation error check here */
    linepos++;
  }
  fclose(fp);
}

int mk_chrom_list (char* file_array[MAXLINESZ], char names[MAXCHRNUM][MAXNAMESZ], int* no_loci, int nolines)
{
  const char delim[] = ":";
  char *token;
  int chrpos=0;
  char buffer[MAXLINESZ];
  for(int i=0; i<nolines; i++)
    {
      strncpy(buffer,file_array[i],MAXLINESZ);
      token = strtok(buffer,delim);
      if(chrpos > 0)
	{
	  if(findstring(names,token,chrpos)==0)
	    {
	      strncpy(names[chrpos],token,MAXNAMESZ);
	      no_loci[chrpos] = no_loci[chrpos] + 1;
	      chrpos++;
	    }
	  else
	    no_loci[chrpos-1] = no_loci[chrpos-1] + 1;
	}
      else
	{
	  strncpy(names[chrpos],token,MAXNAMESZ);
	  no_loci[chrpos] = no_loci[chrpos] + 1;;
	  chrpos++;
	}
    }
  return chrpos;
}

int get_genotypes(char* gstring, struct genotype* chr_locus)
{
  char *token;
  const char delim[] = ": \t";
  int numindivs = 0;
  struct genotype geno[MAXINDIV];
  int g1, g2;
  char phase;

  char *str_copy = strdup(gstring);
  token = strtok(str_copy,delim); // Tokenize based on spaces and :
  token = strtok(NULL,delim);
  token = strtok(NULL,delim);
  while (token != NULL) {
    g1 = token[0] - '0';
    phase = token[1];
    g2 = token[2] - '0';
    geno[numindivs].g1 = g1;
    geno[numindivs].g2 = g2;
    geno[numindivs].phase = phase;
    numindivs++;
    token = strtok(NULL,delim);
  }
  for(int i=0; i<numindivs; i++)
    {
      chr_locus[i].g1 = geno[i].g1;
      chr_locus[i].g2 = geno[i].g2;
      chr_locus[i].phase = geno[i].phase;
      //      printf("%d %d%c%d \n",i,chr_locus[i].g1,chr_locus[i].phase,chr_locus[i].g2);
    }
  free(str_copy); // Free the copy of the string
  return numindivs;
}

void get_haplotypes(unsigned int** haplotypes, struct genotype*** gdata, int nChr, int noInd, int* nloci)
{
  int i,j,k;
  char hapstring_g1[30];
  char hapstring_g2[30];
  for(i=0; i<nChr; i++)
    for(j=0; j<noInd; j++)
      {
	for(k=0; k<nloci[i]; k++)
	  {
	    if(gdata[i][k][j].g1==0)
	      hapstring_g1[k] = '0';
	    else
	      hapstring_g1[k] = '1';
	    if(gdata[i][k][j].g2==0)
	      hapstring_g2[k] = '0';
	    else
	      hapstring_g2[k] = '1';
	  }
	hapstring_g1[nloci[i]]='\0';
	hapstring_g2[nloci[i]]='\0';
	haplotypes[i][j]=binaryToDecimal(hapstring_g1);
	haplotypes[i][j+noInd]=binaryToDecimal(hapstring_g2);
      }
}

  void get_positions(int noChr,int* no_loci, char** raw_data, unsigned long** positions)
  {
    int linepos=0;
    char* token;
    const char delim[] = ":";
    char buffer[MAXLINESZ];
    for(int i=0; i<noChr; i++)
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
	  strncpy(buffer,raw_data[linepos],MAXLINESZ);
	  token = strtok(buffer,delim);
	  token = strtok(NULL,delim);
	  positions[i][j] = strtol(token, NULL, 10);
	}
  }

void get_hap_counts(unsigned int** haplotypes, int** hap_counts, int noChr, int no_indiv)
{
  for(int i=0; i<noChr; i++)
    for(int j=0; j<2*no_indiv; j++)
      hap_counts[i][haplotypes[i][j]] = hap_counts[i][haplotypes[i][j]] + 1;
}

/* get number of polymorphic sites and number of haplotypes for hybrid indivs */
unsigned int count_haplotypes(unsigned int hap1, unsigned int hap2, int noloci)
{
  unsigned int hap1hap2;
  unsigned int no1bits=0;
  hap1hap2 = hap1^hap2; /* 1 if 1^0 or 0^1, 0 otherwise */
  for (int i = 0; i < 32 && i < noloci; ++i)
    if (hap1hap2 >> i & 0x1) no1bits++; /* if ith bit of hap1hap2 = 1 is polymorphic */
  if (no1bits == 0) return 2; /* if no polymorphic sites, only 2 haplotypes */
  return pow(2,no1bits); /* no haplotypes = 2^no1bits */
}

/* find all unique haplotypes compatible with genotypes */
void compatible_haps(unsigned int *hapvec, unsigned int genotype1, unsigned int genotype2)
{
  long int g1g2 = genotype1 ^ genotype2;
  int hapNo=1;
  unsigned int mask=1;
  for(int j=0; j<MAXLOCI; j++)
    {
      if((mask << j) & g1g2)
	{
	  hapNo = hapNo*2;
	  for(int k=(hapNo/2);k<hapNo;k++)
	    hapvec[k]=hapvec[k-(hapNo/2)];
	  for(int i=0;i<(hapNo/2);i++)
	    {
	      int allele=0; int pos=0;
	      while(pos<hapNo && allele < 2)
		{
		  if(hapvec[i] == hapvec[pos])
		    {
		      if(allele == 1)
			  hapvec[pos] = (mask << j) ^ hapvec[pos];
		      allele++;
		    }
		  pos++;
		}
	    }
	}
      else
	{
	  if((mask << j) & genotype1)
	    for(int i=0;i<hapNo;i++)
	      hapvec[i] = (mask << j) ^ hapvec[i];
	}
    }
}

/* arrange haplotypes into successive pairs (hap1,hap2) of diplotypes: even = hap1, odd = hap2 */
void sortDiplotypes(struct indiv sampleInd)
{
  int j;
  int i=0;
  if(sampleInd.numHaps > 2)
    {
      while(i < sampleInd.numHaps)
	{
	  j=1;
	  while(j < sampleInd.numHaps)
	    {
	      if((sampleInd.compHaps[i] ^ sampleInd.compHaps[j]) ==
		 (sampleInd.genotype1 ^ sampleInd.genotype2))
		{
		  unsigned int temp = sampleInd.compHaps[i+1];
		  sampleInd.compHaps[i+1] = sampleInd.compHaps[j];
		  sampleInd.compHaps[j] = temp;
		  i = i+2;
		  break;
		}
	      else
		j++;
	    }
	}
    }
  else
    {
    sampleInd.compHaps[0] = sampleInd.genotype1;
    sampleInd.compHaps[1] = sampleInd.genotype2;
    }
}


