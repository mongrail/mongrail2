#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mongrail.h"

  

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
	  /* if((mask << j) & genotype1) */
	  /*   for(int i=0;i<hapNo;i++) */
	  /*     hapvec[i] = (mask << j) ^ hapvec[i]; */
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
}


