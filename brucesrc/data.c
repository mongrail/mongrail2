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

unsigned int binaryToDecimal(const char *binaryStr) {
    unsigned int decimal = 0;
    int len = strlen(binaryStr);
    for (int i = 0; i < len; i++) {
        if (binaryStr[i] == '1') {
            decimal += pow(2, len - 1 - i);
        }
    }
    return decimal;
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
    if (len > 0 && buffer[len - 1] == '\n') {
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
	hapstring_g1[nloci[i]]='\0';
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

void pr_chr_SNP_stats(int noChr, char chr_nm[MAXCHRNUM][MAXNAMESZ],int no_loci[MAXCHRNUM], unsigned long** positions )
 {
   printf("Chromosome/Scaffold => No. SNPs => Positions\n");
   printf("------------------------------------------\n");
   for(int i=0; i<noChr; i++)
     {
       printf("%s => %d =>",chr_nm[i],no_loci[i]);
       for(int j=0; j<no_loci[i]; j++)
	 printf(" %ld",positions[i][j]);
       printf("\n");
     }
 }
 
