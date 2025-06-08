#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define MAXINDIV 1000
#define MAXLINESZ 100000
#define MAXNAMESZ 1000
#define MAXCHRNUM 1000
#define MAXLOCI 20

struct genotype {
  int g1, g2;
  char phase;
};

int verbose=1;
int findstring (char strarray[MAXCHRNUM][MAXNAMESZ], char *strtarget, int len);
unsigned int binaryToDecimal (const char *binaryStr);
int countfilelines (char filename[]);
void file_to_array (char *target_array[MAXLINESZ], char filename[]);
int mk_chrom_list (char* file_array[MAXLINESZ], char names[MAXCHRNUM][MAXNAMESZ], int* no_loci, int nolines);
int get_genotypes (char* gstring, struct genotype* chr_locus);
void get_haplotypes (unsigned int** haplotypes, struct genotype*** gdata, int nChr, int noInd, int* nloci);
  
int main(int argc, char *argv[])
{
  const char delim[] = ":";
  char* token;
  int no_file_lines = 0;
  int noChrom = 0;
  int noSamplesPopA = 0;
  int noSamplesPopB = 0;
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
  unsigned long** marker_positions = NULL;
  char chr_names[MAXCHRNUM][MAXNAMESZ];
  int no_loci[MAXCHRNUM] = {0};
 
  if (argc != 4)
    {
      printf("Usage: %s popA popB hybrids\n", argv[0]);
      exit(1);
    }

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
  int linepos=0;
  char buffer[MAXLINESZ];
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
      strncpy(buffer,raw_data_popA[linepos],MAXLINESZ);
      token = strtok(buffer,delim);
      token = strtok(NULL,delim);
      marker_positions[i][j] = strtol(token, NULL, 10);
    }

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

  popA_noIndivs = (int**)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    popA_noIndivs[i] = (int*)malloc(no_loci[i]*sizeof(int));

  popB_noIndivs = (int**)malloc(noChrom*sizeof(int *));
  for(int i=0; i<noChrom; i++)
    popB_noIndivs[i] = (int*)malloc(no_loci[i]*sizeof(int));

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
  
  /* get population haplotypes in binary form */
  
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
      /* print haplotypes for population A */
      printf("\nPopulation\tChromosome\tHaplotypes\n");
      printf("-------------------------------------------\n");
      for(int i=0; i<noChrom; i++)
	{
	  printf("\nA\t%s ",chr_names[i]);
	  for(int j=0; j<2*noSamplesPopA; j++)
	    printf(" %u",popA_haplotypes[i][j]);
	}

      /* print haplotypes for population B */
      for(int i=0; i<noChrom; i++)
	{
	  printf("\nB\t%s ",chr_names[i]);
	  for(int j=0; j<2*noSamplesPopB; j++)
	    printf(" %u",popB_haplotypes[i][j]);
	}

    }
  
  if(verbose==2)
    {
      /* print chromosome names, marker positions and samples sizes to screen */
      /* population A */
      printf("Population A:\n");
      for(int i=0; i<noChrom; i++)
	{
	  printf("Chr:%s Position(NoInds)\n",chr_names[i]);
	  for(int j=0; j<no_loci[i]; j++)
	    {
	      printf("%ld",marker_positions[i][j]);
	      printf("(%d) ",popA_noIndivs[i][j]);
	    }
	  printf("\n");
	}
      /* population B */
      printf("Population B:\n");
      for(int i=0; i<noChrom; i++)
	{
	  printf("Chr:%s Position(NoInds)\n",chr_names[i]);
	  for(int j=0; j<no_loci[i]; j++)
	    {
	      printf("%ld",marker_positions[i][j]);
	      printf("(%d) ",popB_noIndivs[i][j]);
	    }
	  printf("\n");
	}
    }
  if(verbose==3)
    {
      /* print genotypes to screen */
      for(int i=0; i<noChrom; i++)
	{
	  printf("Chrom: %s \n",chr_names[i]);
	  printf("-------------------------\n");
	  printf("Position\tGenotypes\n");
	  for(int j=0; j<no_loci[i]; j++)
	    {
	      printf("%ld ",marker_positions[i][j]);
	      for(int k=0; k<=3; k++)
		printf(" %d%c%d ",popA_genotypes[i][j][k].g1,popA_genotypes[i][j][k].phase,popA_genotypes[i][j][k].g2);
	      printf("\n");
	    }
	  printf("\n");
	}
      
      
      printf("Chromosome/Scaffold => No. SNPs => Positions\n");
      printf("------------------------------------------\n");
      for(int i=0; i<noChrom; i++)
	{
	  printf("%s => %d =>",chr_names[i],no_loci[i]);
	  for(int j=0; j<no_loci[i]; j++)
	    printf(" %ld",marker_positions[i][j]);
	  printf("\n");
	}
    } 
      
      return 0;
    }

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
