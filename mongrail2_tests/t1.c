

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define MAXLINESZ 100000
#define MAXNAMESZ 1000
#define MAXCHRNUM 1000
#define MAXLOCI 20

int findstring (char strarray[MAXCHRNUM][MAXNAMESZ], char *strtarget, int len);
int countfilelines (char filename[]);
void file_to_array (char *target_array[MAXNAMESZ], char filename[]);
int mk_chrom_list (char* file_array[MAXNAMESZ], char names[MAXCHRNUM][MAXNAMESZ], int* no_loci, int nolines);

int main(int argc, char *argv[])
{
  char *ln_buffer;
  int no_file_lines=0;
  int noChrom=0;
  char ch;
  char** raw_data_popA = NULL; 
  char** raw_data_popB = NULL;
  char** raw_data_hybrids = NULL;
  char chr_names[MAXCHRNUM][MAXNAMESZ];
  int no_loci[MAXCHRNUM] = {0};
  FILE *fp;
  long count = 0;
 
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
  /* for(int i=0; i<no_file_lines; i++) */
  /*   { */
  /*     printf("popA: %s\n",raw_data_popA[i]); */
  /*   } */
  file_to_array(raw_data_popB,argv[2]);
  /* for(int i=0; i<no_file_lines; i++) */
  /*   { */
  /*     printf("popB: %s\n",raw_data_popB[i]); */
  /*   } */
  file_to_array(raw_data_hybrids,argv[3]);
  /* for(int i=0; i<no_file_lines; i++) */
  /*   { */
  /*     printf("hybrids: %s\n",raw_data_hybrids[i]); */
  /*   } */
  noChrom = mk_chrom_list (raw_data_popA,chr_names,no_loci,no_file_lines);
  printf("Chromosome/Scaffold => No. SNP Positions\n");
  printf("------------------------------------------\n");
  for(int i=0; i<noChrom; i++)
    printf("%s => %d\n",chr_names[i],no_loci[i]);
  /* ln_buffer = malloc(MAXLINESZ); */
  /* no_file_lines = countfilelines(argv[1]); */
  /* raw_data_popA = (char**)malloc(no_file_lines*sizeof(char *)); */

  /* if (raw_data_popA == NULL) { */
  /*   perror("Memory allocation error"); */
  /*   return 1; */
  /* } */
  
  /* if ((fp = fopen(argv[1], "r")) == NULL) */
  /*   { */
  /*     printf("Can't open %s\n", argv[1]); */
  /*     exit(1); */
  /*   } */

  /* int linepos = 0; */

  /* while (fgets(ln_buffer, MAXLINESZ, fp) != NULL) { */
  /*   size_t len = strlen(ln_buffer); */
  /*   if (len > 0 && ln_buffer[len - 1] == '\n') { */
  /*     ln_buffer[len - 1] = '\0'; */
  /*   } */
  /*   raw_data_popA[linepos] = strdup(ln_buffer); */
  /*   /\* add memory allocation error check here *\/ */
  /*   linepos++; */
  /* } */
  /* fclose(fp); */

  
  /* if ((fp = fopen(argv[1], "r")) == NULL) */
  /*   { */
  /*     printf("Can't open %s\n", argv[1]); */
  /*     exit(1); */
  /*   } */

  
  /* int chrpos = 0; */
  /* int strpos = 0; */
  
  /* while (fgets(ln_buffer, MAXLINESZ, fp)) { */
  /*   assert(strlen(ln_buffer)!= MAXLINESZ-1); /\* check that line read from input file does not exceed MAXLINESZ characters *\/ */
  /*   token = strtok(ln_buffer,delim); */
  /*   if(chrpos > 0) */
  /*     { */
  /* 	if(findstring(chr_names,token,chrpos)==0) */
  /* 	  { */
  /* 	    strncpy(chr_names[chrpos],token,MAXNAMESZ); */
  /* 	    chrpos++; */
  /* 	  } */
  /* 	no_loci[chrpos] = no_loci[chrpos] + 1;  */
  /*     } */
  /*   else */
  /*     { */
  /* 	strncpy(chr_names[chrpos],token,MAXNAMESZ); */
  /* 	no_loci[chrpos] = no_loci[chrpos] + 1;; */
  /* 	chrpos++; */
  /*     } */
  /*   while(token != NULL) { */
  /*     printf("T: %s\n",token); */
  /*     token = strtok(NULL,delim); */
  /*   } */
  /*   // Print each line to the standard output. */
  /*   //  printf("%s ln:%ld", ln_buffer, strlen(ln_buffer)); */
  /* } */
  /* noChrom = chrpos; */
  
  /* /\* while ((ch = getc(fp)) != EOF) *\/ */
  /* /\*   { *\/ */
  /* /\*     putc(ch,stdout); *\/ */
  /* /\*     count++; *\/ */
  /* /\*   } *\/ */

  /* fclose(fp); */
  /* //  printf("File %s has %d lines\n", argv[1], chrpos); */
  /* printf("No Chrom: %d\n",noChrom); */
  /* for(int i=0; i<noChrom; i++) { */
  /*   printf("Chr: %s Noloci: %d\n",chr_names[i],no_loci[i]); */
  /* } */
  return 0;
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

void file_to_array (char *target_array[MAXNAMESZ], char filename[])
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

int mk_chrom_list (char* file_array[MAXNAMESZ], char names[MAXCHRNUM][MAXNAMESZ], int* no_loci, int nolines)
{
  const char delim[] = ":";
  char *token;
  int chrpos=0;
  for(int i=0; i<nolines; i++)
    {
      token = strtok(file_array[i],delim);
      if(chrpos > 0)
	{
	  if(findstring(names,token,chrpos)==0)
	    {
	      strncpy(names[chrpos],token,MAXNAMESZ);
	      chrpos++;
	    }
	  no_loci[chrpos] = no_loci[chrpos] + 1;
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
