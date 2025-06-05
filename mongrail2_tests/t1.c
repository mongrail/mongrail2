#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define MAXLINESZ 100000
#define MAXNAMESZ 1000

int findstring(char strarray[1000][MAXNAMESZ], char *strtarget, int len);

int main(int argc, char *argv[])
{
  char *aline;
  int noChrom=0;
  char ch;
  char str_array[1000][MAXNAMESZ];
  FILE *fp;
  long count = 0;
  const char delim[] = ":";
  char *token;
  
  aline = malloc(MAXLINESZ);
  
  if (argc != 2)
    {
      printf("Usage: %s filename\n", argv[0]);
      exit(1);
    }
  if ((fp = fopen(argv[1], "r")) == NULL)
    {
      printf("Can't open %s\n", argv[1]);
      exit(1);
    }
  int chrpos = 0;
  int strpos = 0;
  while (fgets(aline, MAXLINESZ, fp)) {
    assert(strlen(aline)!= MAXLINESZ-1); /* check that line read from input file does not exceed MAXLINESZ characters */
    token = strtok(aline,delim);
    if(chrpos > 0)
      {
	if(findstring(str_array,token,chrpos)==0) {
	  strncpy(str_array[chrpos],token,MAXNAMESZ);
	  chrpos++;
	}
      }
    else
      {
	strncpy(str_array[chrpos],token,MAXNAMESZ);
	chrpos++;
      }
    while(token != NULL) {
      printf("T: %s\n",token);
      token = strtok(NULL,delim);
    }
    // Print each line to the standard output.
    //  printf("%s ln:%ld", aline, strlen(aline));
  }
  noChrom = chrpos;
  
  /* while ((ch = getc(fp)) != EOF) */
  /*   { */
  /*     putc(ch,stdout); */
  /*     count++; */
  /*   } */

  fclose(fp);
  //  printf("File %s has %d lines\n", argv[1], chrpos);
  printf("No Chrom: %d\n",noChrom);
  for(int i=0; i<chrpos; i++)
    printf("Chr: %s\n",str_array[i]);
  return 0;
}

int findstring(char strarray[1000][MAXNAMESZ], char *strtarget, int len)
{
  for(int i=0; i<len; i++)
    if(strcmp(strarray[i],strtarget)==0)
      return(1);
  return(0);
}

