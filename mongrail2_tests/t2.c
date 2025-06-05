#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define MAXLINESZ 100000
#define MAXNAMESZ 1000

int findstring(char strarray[1000][MAXNAMESZ], char *strtarget, int len);

int main()
{
  char str_array[1000][MAXNAMESZ];
  char *str1 = "chrom.1";
  char *str2 = "chrom.3";
  strncpy(str_array[0],"chrom.1",MAXNAMESZ);
  strncpy(str_array[1],"chrom.2",MAXNAMESZ);
  strncpy(str_array[2],"chrom.4",MAXNAMESZ);
  printf("Check T: %d",findstring(str_array,str1,3));
  printf("Check F: %d",findstring(str_array,str2,3));
}

int findstring(char strarray[1000][MAXNAMESZ], char *strtarget, int len)
{
  for(int i=0; i<len; i++)
    if(strcmp(strarray[i],strtarget)==0)
      return(1);
  return(0);
}
