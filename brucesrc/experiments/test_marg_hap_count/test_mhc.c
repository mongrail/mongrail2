#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

// function to reverse the string
void reverse(char *bin, int left, int right) {
    while (left < right) {
        char temp = bin[left];
        bin[left] = bin[right];
        bin[right] = temp;
        left++;
        right--;
    }
}

char* decToBinary(unsigned int n, int noloci) {
    int index = 0;
    char* bin = (char*) malloc(32 * sizeof(char));
  
    while (n > 0) {
        int bit = n % 2;
        bin[index++] = '0' + bit;
        n /= 2;
    }
    while(index < noloci)
      {
	bin[index] = '0';
	index++;
      }
    bin[index] = '\0';
  
	// Reverse the binary string
	reverse(bin, 0, index-1);
  	return bin;
}

unsigned int intpow(unsigned int base, unsigned int exponent)
{
  unsigned int result = 1;
  for (unsigned int i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

int match_sub_hap(unsigned int hyb_hap, unsigned int pop_hap, unsigned int ancvec)
{
  if((hyb_hap & ancvec)==(pop_hap & ancvec))
    return 1;
  else
    return 0;
}

int get_marginal_hap_counts(unsigned int hyb_hap, unsigned int ancvec, int** popY_hap_counts, unsigned int** haplist, int* no_haps, int chrom)
{ 
  int mhcount=0;
  for(int i=0; i<no_haps[chrom]; i++)
    if(match_sub_hap(hyb_hap,haplist[chrom][i],ancvec))
      mhcount += popY_hap_counts[chrom][haplist[chrom][i]];
  return mhcount;
}

unsigned int get_anc_complement(unsigned int ancvec,int noloci)
{
  unsigned int cav = ~ancvec; 
  return cav - intpow(2,32) + intpow(2,noloci);;
}

int main()
{
  int noloci=3;
  unsigned int maxhaps;
  unsigned int nohaps = 3;
  int nochrom = 1;
  int ** pophapCounts;
  unsigned int ** hlist;
  char * bin_number;
  unsigned int ancvec=1;  
  unsigned int hybrid=5;
  int nhapvec[] = {3,3}; 
  maxhaps=intpow(2,noloci);
  pophapCounts = (int **)malloc(nochrom*sizeof(int *));
  for(int i=0; i<nochrom; i++)
    pophapCounts[i] = (int *)calloc(4,maxhaps*sizeof(int));
  hlist = (unsigned int **)malloc(nochrom*sizeof(unsigned int *));
  for(int i=0; i<nochrom; i++)
    hlist[i] = (unsigned int *)calloc(4,nohaps*sizeof(unsigned int));
  
  pophapCounts[0][0] = 5;
  pophapCounts[0][1] = 1;
  pophapCounts[0][7] = 4;
  hlist[0][0]=0;
  hlist[0][1]=7;
  hlist[0][2]=1;

  unsigned int comp_ancvec = get_anc_complement(ancvec,noloci);  
  printf("pop_counts: ");
  for(unsigned int i=0; i<maxhaps; i++)
    {
      bin_number = decToBinary(i,noloci); 
      printf("%s: %d ", bin_number,pophapCounts[0][i]);
      free(bin_number);
    }
  printf("\n");

  printf("hap_list: ");
  for(unsigned int i=0; i<nohaps; i++)
    {
      bin_number = decToBinary(hlist[0][i],noloci); 
      printf("%s ", bin_number);
      free(bin_number);
    }
  printf("\n");
  printf("hvbrid: %s\n",decToBinary(hybrid,noloci)); 
  printf("ancvec: %s marginal counts: %d\n",decToBinary(ancvec,noloci),get_marginal_hap_counts(hybrid,ancvec,pophapCounts,hlist,nhapvec,0));
  printf("comp_ancvec: %s marginal counts: %d\n",decToBinary(comp_ancvec,noloci),get_marginal_hap_counts(hybrid,comp_ancvec,pophapCounts,hlist,nhapvec,0));  
}

