#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct m_haplotype {
  unsigned int m_haplo;
  unsigned int count;
  struct m_haplotype* next;
};

/* convert unsigned int to binary string */

void uint_to_binary_string(unsigned int n, char *out, int out_size) {
    // Fill with zeros for fixed-width output
    int i, bit;
    if (out_size < 2) {
        // At minimum, need space for at least '0' and null terminator!
        if (out_size > 0) out[0] = '\0';
        return;
    }
    out[out_size - 1] = '\0'; // Null-terminate string
    for (i = out_size - 2; i >= 0; i--) {
        bit = n & 1;
        out[i] = bit ? '1' : '0';
        n >>= 1;
    }
    // Output will have leading zeros; if you want variable-width, find first '1'
}


int main()
{
  int maxloci = 12;
  char binary_string[maxloci+1];
  int maxhaps = pow(2,maxloci);
  struct m_haplotype* m_haps = NULL;
  m_haps = malloc(maxhaps*sizeof(struct m_haplotype*));
  unsigned int* hapcountA;
  unsigned int* hapcountB;
  int noHapsA=3;
  int noHapsB=3;
  unsigned int* haplistA;
  unsigned int* haplistB;
  hapcountA = calloc(maxhaps,sizeof(unsigned int));
  hapcountB = calloc(maxhaps,sizeof(unsigned int));
  haplistA = malloc(maxhaps*sizeof(unsigned int));
  haplistB = malloc(maxhaps*sizeof(unsigned int));
  hapcountA[6] = 5;
  hapcountA[4] = 10;
  hapcountA[2] = 5;
  hapcountB[0] = 10;
  hapcountB[4] = 5;
  hapcountB[7] = 5;
  haplistA[0] = 2;
  haplistA[1] = 4;
  haplistA[2] = 6;
  haplistB[0] = 0;
  haplistB[1] = 4;
  haplistB[2] = 7;
  
  for(int i=0; i<noHapsA; i++)
    {
      uint_to_binary_string(haplistA[i], binary_string, sizeof(binary_string));
      printf("haplotype popA: %s counts: %u\n",binary_string,hapcountA[haplistA[i]]);
    }
    for(int i=0; i<noHapsB; i++)
    {
      uint_to_binary_string(haplistB[i], binary_string, sizeof(binary_string));
      printf("haplotype popB: %s counts: %u\n",binary_string,hapcountB[haplistB[i]]);
    }
    uint_to_binary_string(4095, binary_string, sizeof(binary_string));
    printf("4095 in binary: %s\n",binary_string);

}
