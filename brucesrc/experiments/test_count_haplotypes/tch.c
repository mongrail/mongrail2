#include<math.h>
#include<stdio.h>

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

int main()
{
	unsigned int hap1, hap2;
	int noloci;
	printf("Enter haplotype 1 (unsigned int): ");
	scanf("%u", &hap1);
	printf("Enter haplotype 2 (unsigned int): ");
	scanf("%u", &hap2);
	printf("Enter number of loci (int): ");
	scanf("%d", &noloci);
	unsigned int num_haplotypes = count_haplotypes(hap1, hap2, noloci);
	printf("Number of haplotypes for hybrid individual: %u\n", num_haplotypes);
	return 0;
}
