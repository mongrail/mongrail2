#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>


void prn_binary(unsigned int x, int noloci)
{
for (int i = 0; i < 32 && i < noloci; ++i) {
  if (x >> i & 0x1) putchar('1');
  else putchar('0');
}
}

unsigned int intpow(unsigned int base, unsigned int exponent)
{
  unsigned int result = 1;
  for (unsigned int i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

double prPopSwitch(double r, unsigned long d)
{
  return 1.0 - cosh(r*d)*exp(-r*d); 
}

double prQ(unsigned int anc, int noloci, double rec_rate_bp, unsigned long* position)
{
  double t1 = 0.0;
  double t2 = 0.0;
  double Q = 0.0;
  double p = 0.0;
  double pprod = 1.0;
  unsigned int mask1 = intpow(2,noloci) - 2;
  p = prPopSwitch(rec_rate_bp,position[0]);
  if((anc & mask1) == 0)
    {
      t1 = 0.5*(1.0 - p);
      t2 = 0.5*p;
    }
  else
    {
      t1 = 0.5*p;
      t2 = 0.5*(1.0 - p);
    }
  int curr_anc = 0;
  int prev_anc = ((anc >> 0)&1);
  for(int i=1; i<noloci; i++)
    {
      curr_anc = ((anc >> i)&1);
      p = prPopSwitch(rec_rate_bp,position[i]-position[i-1]); 
      if(curr_anc != prev_anc)
	pprod *= p;
      else
	pprod *= 1.0 - p;
      prev_anc = ((anc >> i)&1);
    }
  Q = pprod*(t1 + t2);
  return Q;
}

int main()
{
  static unsigned long position[] = {100,200,300};
  double prob_sum=0.0;
  double curr_pr;
  for(unsigned int i=0; i<8; i++)
    {
      printf("prQ(");
      curr_pr = prQ(i,3,0.01,position);
      prob_sum += curr_pr;
      prn_binary(i,3);
      printf("): %f\n",curr_pr);
    }
  printf("Total pr: %f\n",prob_sum);
}
