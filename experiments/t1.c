#include <stdio.h>
#include <math.h>

unsigned int intpow(unsigned int base, unsigned int exponent) {
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
  unsigned long position[16]={100,200,300,400,500,1000,2000,3000,4000,5000,8000,10000,11000,12000,13000,13500};
  int noloci=16;
  unsigned int anc=0;
  double rec_rate_bp = 0.0001;
  double sum = 0;
  for(unsigned int i=0; i<intpow(2,16); i++)
    {
      printf("anc: %u Q: %f\n",i,prQ(i,noloci,rec_rate_bp,position));
      sum+=prQ(i,noloci,rec_rate_bp,position);
    }
  printf("sum: %f\n",sum);
}
  
