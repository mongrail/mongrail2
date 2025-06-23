#include <math.h>
#include <string.h>
#include <stdio.h>

/* Returns the value ln[Γ(xx)] for xx > 0. */
float gammln(float xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
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

void prn_binary(unsigned int x, int noloci)
{
for (int i = 0; i < 32 && i < noloci; ++i) {
  if (x >> i & 0x1) putchar('1');
  else putchar('0');
}
}
