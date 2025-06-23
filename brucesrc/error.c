#include <stdio.h>
#include <stdlib.h>
#include "mongrail.h"

void err_line_n(char pAfn[MAXFILENMSZ], char pBfn[MAXFILENMSZ], char hybfn[MAXFILENMSZ])
{
  int noflA = countfilelines(pAfn);
  if(noflA != countfilelines(pBfn))
    {
      fprintf(stderr,"Error: number of markers in files %s and %s must be equal.\n",pAfn,pBfn);
      exit(1);
    }
  if(noflA != countfilelines(hybfn))
    {
      fprintf(stderr,"Error: number of markers in files %s and %s must be equal.\n",pAfn,hybfn);
      exit(1);
    }
}
