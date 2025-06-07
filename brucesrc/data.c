#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<regex.h>

void print_result(int value);

int main()
{

  regex_t reexp;
  int rexp_result, rexp_match;

  rexp_result = regcomp( &reexp, "[Aa][Bb][Cc]",0);

  if(rexp_result == 0) {
    printf("RegEx compiled successfully\n");
  }
  else {
    printf("compilation error\n");
  }

  rexp_match = regexec( &reexp, "ABg in a string", 0, NULL, 0);
    print_result(rexp_match);
}

void print_result(int value)
{

    // If pattern found
    if (value == 0) {
        printf("Pattern found.\n");
    }

    // If pattern not found
    else if (value == REG_NOMATCH) {
        printf("Pattern not found.\n");
    }

    // If error occurred during Pattern
    // matching
    else {
        printf("An error occurred.\n");
    }
}
