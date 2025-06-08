#include <stdlib.h>
#include <string.h>
#include <stdio.h> 


int main() {
  char str[] = "This is a  sample\tstring";
    char *token;
    char *words[100]; // Array to store words, adjust size as needed
    int i = 0;

    // Create a copy of the string
    char *str_copy = strdup(str);

    token = strtok(str_copy, " \t"); // Tokenize based on spaces

    while (token != NULL) {
        words[i] = strdup(token); // Allocate memory for the token
        i++;
        token = strtok(NULL, " \t");
    }
    // Print the words
    for (int j = 0; j < i; j++) {
        printf("%s\n", words[j]);
        free(words[j]); // Free allocated memory
    }
    free(str_copy); // Free the copy of the string

    return 0;
}
