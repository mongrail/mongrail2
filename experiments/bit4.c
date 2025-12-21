#include <stdio.h>
#include <limits.h> // For CHAR_BIT

int main() {
    unsigned int num = 12; // Example: 1100 in binary
    int num_bits = sizeof(unsigned int) * CHAR_BIT; // Total bits in unsigned int

    printf("Bits of %u:\n", num);
    for (int i = 0; i < 14; i++) {
        // Check if the current bit is set (1) or not (0)
        if ((num >> i) & 1) {
            printf("Bit %d: 1\n", i);
        } else {
            printf("Bit %d: 0\n", i);
        }
    }
    return 0;
}
