#include <stdio.h>
#include <math.h> // For fabs()

// Function to perform Kahan summation
double kahanSum(double arr[], int n) {
    double sum = 0.0;
    double c = 0.0; // A running compensation for lost low-order bits

    for (int i = 0; i < n; i++) {
        double y = arr[i] - c;
        double t = sum + y;
        c = (t - sum) - y; // Calculate the lost low-order bits
        sum = t;
    }
    return sum;
}

// Function to perform naive summation for comparison
double naiveSum(double arr[], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += arr[i];
    }
    return sum;
}

int main() {
    // Example usage
    int n = 100000;
    double arr[n];

    // Populate array with values that highlight the issue (e.g., small values added to a large sum)
    for (int i = 0; i < n; i++) {
        arr[i] = 1.0 / (i + 1); // Reciprocals of integers
    }

    double kahan_result = kahanSum(arr, n);
    double naive_result = naiveSum(arr, n);

    printf("Kahan Summation Result: %.15f\n", kahan_result);
    printf("Naive Summation Result: %.15f\n", naive_result);
    printf("Difference: %.15e\n", fabs(kahan_result - naive_result));

    return 0;
}
