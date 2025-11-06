#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 1000

// Function to calculate least squares using direct formula
double least_squares_direct(double x[], double y[], int N) {
    double sum_xy = 0.0;
    double sum_xx = 0.0;
    
    for (int ct = 0; ct < N; ct++) {
        sum_xy += x[ct] * y[ct];
        sum_xx += x[ct] * x[ct];
    }
    
    return sum_xy / sum_xx;
}

// Function to calculate least squares using linear algebra method
double least_squares_linalg(double x[], double y[], int N) {
    // Calculate A'*A (which is sum of x^2)
    double AtA = 0.0;
    for (int ct = 0; ct < N; ct++) {
        AtA += x[ct] * x[ct];
    }
    
    // Calculate A'*y (which is sum of x*y)
    double Aty = 0.0;
    for (int ct = 0; ct < N; ct++) {
        Aty += x[ct] * y[ct];
    }
    
    // Calculate (A'*A)^-1 * A'*y
    return Aty / AtA;
}

// Function to calculate sum of squared errors
double calculate_error(double x[], double y[], int N, double b) {
    double error = 0.0;
    
    for (int ct = 0; ct < N; ct++) {
        double predicted_y = b * x[ct];
        double diff = y[ct] - predicted_y;
        error += diff * diff;
    }
    
    return error;
}

// Function to save data for plotting
void save_to_csv(const char* filename, double x[], double y[], double pred_y[], 
                 double alt_pred_y[], double not_good_y[], int N) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return;
    }
    
    fprintf(fp, "x,y,predicted_y,alt_predicted_y,not_as_good_y\n");
    for (int i = 0; i < N; i++) {
        fprintf(fp, "%.6f,%.6f,%.6f,%.6f,%.6f\n", 
                x[i], y[i], pred_y[i], alt_pred_y[i], not_good_y[i]);
    }
    
    fclose(fp);
    printf("\nData saved to %s for plotting\n", filename);
}

int main() {
    double x[MAX_POINTS], y[MAX_POINTS];

    // Load data from file
    FILE *fp = fopen("data.res", "r");
    if (fp == NULL) {
        printf("Error opening data.res\n");
        return 1;
    }

    int N = 0;
    while (fscanf(fp, "%lf %lf", &x[N], &y[N]) == 2) {
        N++;
        if (N >= MAX_POINTS) break;
    }
    fclose(fp);
    
    printf("=== Simple Least Squares Regression: y = b*x ===\n");
    printf("Number of data points: %d\n\n", N);
    
    // Method 1: Direct formula
    printf("Method 1: Direct Formula\n");
    double b_tilde = least_squares_direct(x, y, N);
    printf("b_tilde = %.10f\n\n", b_tilde);
    
    // Method 2: Linear algebra approach
    printf("Method 2: Linear Algebra (A'*A)^-1 * A'*y\n");
    double alt_b = least_squares_linalg(x, y, N);
    printf("alt_b = %.10f\n\n", alt_b);
    
    // Verify both methods give same result
    printf("Difference between methods: %.10e\n\n", fabs(b_tilde - alt_b));
    
    // Define a suboptimal slope for comparison
    double not_as_good_b = b_tilde * 1.05;  // 5% higher than optimal
    printf("Suboptimal slope (5%% higher): %.10f\n\n", not_as_good_b);
    
    // Calculate predictions and errors
    double predicted_y[MAX_POINTS];
    double alt_predicted_y[MAX_POINTS];
    double not_as_good_y[MAX_POINTS];
    
    for (int ct = 0; ct < N; ct++) {
        predicted_y[ct] = b_tilde * x[ct];
        alt_predicted_y[ct] = alt_b * x[ct];
        not_as_good_y[ct] = not_as_good_b * x[ct];
    }
    
    // Calculate sum of squared errors
    double error = calculate_error(x, y, N, b_tilde);
    double not_as_good_error = calculate_error(x, y, N, not_as_good_b);
    
    printf("=== Error Analysis ===\n");
    printf("Sum of squared errors (optimal b_tilde):  %.10f\n", error);
    printf("Sum of squared errors (suboptimal b):     %.10f\n", not_as_good_error);
    printf("Error difference:                         %.10f\n\n", not_as_good_error - error);
    
    // Display data table
    printf("=== Data Points and Predictions ===\n");
    printf("  i     x         y      pred_y   alt_pred  not_good\n");
    printf("----  ------  ------  --------  --------  --------\n");
    for (int ct = 0; ct < N; ct++) {
        printf("%3d  %6.2f  %6.2f  %8.4f  %8.4f  %8.4f\n", 
               ct+1, x[ct], y[ct], predicted_y[ct], alt_predicted_y[ct], not_as_good_y[ct]);
    }
    
    // Save data to CSV for external plotting
    save_to_csv("least_squares_data.csv", x, y, predicted_y, alt_predicted_y, not_as_good_y, N);
    
    printf("\nTo visualize the results, use the CSV file with plotting software.\n");
    
    return 0;
}
