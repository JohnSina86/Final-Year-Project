#include <stdio.h>
#include <math.h>

#define VERTICES 5
#define ROWS 3

// Function to multiply 3x3 matrix with 3xN matrix
void matrix_multiply(double matrix[ROWS][ROWS], double polygon[ROWS][VERTICES], double result[ROWS][VERTICES]) {
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < VERTICES; j++) {
            result[i][j] = 0;
            for (int k = 0; k < ROWS; k++) {
                result[i][j] += matrix[i][k] * polygon[k][j];
            }
        }
    }
}

// Function to print polygon vertices
void print_polygon(const char* name, double polygon[ROWS][VERTICES]) {
    printf("\n%s:\n", name);
    for (int i = 0; i < 2; i++) {  // Only print x and y coordinates
        printf("%c: ", i == 0 ? 'x' : 'y');
        for (int j = 0; j < VERTICES; j++) {
            printf("%8.2f ", polygon[i][j]);
        }
        printf("\n");
    }
}

int main() {
    // Define original polygon (square)
    double polygon[ROWS][VERTICES] = {
        {0, 10, 10, 0, 0},    // x coordinates
        {0, 0, 10, 10, 0},    // y coordinates
        {1, 1, 1, 1, 1}       // homogeneous coordinate
    };

    // Calculate centroid
    double centre_x = (polygon[0][0] + polygon[0][1] + polygon[0][2] + polygon[0][3]) / 4.0;
    double centre_y = (polygon[1][0] + polygon[1][1] + polygon[1][2] + polygon[1][3]) / 4.0;
    printf("Centroid: (%.2f, %.2f)\n", centre_x, centre_y);

    // Rotation matrix (135 degrees)
    double theta = 3.0 * M_PI / 4.0;
    double rot_matrix[ROWS][ROWS] = {
        {cos(theta), -sin(theta), 0},
        {sin(theta), cos(theta), 0},
        {0, 0, 1}
    };
    double rot_polygon[ROWS][VERTICES];
    matrix_multiply(rot_matrix, polygon, rot_polygon);

    // Translation matrix (tx=2, ty=4)
    double trans_matrix[ROWS][ROWS] = {
        {1, 0, 2},
        {0, 1, 4},
        {0, 0, 1}
    };
    double trans_polygon[ROWS][VERTICES];
    matrix_multiply(trans_matrix, polygon, trans_polygon);

    // Reflection matrix (reflect about y-axis)
    double ref_matrix[ROWS][ROWS] = {
        {-1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    double ref_polygon[ROWS][VERTICES];
    matrix_multiply(ref_matrix, polygon, ref_polygon);

    // Scaling matrix (scale by 0.5)
    double scaling_matrix[ROWS][ROWS] = {
        {0.5, 0, 0},
        {0, 0.5, 0},
        {0, 0, 1}
    };
    double scaled_polygon[ROWS][VERTICES];
    matrix_multiply(scaling_matrix, polygon, scaled_polygon);

    // Shear matrix (shear factor 0.5 in x-direction)
    double shear_matrix[ROWS][ROWS] = {
        {1, 0.5, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    double shear_polygon[ROWS][VERTICES];
    matrix_multiply(shear_matrix, polygon, shear_polygon);

    // Combined: scaling then translation
    double temp[ROWS][ROWS];
    double combined_matrix[ROWS][ROWS];
    
    // Multiply trans_matrix * scaling_matrix
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < ROWS; j++) {
            combined_matrix[i][j] = 0;
            for (int k = 0; k < ROWS; k++) {
                combined_matrix[i][j] += trans_matrix[i][k] * scaling_matrix[k][j];
            }
        }
    }
    
    double scaled_then_translated[ROWS][VERTICES];
    matrix_multiply(combined_matrix, polygon, scaled_then_translated);

    // Print all results
    print_polygon("Original Polygon", polygon);
    print_polygon("Rotated Polygon", rot_polygon);
    print_polygon("Translated Polygon", trans_polygon);
    print_polygon("Reflected Polygon", ref_polygon);
    print_polygon("Scaled Polygon", scaled_polygon);
    print_polygon("Sheared Polygon", shear_polygon);
    print_polygon("Scaled then Translated", scaled_then_translated);

    return 0;
}
