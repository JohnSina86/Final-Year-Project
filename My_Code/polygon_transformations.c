#include <stdio.h>
#include <math.h>

#define VERTICES 14
#define ROWS 3

// Function to multiply 3x3 matrix with 3x3 matrix
void matrix_mult_3x3(double A[ROWS][ROWS], double B[ROWS][ROWS], double result[ROWS][ROWS]) {
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < ROWS; j++) {
            result[i][j] = 0;
            for (int k = 0; k < ROWS; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Function to multiply 3x3 matrix with 3xN polygon matrix
void matrix_mult_polygon(double matrix[ROWS][ROWS], double polygon[ROWS][VERTICES], double result[ROWS][VERTICES]) {
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
    printf("x: ");
    for (int j = 0; j < VERTICES; j++) {
        printf("%7.2f ", polygon[0][j]);
    }
    printf("\ny: ");
    for (int j = 0; j < VERTICES; j++) {
        printf("%7.2f ", polygon[1][j]);
    }
    printf("\n");
}

// Function to save polygon to CSV for plotting
void save_to_csv(const char* filename, double original[ROWS][VERTICES], double transformed[ROWS][VERTICES]) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return;
    }
    
    fprintf(fp, "orig_x,orig_y,trans_x,trans_y\n");
    for (int i = 0; i < VERTICES; i++) {
        fprintf(fp, "%.2f,%.2f,%.2f,%.2f\n", 
                original[0][i], original[1][i], 
                transformed[0][i], transformed[1][i]);
    }
    
    fclose(fp);
    printf("\nPolygon data saved to %s\n", filename);
}

int main() {
    // Define dimensions for letter 'E'
    double e_height = 10.0;
    double e_thickness = 2.0;
    double e_width = 6.0;
    
    // Define letter 'E' polygon (14 vertices in homogeneous coordinates)
    double letter[ROWS][VERTICES] = {
        {0, e_thickness, e_thickness, e_width, e_width, e_thickness, e_thickness, e_width, e_width, e_thickness, e_thickness, e_width, e_width, 0},
        {0, 0, e_thickness, e_thickness, e_thickness + e_thickness, e_thickness + e_thickness, e_height - e_thickness - e_thickness, e_height - e_thickness - e_thickness, e_height - e_thickness, e_height - e_thickness, e_height, e_height, 0, 0},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
    };
    
    // Define rotation matrix (45 degrees)
    double theta = M_PI / 4.0;
    double rotation[ROWS][ROWS] = {
        {cos(theta), -sin(theta), 0},
        {sin(theta), cos(theta), 0},
        {0, 0, 1}
    };
    
    // Define translation matrix (tx=20, ty=-5)
    double translation[ROWS][ROWS] = {
        {1, 0, 20},
        {0, 1, -5},
        {0, 0, 1}
    };
    
    // Define shear matrix (shear factor 2 in x-direction)
    double shear[ROWS][ROWS] = {
        {1, 2, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    
    // Define scale matrix (sx=0.5, sy=2)
    double scale[ROWS][ROWS] = {
        {0.5, 0, 0},
        {0, 2, 0},
        {0, 0, 1}
    };
    
    // Define reflection matrix (reflect about y-axis)
    double reflection[ROWS][ROWS] = {
        {-1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    
    // Combine transformations: reflection * scale * shear * translation * rotation
    double temp1[ROWS][ROWS], temp2[ROWS][ROWS], temp3[ROWS][ROWS], temp4[ROWS][ROWS];
    double transformation[ROWS][ROWS];
    
    matrix_mult_3x3(translation, rotation, temp1);      // translation * rotation
    matrix_mult_3x3(shear, temp1, temp2);               // shear * (translation * rotation)
    matrix_mult_3x3(scale, temp2, temp3);               // scale * (shear * translation * rotation)
    matrix_mult_3x3(reflection, temp3, transformation); // reflection * (scale * shear * translation * rotation)
    
    // Apply transformation to letter 'E'
    double new_poly[ROWS][VERTICES];
    matrix_mult_polygon(transformation, letter, new_poly);
    
    // Print results
    printf("=== Letter 'E' Polygon Transformation ===\n");
    printf("\nDimensions:\n");
    printf("  Height: %.2f\n", e_height);
    printf("  Thickness: %.2f\n", e_thickness);
    printf("  Width: %.2f\n", e_width);
    
    print_polygon("Original Letter 'E'", letter);
    print_polygon("Transformed Letter 'E'", new_poly);
    
    // Print transformation details
    printf("\nTransformations applied (in order):\n");
    printf("  1. Rotation: %.0f degrees\n", theta * 180.0 / M_PI);
    printf("  2. Translation: (20, -5)\n");
    printf("  3. Shear: factor 2 in x-direction\n");
    printf("  4. Scale: (0.5, 2)\n");
    printf("  5. Reflection: about y-axis\n");
    
    // Save to CSV for external plotting (optional)
    save_to_csv("polygon_data.csv", letter, new_poly);
    
    printf("\nTo visualize the polygons, you can use the CSV file with plotting software.\n");
    
    return 0;
}
