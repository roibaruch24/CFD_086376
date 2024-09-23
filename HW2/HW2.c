#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI  3.14159265359
#define gamma 1.4
#define MAX(a, b) ((a) > (b) ? (a) : (b))
// Makes a 2D array to a vector
int offset2d(int i, int j, double ni) {
    return j * ni + i;
}
int offset3d(int const i, int const j,int const k, int const ni ,int const nj)
{
    return (k * nj +j) * ni + i;
}
// Reads the input file for parameters
int input(int *ni, int *nj, int *teu, int *tel, double *M, double *alpha, double *p, double *density) 
{
    const char *parameters_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\parametrs.txt";
    const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\Grid_File2.txt";

    // Open parameters file
    FILE *input_param = fopen(parameters_file_path, "rt");
    if (input_param == NULL) 
    {
        perror("Error opening parameters file");
        return -1; // Return an error code
    }

    // Read Mach number, alpha, pressure, and density from parameters file
    if (fscanf(input_param, "%f %f %f %f", M, alpha, p, density) != 4) 
    {
        printf("Error reading M, alpha, p, density from the parameters file.\n");
        fclose(input_param);
        return -1;
    }

    fclose(input_param); // Close parameters file

    // Open mesh file
    FILE *input = fopen(input_file_path, "rt");
    if (input == NULL) 
    {
        perror("Error opening mesh file");
        return -1; // Return an error code
    }

    // Read ni, nj, teu, tel from mesh file
    if (fscanf(input, "%d %d %d %d", ni, nj, teu, tel) != 4) 
    {
        printf("Error reading ni, nj, teu, tel from the mesh file.\n");
        fclose(input);
        return -1;
    }

    fclose(input); // Close mesh file
    return 0; // Success
}

// Reads the mesh file (but assumes memory for x and y is already allocated)
int input_mesh(double *x, double *y, int ni, int nj) 
{
    const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\Grid_File2.txt";

    // Open mesh input file
    FILE *input = fopen(input_file_path, "rt");
    if (input == NULL) 
    {
        perror("Error opening mesh file");
        return -1; // Return an error code
    }

    // Skip the first line that contains ni, nj, teu, tel
    fscanf(input, "%*d %*d %*d %*d");

    // Read the x and y values
    int num_points = ni * nj;
    for (int i = 0; i < num_points; i++) 
    {
        if (fscanf(input, "%lf %lf", &x[i], &y[i]) != 2) 
        {
            printf("Error reading x and y values at index %d.\n", i);
            fclose(input);
            return -1;
        }
    }

    fclose(input); // Close mesh file
    return 0; // Success
}
void matrics(int imax, int jmax, int ni, double *x, double *y, double *y_ksi, double *x_ksi, double *x_etta, double *y_etta, double *jacobian)
{
    for (int j = 0; j <= jmax; j++){
        for (int i = 1; i < imax; i++){
            x_ksi[offset2d(i,j,ni)] = ((x[offset2d(i+1,j,ni)]-x[offset2d(i-1,j,ni)])/2);
            y_ksi[offset2d(i,j,ni)] = ((y[offset2d(i+1,j,ni)]-y[offset2d(i-1,j,ni)])/2);    
        }
    }
    for (int i = 0; i <= imax; i++){
        for (int j = 1; j < jmax; j++){
            x_etta[offset2d(i,j,ni)] = ((x[offset2d(i,j+1,ni)]-x[offset2d(i,j-1,ni)])/2);
            y_etta[offset2d(i,j,ni)] = ((y[offset2d(i,j+1,ni)]-y[offset2d(i,j-1,ni)])/2);
        }
    }
    for (int i =1; i< imax; i++){
        for (int j = 1; j<jmax; j++)
        {
            jacobian[offset2d(i,j,ni)] = 1/(x_ksi[offset2d(i,j,ni)]*y_etta[offset2d(i,j,ni)] - y_ksi[offset2d(i,j,ni)]*x_etta[offset2d(i,j,ni)]);
        }
    }
}
void freestream (double density, double p, double M, double alpha, double *u, double *v){
    double a = sqrt(gamma*density/p);
    *u = M*a*cos(alpha);
    *v = M*a*sin(alpha);
}
void bc_wall(double x, double y, double u, double v, int tel, int teu, double *density, double *p, double *e, int ni,double y_ksi, double x_ksi, double x_etta, double y_etta){
    for (int i = tel; i<=teu; i++){
            
            density[offset2d(i,0,ni)] = density[offset2d(i,1,ni)];
            p[offset2d(i,0,ni)] = p[offset2d(i,1,ni)];
            e[offset2d(i,0,ni)] = p[offset2d(i,0,ni)]/(gamma-1) + 0.5*density[offset2d(i,0,ni)]*(u+v);

    }

}
int print_output(int ni,int nj,double *x, double *y)
{
    char *output_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\output.txt";
    FILE *output= fopen(output_file_path, "wt"); //crating output file
    for (int i = 0; i < ni*nj; i++){
            fprintf(output, "%f %f\n", x[i], y[i]);
        }
    
    fclose(output);
    return 0;
}

int main() {
    int ni, nj, teu, tel;
    double M_0, alpha, p_0, density_0;
    // Read the parameters (ni, nj, teu, tel, Mach number, etc.)
    input(&ni, &nj, &teu, &tel, &M_0, &alpha, &p_0, &density_0);
    // Allocate memory for x and y in main
    int num_points = ni * nj;
    int jmax = nj-1;
    int imax = ni-1;
    double *x = (double *)malloc(num_points * sizeof(double));
    double *y = (double *)malloc(num_points * sizeof(double));
    double *x_ksi = (double *)malloc(num_points * sizeof(double));
    double *y_ksi = (double *)malloc(num_points * sizeof(double));
    double *x_etta = (double *)malloc(num_points * sizeof(double));
    double *y_etta = (double *)malloc(num_points * sizeof(double));
    double *jacobian = (double *)malloc(num_points * sizeof(double));
    double *u = (double *)malloc(num_points * sizeof(double));
    double *v = (double *)malloc(num_points * sizeof(double));
    double *p = (double *)malloc(num_points * sizeof(double));
    double *density = (double *)malloc(num_points * sizeof(double));

    // Read the mesh data
    input_mesh(x, y, ni, nj);
    matrics(imax, jmax,ni, x, y, y_ksi, x_ksi, x_etta, y_etta, jacobian);
    print_output(ni, nj, jacobian, y_etta);
    // Free the dynamically allocated memory
    free(x);
    free(y);
    

    return 0; 
}
