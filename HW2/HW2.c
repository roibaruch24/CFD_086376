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
    //const char *parameters_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\parametrs.txt";
    const char *parameters_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\parametrs.txt";
    //const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\Grid_File2.txt";
    const char *input_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\Grid_File2.txt";

    // Open parameters file
    FILE *input_param = fopen(parameters_file_path, "rt");
    if (input_param == NULL) 
    {
        perror("Error opening parameters file");
        return -1; // Return an error code
    }

    // Read Mach number, alpha, pressure, and density from parameters file
    if (fscanf(input_param, "%lf %lf %lf %lf", M, alpha, p, density) != 4) 
    {
        printf("Error reading M, alpha, p, density from the parameters file.\n");
        fclose(input_param);
        return -1;
    }
    *alpha = *alpha/57.3;

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
    //const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\Grid_File2.txt";
    const char *input_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\Grid_File2.txt";

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
    for (int i=0; i<=imax; i++) //first order
    {
        for (int j=0; j<=jmax; j+=jmax)
        {
            if (j==0) //forward
            {
                x_etta[offset2d(i,j,ni)] = (x[offset2d(i,j+1,ni)]-x[offset2d(i,j,ni)]);
                y_etta[offset2d(i,j,ni)] = (y[offset2d(i,j+1,ni)]-y[offset2d(i,j,ni)]);
            }
            else //backward
            {
                x_etta[offset2d(i,j,ni)] = (x[offset2d(i,j,ni)]-x[offset2d(i,j-1,ni)]);
                y_etta[offset2d(i,j,ni)] = (y[offset2d(i,j,ni)]-y[offset2d(i,j-1,ni)]);
            }
        }
    }
    for (int j=0; j<=jmax; j++) //first order
    {
        for (int i=0; i<=imax; i+=imax) //xi derivatives arent defined on i min and i max so wont be calculated
        {
            if (i==0)
            {
                x_ksi[offset2d(i,j,ni)] = (x[offset2d(i+1,j,ni)]-x[offset2d(i,j,ni)]);
                y_ksi[offset2d(i,j,ni)] = (y[offset2d(i+1,j,ni)]-y[offset2d(i,j,ni)]);
            }
            else
            {
                x_ksi[offset2d(i,j,ni)] = (x[offset2d(i,j,ni)]-x[offset2d(i-1,j,ni)]);
                y_ksi[offset2d(i,j,ni)] = (y[offset2d(i,j,ni)]-y[offset2d(i-1,j,ni)]);
            }
        }
    }
}
void freestream (double density, double p, double M, double alpha, double *Q, int imax, int jmax, int ni, int nj){
    double a = sqrt(gamma*density/p);
    double u = M*a*cos(alpha/57.3);
    double v = M*a*sin(alpha/57.3);
    for(int i = 0; i<= imax; i++){
        for( int j=0; j<= jmax; j++){
            Q[offset3d(i, j, 0, ni, nj)] = density;
            Q[offset3d(i, j, 1, ni, nj)] = density*u;
            Q[offset3d(i, j, 2, ni, nj)] = density*v;
            Q[offset3d(i, j, 3, ni, nj)] = p/(gamma-1) + 0.5*density*(pow(u,2)+pow(v,2));
        }
    }

}
void bc_wall(double x, double y, double u, double v, int tel, int teu, double *Q, int ni, int nj, double *y_ksi, double *x_ksi, double *x_etta, double *y_etta){
    for (int i = tel; i<=teu; i++){
            double u_1 = Q[offset3d(i, 1, 1,ni, nj)]/Q[offset3d(i, 0, 1, ni, nj)];
            double v_1 = Q[offset3d(i, 1, 2,ni, nj)]/Q[offset3d(i, 1, 0, ni, nj)];
            double u = (x_ksi[offset2d(i,0,ni)]*((y_etta[offset2d(i,1,ni)])*u_1-(x_etta[offset2d(i,1,ni)])*v_1))/
            ((x_ksi[offset2d(i,0,ni)]*(y_etta[offset2d(i,0,ni)]))-(x_etta[offset2d(i,0,ni)])*(y_ksi[offset2d(i,0,ni)]));
            double v =(y_ksi[offset2d(i,0,ni)])/(x_ksi[offset2d(i,0,ni)])*u;
            double density = Q[offset3d(i,1,0,ni,nj)];
            double e1 = Q[offset3d(i,1,3,ni,nj)];
            double P = (gamma-1)*(e1-0.5*density*(pow(u_1,2)+(pow(v_1,2))));
            double e = P/(gamma-1)+0.5*density*(pow(u,2)+(pow(v,2)));
            Q[offset3d(i,0,0,ni,nj)] = density;
            Q[offset3d(i,0,1,ni,nj)] = density*u;
            Q[offset3d(i,0,2,ni,nj)] = density*v;
            Q[offset3d(i,0,3,ni,nj)] = e;
    }

}
void RHS(int ni, int nj, int imax, int jamx, double *Jacobian, double *y_ksi, double *x_ksi, double *x_etta, double *y_etta, double *Q, double *W, double *S){   
    /* ksi direction */ 
    for (int j = 1; j < jamx; j++){
        for ( int i = 0; i <= imax; i++){
            double u = Q[offset3d(i, j, 1, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double v = Q[offset3d(i, j, 2, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double U_conv = Jacobian[offset2d(i, j, ni)]*(y_etta[offset2d(i, j, ni)]*u - x_etta[offset2d(i, j, ni)]*v);
            double P = ((gamma-1)*Q[offset3d(i, j, 3, ni, nj)]-0.5*Q[offset3d(i, j, 0, ni, nj)]*(pow(u,2)+pow(v,2)))/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 0, ni)] = (Q[offset3d(i, j, 0, ni, nj)]*U_conv)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 1, ni)] = (Q[offset3d(i, j, 1, ni, nj)]*U_conv + Jacobian[offset2d(i, j, ni)]*y_etta[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 2, ni)] = (Q[offset3d(i, j, 2, ni, nj)]*U_conv - Jacobian[offset2d(i, j, ni)]*x_etta[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 3, ni)] = U_conv*(P+Q[offset3d(i, j, 3, ni, nj)])/Jacobian[offset2d(i, j, ni)];
    }
    for(int i = 1; i< imax; i++){
        for( int k = 0; k < 4; k++ ){
            S[offset3d(i, j, k, ni, nj)] = -0.5*(W[offset2d(i+1,k, ni)]-W[offset2d(i-1,k,ni)]);
        }
    }

}

/* etta direction */ 
    for (int i = 1; i < imax; i++){
        for ( int j = 0; j <= jamx; j++){
            double u = Q[offset3d(i, j, 1, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double v = Q[offset3d(i, j, 2, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double V_conv = Jacobian[offset2d(i, j, ni)]*(-y_ksi[offset2d(i, j, ni)]*u + x_ksi[offset2d(i, j, ni)]*v);
            double P = ((gamma-1)*Q[offset3d(i, j, 3, ni, nj)]-0.5*Q[offset3d(i, j, 0, ni, nj)]*(pow(u,2)+pow(v,2)))/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 0, ni)] = (Q[offset3d(i, j, 0, ni, nj)]*V_conv)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 1, ni)] = (Q[offset3d(i, j, 1, ni, nj)]*V_conv - Jacobian[offset2d(i, j, ni)]*y_ksi[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 2, ni)] = (Q[offset3d(i, j, 2, ni, nj)]*V_conv + Jacobian[offset2d(i, j, ni)]*x_ksi[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 3, ni)] = V_conv*(P+Q[offset3d(i, j, 3, ni, nj)])/Jacobian[offset2d(i, j, ni)];
        }
        for(int j = 1; j < jamx ; j++){
            for( int k = 0; k < 4; k++ ){
                S[offset3d(i, j, k, ni, nj)] = -0.5*(W[offset2d(j+1,k, ni)]-W[offset2d(j-1,k,ni)]); 
            }
        }

    }

}
int print_output(int ni,int nj,
double *x, double *y)
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
    double M_0, alpha, p_0, density_0; // Ensure these are doubles
    // Read the parameters (ni, nj, teu, tel, Mach number, etc.)
    if (input(&ni, &nj, &teu, &tel, &M_0, &alpha, &p_0, &density_0) != 0) {
        printf("Failed to read input parameters.\n");
        return -1; // Exit if input fails
    }
    printf("ni: %d, nj: %d, teu: %d, tel: %d\n", ni, nj, teu, tel);
    printf("Mach number: %f, Alpha: %f, Pressure: %f, Density: %f\n", M_0, alpha, p_0, density_0);

    // Allocate memory for x and y in main
    int num_points = ni * nj;
    int jmax = nj - 1;
    int imax = ni - 1;

    double *x = (double *)malloc(num_points * sizeof(double));
    double *y = (double *)malloc(num_points * sizeof(double));

    // Check memory allocation
    if (x == NULL || y == NULL) {
        printf("Memory allocation failed!\n");
        return -1;
    }

    // Read the mesh data
    if (input_mesh(x, y, ni, nj) != 0) {
        printf("Failed to read mesh data.\n");
        free(x);
        free(y);
        return -1; // Exit if input fails
    }

    // Example of printing mesh data for verification
    printf("Mesh coordinates (x, y):\n");
    for (int i = 0; i < num_points; i++) {
        printf("Point %d: x = %f, y = %f\n", i, x[i], y[i]);
    }

    // Perform calculations
    double *y_ksi = (double *)malloc(num_points * sizeof(double));
    double *x_ksi = (double *)malloc(num_points * sizeof(double));
    double *x_etta = (double *)malloc(num_points * sizeof(double));
    double *y_etta = (double *)malloc(num_points * sizeof(double));
    double *jacobian = (double *)malloc(num_points * sizeof(double));
    //double *u = (double *)malloc(num_points * sizeof(double));
    //double *v = (double *)malloc(num_points * sizeof(double));
    double *p = (double *)malloc(num_points * sizeof(double));
    double *density = (double *)malloc(num_points * sizeof(double));
    double *Q = (double *)malloc(num_points * 4 * sizeof(double));
    double *S = (double *)malloc(num_points * 4 * sizeof(double));
    double *W = (double *)malloc(MAX(ni, nj) * 4 * sizeof(double));
    double *S = (double *)malloc(MAX(ni,nj) * 4 * sizeof(double));
    
    

    // Verify memory allocation for all arrays
    if (!y_ksi || !x_ksi || !x_etta || !y_etta || !jacobian || !p || !density || !Q) {
        printf("Memory allocation failed!\n");
        free(x);
        free(y);
        return -1;
    }

    matrics(imax, jmax, ni, x, y, y_ksi, x_ksi, x_etta, y_etta, jacobian);
    freestream (density_0, p_0, M_0, alpha, Q, imax, jmax, ni, nj);
    print_output(ni, nj, Q, jacobian);

    // Free the dynamically allocated memory
    free(x);
    free(y);
    free(y_ksi);
    free(x_ksi);
    free(x_etta);
    free(y_etta);
    free(jacobian);
    free(p);
    free(density);
    free(Q);

    printf("Memory freed. Program finished successfully.\n");
    return 0; 
}

