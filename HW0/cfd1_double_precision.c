/* 7/6/24, Roi Baruch
This code id HW0 in Computational Fluid Dynamics course No. 086376
link to Git repository: https://github.com/roibaruch24/CFD_086376
*/
// Including all the relevent headers:
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
// Init - reads the input file and initializes all the relevent valuse
int Init(int *N, char *boundary_condition, double *start,double *end,double *bc_0,double *bc_n){
    const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW0\\input.txt";
    FILE *input = fopen(input_file_path, "rt");
    if (input == NULL) {
        return 1;
    }
    fscanf(input,"%d %lf %lf %c %lf %lf",N,start,end,boundary_condition,bc_0,bc_n);
    fclose(input);
    return 0;
}
// calc_h - calculates the step value h
void calc_h(double *h, double start, double end, int N)
{
    *h=(end - start) / N;
}
// calc_a - initialize array 'a' to 1
void calc_a(double *a, int N) {
    for (int i = 0; i <= N; i++) {
        a[i] = 1;
    }
}
// calc_b - initialize array 'b' to x
void calc_b(double *b,double h,int N)
{
    for (int i = 0;i<=N; i++)
    {
        b[i]=i*h;
    }
}
// calc_c - initialize array 'c' to x^2
void calc_c(double *c,double h,int N)
{
    for (int i = 0;i<=N; i++)
    {
        c[i]=(i*h)*(i*h);
    }
}
// calc_A - initialize the lower diagonal array 'A'
void calc_A(double *A_dig, double *a,double *b,double h,int N)
{
    for (int i = 0; i<=N; i++)
    {
        A_dig[i]=a[i]/(h*h) - b[i]/(2*h);
    }
}
// calc_B - initialize the middle diagonal array 'B'
void calc_B(double *B_dig,double *a,double *c,double h,int N)
{
    for (int i = 0; i<=N; i++)
    {
        B_dig[i]=-2*a[i]/(h*h) +c[i];
    }
}
// calc_C - initialize the upper diagonal array 'C'
void calc_C(double *C_dig,double *a,double *b,double h,int N)
{
    for (int i = 0; i<=N; i++)
    {
        C_dig[i]=a[i]/(h*h) +b[i]/(2*h);
    }
}
// LHS - calls the the diagonals for the LHS of the function
void LHS (double *a, double *b, double *c, double *A_dig, double *B_dig, double *C_dig, double h, int N)
{
    calc_A(A_dig, a, b, h, N);
    calc_B(B_dig, a, c, h, N);
    calc_C(C_dig, a, b, h, N);
}
// RHS - initialize arry 'd' for the RHS of the function
void RHS(double *d,double h,int N)
{
    for (int i = 0; i<=N; i++)
    {
        double x = 2*3.1415*i*h;
        d[i] = sin(x) + cos(x);
    }
}
// Boundary_conditions - applies the boundary conditions
void Boundary_conditions(double *A_dig, double *C_dig, double *d, int N, double h, int *is, int *ie, double *u, char boundary_condition, double bc_0, double bc_N) {
    if (boundary_condition == 'D')
    {
        *is = 1;
        *ie = N-1;
        d[*is] -= A_dig[*is]*bc_0;
        d[*ie] -= C_dig[*ie]*bc_N;
        u[0] = bc_0;
        u[N] = bc_N;

    }
    else if (boundary_condition == 'N')
    {
    *is = 0;
    *ie = N;
    C_dig[*is] += A_dig[*is];
    A_dig[*ie] += C_dig[*ie];
    d[*is] += 2*h*A_dig[*is]*bc_0;
    d[*ie] -= 2*h*C_dig[*ie]*bc_N ;
    }
}
// initialize A, B, C vectors for error calculation
void error_calc(double *A_dig, double *B_dig, double *C_dig, double *d, double *A_dig_er, double *B_dig_er, double *C_dig_er, double *d_er, int N ){
    for (int i = 0; i<=N; i++)
    {
    A_dig_er[i] = A_dig[i];
    B_dig_er[i] = B_dig[i];
    C_dig_er[i] = C_dig[i];
    d_er[i] = d[i]; 
    }
}
// Tridiag - is the solver
int tridiag(double *A_dig, double *B_dig, double *C_dig, double *d, double *u, int is, int ie)
{
  int i;
  double beta;
  for (i = is + 1; i <= ie; i++)
    {
      if(B_dig[i-1] == 0.) return(1);
      beta = A_dig[i] / B_dig[i-1];
      B_dig[i] = B_dig[i] - C_dig[i-1] * beta;
      d[i] = d[i] - d[i-1] * beta;
    }

  u[ie] = d[ie] / B_dig[ie];
  for (i = ie - 1; i >= is; i--)
    {
      u[i] = (d[i] - C_dig[i] * u[i+1]) / B_dig[i];
    }
  return(0);
}
// output - writes the output files
/*int output (int N, double h, double *u, double *A_dig, double *B_dig, double *C_dig)
{
    const char *solution_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW0\\solution.dat";
    const char *error_calc_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW0\\error_calc.dat";
    FILE *solution = fopen(solution_file_path, "wt");
    FILE *error_calc = fopen(error_calc_file_path, "wt");
    if (solution == NULL) {
        return 1;
    }
    for (int i = 0; i <= N; i++) {
       fprintf(solution, "%f %f\n", i * h, u[i]);
       fprintf(error_calc, "%f %f %f %f %f\n", i * h, A_dig[i], B_dig[i], C_dig[i], d[i]);
    }
    fclose(solution);
    return 0;
}*/
//int main() {

    void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // declaring all the parametrs:
    double h, bc_0, bc_N, start, end;
    int ie, is, N;
    char boundary_condition;
    
    // call Init function
    Init(&N, &boundary_condition, &start, &end, &bc_0, &bc_N);
    
    // Allocate memory for arrays   
    double *u = malloc((N + 1) * sizeof(double));
    double *a = malloc((N + 1) * sizeof(double));
    double *b = malloc((N + 1) * sizeof(double));
    double *c = malloc((N + 1) * sizeof(double));
    double *d = malloc((N + 1) * sizeof(double));
    double *A_dig = malloc((N + 1) * sizeof(double));
    double *B_dig = malloc((N + 1) * sizeof(double));
    double *C_dig = malloc((N + 1) * sizeof(double));
    double *A_dig_er = malloc((N + 1) * sizeof(double));
    double *B_dig_er = malloc((N + 1) * sizeof(double));
    double *C_dig_er = malloc((N + 1) * sizeof(double));
    double *d_er = malloc((N + 1) * sizeof(double));
    
    // calls all the main functions
    calc_h(&h, start, end, N);
    calc_a(a, N);
    calc_b(b, h, N);
    calc_c(c, h, N);
    RHS(d, h, N);
    LHS(a, b, c, A_dig, B_dig, C_dig, h, N);
    Boundary_conditions(A_dig, C_dig, d, N, h, &is, &ie, u, boundary_condition, bc_0, bc_N);
    error_calc(A_dig, B_dig, C_dig, d, A_dig_er, B_dig_er, C_dig_er, d_er, N);
    tridiag(A_dig, B_dig, C_dig, d, u, is, ie);
    //output(N, h, u, A_dig, B_dig, C_dig, d); IN A COMMENT BECAUSE IM USING MATLAB
    
   // Prepare output mxArrays
    plhs[0] = mxCreateNumericMatrix(1, N + 1, mxDOUBLE_CLASS, mxREAL); 
    plhs[1] = mxCreateNumericMatrix(1, N + 1, mxDOUBLE_CLASS, mxREAL); 
    plhs[2] = mxCreateNumericMatrix(1, N + 1, mxDOUBLE_CLASS, mxREAL); 
    plhs[3] = mxCreateNumericMatrix(1, N + 1, mxDOUBLE_CLASS, mxREAL); 
    plhs[4] = mxCreateNumericMatrix(1, N + 1, mxDOUBLE_CLASS, mxREAL); 
    
    double *output_u = (double *)mxGetData(plhs[0]);
    double *output_A_dig = (double *)mxGetData(plhs[1]);
    double *output_B_dig = (double *)mxGetData(plhs[2]);
    double *output_C_dig = (double *)mxGetData(plhs[3]);
    double *output_d_er = (double *)mxGetData(plhs[4]);
    
    for (int i = 0; i <= N; i++) {
        output_u[i] = (double)u[i];
        output_A_dig[i] = (double)A_dig_er[i];
        output_B_dig[i] = (double)B_dig_er[i];
        output_C_dig[i] = (double)C_dig_er[i];
        output_d_er[i]  = (double)d_er[i];
    }
    // Free allocated memory
    free(a);
    free(b);
    free(c);
    free(d);
    free(A_dig);
    free(B_dig);
    free(C_dig);
    free(A_dig_er);
    free(B_dig_er);
    free(C_dig_er);
    free(d_er);
    free(u);
}

    