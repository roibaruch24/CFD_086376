#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.1415926
void calc_h(float *h, float N)
{
    *h=1/N; // define step size as 1/N, when N is the number of points in the interval
}
void calc_a(float *a,float h,float N)
{
    for (int i = 0;i<=N; i++)
    {
        a[i]=1;
    }
}
void calc_b(float *b,float h,float N)
{
    for (int i = 0;i<=N; i++)
    {
        b[i]=i*h;
    }
}
void calc_c(float *c,float h,float N)
{
    for (int i = 0;i<=N; i++)
    {
        c[i]=(i*h)*(i*h);
    }
}
void calc_A(float *A, float *a,float *b,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        A[i]=a[i]/(h*h) - b[i]/(2*h);
    }
}
void calc_B(float *B,float *a,float *c,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        B[i]=-2*a[i]/(h*h) +c[i]*c[i];
    }
}
void calc_C(float *C,float *a,float *b,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        C[i]=a[i]/(h*h) +b[i]/(2*h);
    }
}
void RHS(float *d,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        d[i]=sin(2*PI*i*h)+cos(2*PI*i*h);
    }
}
typedef enum {
    Dirichlet,
    Neumann
} OperationType;
float LHS(float *A_vec, float *B_vec, float *C_vec, float *A, float *B, float *C,float N, OperationType boundary_conditions){
    switch (boundary_conditions){
        case Dirichlet:
        for (int i = 0; i < N-2; i++) {
                A_vec[i] = A[i+2];
                B_vec[i] = B[i+1];
                C_vec[i] = C[i+1];
            }
            break;
            B_vec[(int)N-1] = B[(int)N-1];  
        case Neumann:
        for (int i = 0; i < N; i++) {
                A_vec[i] = A[i+1];
                B_vec[i] = B[i];
                C_vec[i] = C[i];
            }
            A_vec[(int)N-1] = A[(int)N] + C[(int)N];
             C_vec[0] = A[0] + C[0];
            break;
    }  
}
int tridiag(float *A_vec, float *B_vec, float *C_vec, float *d, float *u, int is, int ie)
{

  int i;
  float beta;

  for (i = is + 1; i <= ie; i++)
    {
      if(B_vec[i-1] == 0.) return(1);
      beta = A_vec[i] / B_vec[i-1];
      B_vec[i] = B_vec[i] - C_vec[i-1] * beta;
      d[i] = d[i] - d[i-1] * beta;
    }

  u[ie] = d[ie] / B_vec[ie];
  for (i = ie - 1; i >= is; i--)
    {
      u[i] = (d[i] - C_vec[i] * u[i+1]) / B_vec[i];
    }
  return(0);
}
int main() {
    float N = 5; // Example value for N
    float h;
    float *u = (float *)malloc((N + 1) * sizeof(float));
    float *a = (float *)malloc((N + 1) * sizeof(float));
    float *b = (float *)malloc((N + 1) * sizeof(float));
    float *c = (float *)malloc((N + 1) * sizeof(float));
    float *A = (float *)malloc((N + 1) * sizeof(float));
    float *B = (float *)malloc((N + 1) * sizeof(float));
    float *C = (float *)malloc((N + 1) * sizeof(float));
    float *d = (float *)malloc((N + 1) * sizeof(float));
    float *A_vec = (float *)malloc((N) * sizeof(float));
    float *B_vec = (float *)malloc((N) * sizeof(float));
    float *C_vec = (float *)malloc((N) * sizeof(float));
    calc_h(&h, N);
    calc_a(a, h, N);
    calc_b(b, h, N);
    calc_c(c, h, N);
    calc_A(A, a, b, h, N);
    calc_B(B, a, c, h, N);
    calc_C(C, a, b, h, N);
    RHS(d, h, N);
    LHS(A_vec, B_vec, C_vec, A, B, C, N, Dirichlet);
    tridiag(A_vec, B_vec, C_vec, d, u, 0, N - 1);
    FILE *f = fopen("C:\\Users\\roiba\\Documents\\CFD_086376\\HW0\\output.dat", "w");
    if (f == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        return 1;
    }
    for (int i = 0; i <= N; i++) {
        fprintf(f, "%f %f\n", i * h, u[i]);
    }
    fclose(f); 
    free(a);
    free(b);
    free(c);
    free(d);
    free(A_vec);
    free(B_vec);
    free(C_vec);
    free(u);
    return 0;
}
    