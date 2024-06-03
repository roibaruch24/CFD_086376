#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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
/*void calc_A(float *A, float *a,float *b,float h,float N)
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
}*/
void LHS (float *a, float *b, float *c, float *A_dig, float *B_dig, float *C_dig, float h, float N)
{
    for (int i = 0; i <= N; i++)
    {
        A_dig[i]=a[i]/(h*h) - b[i]/(2*h);
        B_dig[i]=-2*a[i]/(h*h) +c[i]*c[i];
        C_dig[i]=a[i]/(h*h) +b[i]/(2*h);
    }
}
void RHS(float *d,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        d[i]=sin(2*M_PI*i*h)+cos(2*M_PI*i*h);
    }
}
typedef enum {
    Dirichlet,
    Neumann
} OperationType;
void Boundary_conditions(float *A_dig, float *B_dig, float *C_dig, float *d, float N, float h, int *is, int *ie, OperationType Boundary_conditions)
{
    switch (Boundary_conditions){
        case Dirichlet:
        *is = 1;
        *ie = N-1;
        d[*is] -= A_dig[*is]*0;
        d[*ie] -= C_dig[*ie]*1 ;
        break;
        case Neumann:
        *is = 0;
        *ie = N;
        C_dig[*is] += A_dig[*is];
        A_dig[*ie] += C_dig[*ie];
        d[*is] += 2*h*A_dig[*is]*1;
        d[*ie] -= 2*h*C_dig[*ie]*(-1) ;
        break;
    }
}
/*float LHS(float *A_vec, float *B_vec, float *C_vec, float *A, float *B, float *C,float N, OperationType boundary_conditions){
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
}*/
int tridiag(float *A_dig, float *B_dig, float *C_dig, float *d, float *u, int is, int ie)
{
  int i;
  float beta;

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
int main() {
    float N = 20; // Example value for N
    float h;
    int ie, is;
    float *u = (float *)malloc((N + 1) * sizeof(float));
    float *a = (float *)malloc((N + 1) * sizeof(float));
    float *b = (float *)malloc((N + 1) * sizeof(float));
    float *c = (float *)malloc((N + 1) * sizeof(float));
    float *d = (float *)malloc((N + 1) * sizeof(float));
    float *A_dig = (float *)malloc((N+1) * sizeof(float));
    float *B_dig = (float *)malloc((N+1) * sizeof(float));
    float *C_dig = (float *)malloc((N+1) * sizeof(float));
    calc_h(&h, N);
    calc_a(a, h, N);
    calc_b(b, h, N);
    calc_c(c, h, N);
    RHS(d, h, N);
    LHS(a, b, c, A_dig, B_dig, C_dig, h, N);
    Boundary_conditions(A_dig, B_dig, C_dig, d, N, h, &is, &ie, Dirichlet);
    tridiag(A_dig, B_dig, C_dig, d, u, is, ie);
    for (int i = 0; i <= N; i++) {
        printf("u[%d] = %f\n", i, u[i]);
    }
     const char *file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW0\\solution.dat";
    FILE *file = fopen(file_path, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }
    
    for (int i = 0; i <= N; i++) {
        fprintf(file, "%d %f\n", i, u[i]);
    }
    fclose(file);
    free(a);
    free(b);
    free(c);
    free(d);
    free(A_dig);
    free(B_dig);
    free(C_dig);
    free(u);
    return 0;
}
    