#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.1415926


void calc_h(float *h, float N)
{
    *h=1/N; // define step size as 1/N, when N is the number of points in the interval
}
void calc_A(float *a,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        a[i]=1/(h*h) - i/(2*h);
    }
}
void calc_B(float *b,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        b[i]=-2/(h*h) +i*i;
    }
}
void calc_C(float *c,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        c[i]=1/(h*h) +i/(2*h);
    }
}
void calc_D(float *d,float h,float N)
{
    for (int i = 0; i<=N; i++)
    {
        d[i]=sin(2*PI*i)+cos(2*PI*i);
    }
}
typedef enum {
    Dirichlet,
    Neumann
} OperationType;

double calc_dig_A(float *A_vec,float *a, float N,float *c, OperationType boundary_conditions) {
    switch (boundary_conditions) {
        case Dirichlet:
                for (int i = 2; i < N; i++) {
                A_vec[i-2] = a[i];
            }
            break;       
        case Neumann:
            for (int i = 0; i < N-1; i++) {
                A_vec[i] = a[i+1];
            }
            A_vec[(int)N-1] = a[(int)N] + c[(int)N];
            break;
    }
}
double calc_dig_B(float *B_vec, float *b, float N,float *c, OperationType boundary_conditions) {
    switch (boundary_conditions) {
        case Dirichlet:
                for (int i = 0; i < N-1; i++) {
                B_vec[i] = b[i+1];
            }
            break;       
        case Neumann:
             for (int i = 0; i <= N; i++) {
                B_vec[i] = b[i];
            }
            break;
           
    }
    }
    double calc_dig_C(float *C_vec, float *c, float N,float *a, OperationType boundary_conditions) {
    switch (boundary_conditions) {
        case Dirichlet:
                for (int i = 1; i < N-1; i++) {
                C_vec[i-1] = c[i];
            }
            break;       
        case Neumann:
             for (int i = 1; i < N; i++) {
                C_vec[i] = a[i];
            }
            C_vec[0] = a[0] + c[0];
            break;
           
    }
    }
    int tridiag(float *a, float *b, float *c, float *d, float *u, int is, int ie)
{

  int i;
  float beta;

  for (i = is + 1; i <= ie; i++)
    {
      if(b[i-1] == 0.) return(1);
      beta = a[i] / b[i-1];
      b[i] = b[i] - c[i-1] * beta;
      d[i] = d[i] - d[i-1] * beta;
    }

  u[ie] = d[ie] / b[ie];
  for (i = ie - 1; i >= is; i--)
    {
      u[i] = (d[i] - c[i] * u[i+1]) / b[i];
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
    float *d = (float *)malloc((N + 1) * sizeof(float));
    float *A_vec = (float *)malloc((N - 2) * sizeof(float)); // Allocate memory for A_vec
    float *B_vec = (float *)malloc((N - 1) * sizeof(float)); // Allocate memory for B_vec
    float *C_vec = (float *)malloc((N - 1) * sizeof(float)); // Allocate memory for C_vec

    if (!a || !b || !c || !d || !A_vec || !B_vec || !C_vec) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    calc_h(&h, N);
    calc_A(a, h, N);
    calc_B(b, h, N);
    calc_C(c, h, N);
    calc_D(d, h, N);

    calc_dig_A(A_vec, a, N, c, Neumann);
    calc_dig_B(B_vec, b, N, c, Neumann);
    calc_dig_C(C_vec, c, N, a, Neumann);
    tridiag(a, b, c, d, &u, 0, 1);

    

    // Remember to free allocated memory
    free(a);
    free(b);
    free(c);
    free(d);
    free(A_vec);
    free(B_vec);
    free(C_vec);

    return 0;
}
    
