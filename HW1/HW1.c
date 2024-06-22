// git link: https://github.com/roibaruch24/CFD_086376
// Including all the relevent headers:
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
// Makes a 2D arry to a vector
int offset2d(int i, int j, int ni)
{
    return j * ni + i;
}

// Implements the X vector for the NACA airfoil, Jmin
void naca_x_vec (double *x, int ni,double dx){
    for (int i = 11; i<=25; i++){
        x[offset2d(i, 1, ni)] = 1 - cos(0.5 * 3.1415 * dx * (26-i));
    }
    for(int i = 25; i<=39; i++){
        x[offset2d(i, 1, ni)] = 1 - cos(0.5*3.1415*dx*(26-i));
    }
}

//Implements the Y vector for the lower end of the NACA airfoil, Jmin
void naca_y_lower_vec (double *x,double *y, int ni, double t){
    double x_int = 1.008930411365;
    for(int i = 25; i<=39; i++){
        y[offset2d(i, 1, ni)] = -5*t*0.2969*sqrt(x_int*x[offset2d(i, 1, ni)])-
                                0.1260 * pow( (x_int * x[offset2d(i, 1, ni)]),1) -
                                0.3516 * pow( (x_int * x[offset2d(i, 1, ni)]),2) +
                                0.2843 * pow( (x_int * x[offset2d(i, 1, ni)]),3) -
                                0.1015 * pow( (x_int * x[offset2d(i, 1, ni)]),4);
    }
    
}
//Implements the Y vector for the upper end of the NACA airfoil, Jmin
void naca_y_upper_vec (double *x,double *y, int ni, double t)
{
    double x_int = 1.008930411365;
    for(int i = 11; i<=25; i++){
        y[offset2d(i, 1, ni)] = 5*t*0.2969*sqrt(x_int*x[offset2d(i, 1, ni)])-
        0.1260 * pow( (x_int * x[offset2d(i, 1, ni)]),1) -
        0.3516 * pow( (x_int * x[offset2d(i, 1, ni)]),2) +
        0.2843 * pow( (x_int * x[offset2d(i, 1, ni)]),3) -
        0.1015 * pow( (x_int * x[offset2d(i, 1, ni)]),4);

    }
    
}
// Implements the X vec for the upper part of the wake, Jmin
void x_wake_upper (double *x, int ni, double xsf){
    for (int i = 40; i<=50; i++){
         x[offset2d(i, 1, ni)] =  x[offset2d(i-1, 1, ni)] +xsf*( x[offset2d(i-1, 1, ni)]- x[offset2d(i-2, 1, ni)]);
    }
}
// Implements the X vec for the lower part of the wake, Jmin
void x_wake_lower (double *x, int ni){
    for (int i = 0; i<=10; i++){
         x[offset2d(i, 1, ni)] = x[offset2d(52-i, 1, ni)];
    }
}
// Implements the y vec for the TE upwards, Imax
void y_i_max (double *y, int ni, double ysf){
    y[offset2d(51, 1, ni)] = 0.02;
    for (int j = 2; j<=25; j++){
        y[offset2d(51, j, ni)] = y[offset2d(51, j-1, ni)]+ysf*(y[offset2d(51, j-1, ni)]-y[offset2d(51, j-2, ni)]);
    }
}
// Implements the y vec for the TE downwards, Imin
void y_i_min (double *y, int ni){
    y[offset2d(51, 1, ni)] = 0.02;
    for (int j = 2; j<=25; j++){
        y[offset2d(1, j, ni)] = -y[offset2d(51, j, ni)];
    }
}
// Implements the x vec for the TE upwards, Imax
void x_i_max (double *x, int ni){
    for (int j = 2; j<=25; j++){
        x[offset2d(51, j, ni)] =  x[offset2d(51, j, ni)] +  x[offset2d(51, 1, ni)];
    }
}
// Implements the x vec for the TE downwards, Imin
void x_i_min (double *x, int ni){
    for (int j = 2; j<=25; j++){
        x[offset2d(1, j, ni)] = x[offset2d(1, j, ni)] +  x[offset2d(51, 1, ni)];
    }
}

int main()
{
    int ni = 50; 
    int nj = 25;
    double xsf = 1.15;
    double ysf = 1.15;
    double dx = 0.0714;
    double t = 0.12;
    double *x = (double *)malloc(ni * nj * sizeof(double));
    double *y = (double *)malloc(ni * nj * sizeof(double));
    
    if (!x || !y) {
        printf("Memory allocation failed\n");
        return 1;
    }
    
    naca_x_vec(x, ni, dx);
    naca_y_lower_vec(x, y, ni, t);
    naca_y_upper_vec(x, y, ni, t);
    x_wake_upper(x, ni, xsf);
    x_wake_lower(x, ni);
    y_i_max(y, ni, ysf);
    y_i_min(y, ni);
    x_i_max(x, ni);
    x_i_min(x, ni);

    // Printing some values for verification
    for (int i = 0; i < ni; i++) {
        printf("x[%d] = %f, y[%d] = %f\n", i, x[offset2d(i, 1, ni)], i, y[offset2d(i, 1, ni)]);
    }

    free(x);
    free(y);

    return 0;
}

