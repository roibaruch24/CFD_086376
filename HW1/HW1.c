#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#define PI  3.14159265359

// Makes a 2D array to a vector
int offset2d(int i, int j, double ni) {
    return j * ni + i;
}
int input(double *t,int *imax,int *jmax,int *tel,int *le,int *teu,double *dy,double *xsf,double *ysf)
{
    const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW1\\input.txt";
    FILE *input = fopen(input_file_path, "rt");
    if (input == NULL)
    {
        perror("Error opening file");
        return -1; // Return an error code
    }
    fscanf(input,"%lf %d %d %d %d %d %lf %lf %lf",t,imax,jmax,tel,le,teu,dy,xsf,ysf); //reading the input
    fclose(input);
    return 0;
}
int print_output(int  ni,int  nj,double const *x,double const *y)
{
    const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW1\\output.txt";
    FILE *output = fopen(input_file_path, "wt"); //crating output file
    for (int i = 0; i <= ni*nj; i++)
    {
        fprintf(output, "%f %f\n",x[i],y[i]); //printing values into the file
    }
    fclose(output);
    return 0;
}

double airfoil(double *x, int ni, int i, int up_down, double t) {
    double x_int = 1.008930411365;
    double x_val = x[offset2d(i, 0, ni)];
    return up_down * 5 * t * (0.2969 * sqrt(x_int * x_val) -
                              0.1260 * x_int * x_val -
                              0.3516 * pow(x_int * x_val, 2) +
                              0.2843 * pow(x_int * x_val, 3) -
                              0.1015 * pow(x_int * x_val, 4));
}

void j_min(double *x, double *y, int ni, double dx, double t, int tel, int le, int teu, double xsf) {
    int k = 1;
    for (int i = tel; i <= le; i++) {
        x[offset2d(i, 0, ni)] = 1 - cos(0.5 * 3.14159265359 * dx * (le - i));
        y[offset2d(i, 0, ni)] = airfoil(x, ni, i, -1, t);
    }
    for (int i = le + 1; i <= teu; i++) {
        x[offset2d(i, 0, ni)] = x[offset2d(le - k, 0, ni)];
        y[offset2d(i, 0, ni)] = airfoil(x, ni, le - k, 1, t);
        k++;
    }
    for (int i = teu + 1; i < ni; i++) {
        x[offset2d(i, 0, ni)] = x[offset2d(i - 1, 0, ni)] + (x[offset2d(i - 1, 0, ni)] - x[offset2d(i - 2, 0, ni)]) * xsf;
        y[offset2d(i, 0, ni)] = 0;
    }
    for (int i = 0; i < tel; i++) {
        x[offset2d(i, 0, ni)] = x[offset2d(ni - i - 1, 0, ni)];
        y[offset2d(i, 0, ni)] = 0;
    }
}
void i_max(double *y, double *x, double ni, double ysf, int le, double dy) {
    y[offset2d(50, 0, ni)] = 0;
    y[offset2d(50, 1, ni)] = dy;
    //x[offset2d(51, 0, ni)] = 0;
    for (int j = 2; j <= 25; j++) {
        y[offset2d(50, j, ni)] = y[offset2d(50, j - 1, ni)] + ysf * (y[offset2d(50, j - 1, ni)] - y[offset2d(50, j - 2, ni)]);
    }
    for (int j = 1; j <= le; j++) {
        x[offset2d(50, j, ni)] = x[offset2d(50, 0, ni)];
    }
}

void i_min(double *y, double *x, double ni, int le, int dy) {
    for (int j = 1; j <= 25; j++) {
        y[offset2d(0, j, ni)] = -y[offset2d(50, j, ni)];
        x[offset2d(0, j, ni)] = x[offset2d(0, j, ni)] + x[offset2d(50, 0, ni)];
    }
}
void j_max(double *x,double *y,double const jmax,int const ni)
{
    int i_line_final = 0;
    double const R = y[offset2d(50,jmax,ni)];
    double const L = 0.5*3.14159265359*R + x[offset2d(50,jmax,ni)];
    double const delta = L/25;
    for (int i=0; i<=25; i++)
    {
        x[offset2d(i,jmax,ni)] = x[offset2d(0,jmax,ni)]-i*delta;
        y[offset2d(i,jmax,ni)] = y[offset2d(0,jmax,ni)];
        if (x[offset2d(i,jmax,ni)]<0)
        {
            i_line_final = i-1;
            break;
        }
    }
    double const theta = 0.5*3.14159265359+atan(x[offset2d(i_line_final,jmax,ni)]/y[offset2d(i_line_final,jmax,ni)]);
    double const d_theta = theta/(25-(i_line_final+1));
    for (int i = 25;i>=i_line_final;i--)
    {
        x[offset2d(i,jmax,ni)] = R*cos(3.14159265359+d_theta*(25-i));
        y[offset2d(i,jmax,ni)] = R*sin(3.14159265359+d_theta*(25-i));
    }
    int k = 1;
    for (int i = 26;i<=50;i++)
    {
        x[offset2d(i,jmax,ni)] = x[offset2d(25-k,jmax,ni)];
        y[offset2d(i,jmax,ni)] = -y[offset2d(25-k,jmax,ni)];
        k++;
    }
}
void interp_j(double *x,double *y,int const i,int const j,int const jmax,int const ni)
{
    double  mx = (x[offset2d(i,jmax,ni)]-x[offset2d(i,0,ni)])/jmax;
    double  my = (y[offset2d(i,jmax,ni)]-y[offset2d(i,0,ni)])/jmax;
    x[offset2d(i,j,ni)] = mx*j+x[offset2d(i,0,ni)];
    y[offset2d(i,j,ni)] = my*j+y[offset2d(i,0,ni)];
}
double x_ksi (double *x, int ni, int i, int j){
    return ((x[offset2d(i+1,j,ni)]-x[offset2d(i-1,j,ni)])/2);
}
double x_etta (double *x, int ni, int i, int j){
    return ((x[offset2d(i,j+1,ni)]-x[offset2d(i,j-1,ni)])/2);
}
double y_ksi (double *y, int ni, int i, int j){
    return ((y[offset2d(i+1,j,ni)]-y[offset2d(i-1,j,ni)])/2);
}
double y_etta (double *y, int ni, int i, int j){
    return ((y[offset2d(i,j+1,ni)]-y[offset2d(i,j-1,ni)])/2);
}


double x_ksi_ksi (double *x, int ni, int i, int j){
    return (x[offset2d(i+1,j,ni)]-2*x[offset2d(i,j,ni)]+x[offset2d(i-1,j,ni)]);
}
double x_etta_etta (double *x, int ni, int i, int j){
    return (x[offset2d(i,j+1,ni)]-2*x[offset2d(i,j,ni)]+x[offset2d(i,j-1,ni)]);
}
double y_ksi_ksi (double *y, int ni, int i, int j){
    return (y[offset2d(i+1,j,ni)]-2*y[offset2d(i,j,ni)]+y[offset2d(i-1,j,ni)]);
}
double y_etta_etta (double *y, int ni, int i, int j){
    return (y[offset2d(i,j+1,ni)]-2*y[offset2d(i,j,ni)]+y[offset2d(i,j-1,ni)]);
}

double calc_alpha (double *x ,double *y, int ni, int i, int j){
    double x_val_ksi = x_ksi(x, ni, i, j);
    double y_val_ksi = y_ksi(x, ni, i, j);
    return (pow(x_val_ksi,2)+pow(y_val_ksi,2));
}
double calc_beta (double *x ,double *y, int ni, int i, int j){
    double x_val_ksi = x_ksi(x, ni, i, j);
    double y_val_ksi = y_ksi(x, ni, i, j);
    double x_val_etta = x_etta(x, ni, i, j);
    double y_val_etta = y_etta(x, ni, i, j);
    return (x_val_etta*x_val_ksi+y_val_etta*y_val_ksi);
}
double calc_gamma (double *x ,double *y, int ni, int i, int j){
    double x_val_etta = x_etta(x, ni, i, j);
    double y_val_etta = y_etta(x, ni, i, j);
    return (pow(x_val_etta,2)+pow(y_val_etta,2));
}   

void Calc_A (double *x, double *y, int ni,int sweep, int imax, double *A, int j, int jmax, int i){
    if (sweep  == 1){
        for (int i = 0; i <= imax; i++){
            A[i] = -1*calc_alpha(x, y, ni, i, j);
        }
    }
    else{
        for (int j = 0; j <= jmax; j++)
        {
            A[j] = -1*calc_gamma(x, y, ni, i, j);
        }
    }

}
void Calc_B (double *x, double *y, int ni,int sweep, int imax, double *B, int j, int jmax, int i, double r){
    if (sweep  == 1){
        for (int i = 0; i <= imax; i++){
            B[i] = r + 2*calc_alpha(x, y, ni, i, j);
        }
    }
    else{
        for (int j = 0; j <= jmax; j++)
        {
             B[j] = r + 2*calc_gamma(x, y, ni, i, j);
        }
    }
}
void Calc_C (double *x, double *y, int ni,int sweep, int imax, double *C, int j, int jmax, int i){
    if (sweep  == 1){
        for (int i = 0; i <= imax; i++){
            C[i] = -1*calc_alpha(x, y, ni, i, j);
        }
    }
    else{
        for (int j = 0; j <= jmax; j++)
        {
            C[j] = -1*calc_gamma(x, y, ni, i, j);
        }
    }

}




int main() {
    double t;int imax;int jmax;int tel;int le;int teu;double dy;double xsf;double ysf;
    double dx = 0.0714;
    input(&t,&imax,&jmax,&tel,&le,&teu,&dy,&xsf,&ysf);
    int const ni = imax+1;int const nj = jmax+1;
    double *x =malloc((ni*nj) *sizeof(double));
    double *y =malloc((ni*nj) *sizeof(double));
    double *A =malloc((max(ni,nj)) *sizeof(double));
    double *B =malloc((max(ni,nj)) *sizeof(double));
    double *C =malloc((max(ni,nj)) *sizeof(double));
    double *D =malloc((max(ni,nj)) *sizeof(double));
    // Check if memory allocation succeeded
    if (x == NULL || y == NULL)
    {
        fprintf(stderr, "Memory allocation failed. Exiting...\n");
        return EXIT_FAILURE;
    }

    // Call the functions
    j_min(x, y, ni, dx, t, tel, le, teu, xsf);
    i_max(y, x, ni, ysf, le, dy);
    i_min(y, x, ni, le, dy);
    j_max(x, y, jmax, ni);
    
    
    for (int i = 1; i<=imax-1; i++) {
        for(int j = 1; j<= jmax-1; j++){
            interp_j(x, y, i, j, jmax,ni);
        }
    }

    print_output(ni,nj,x,y);
    free(x);
    free(y);

    return 0;
}
