#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI  3.14159265359
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Makes a 2D array to a vector
int offset2d(int i, int j, double ni) {
    return j * ni + i;
}
int find_fabs_Max(double *vector, double size) {
    double max = fabs(vector[0]);  // Assume the first element is the maximum
    for (int i = 1; i < size; i++) {
        if (fabs(vector[i]) > max) {
            max = fabs(vector[i]);
        }
    }
    return max;
}

int input(double *t,int *imax,int *jmax,int *tel,int *le,int *teu,double *dy,double *xsf,double *ysf)
{
     char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW1\\input.txt";
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

int print_output(int ni, int nj, double *x, double *y) {
    char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW1\\output.txt";
    FILE *output = fopen(input_file_path, "wt"); // creating output file

    if (output == NULL) {
        perror("Error opening file");
        return 1;
    }

    for (int i = 0; i < ni * nj; i++) {
        fprintf(output, "%f %f\n", x[i], y[i]); // printing values into the file
    }

    fclose(output); // closing the file
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
void j_max(double *x,double *y,double  jmax,int  ni)
{
    int i_line_final = 0;
    double  R = y[offset2d(50,jmax,ni)];
    double  L = 0.5*3.14159265359*R + x[offset2d(50,jmax,ni)];
    double  delta = L/25;
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
    double  theta = 0.5*3.14159265359+atan(x[offset2d(i_line_final,jmax,ni)]/y[offset2d(i_line_final,jmax,ni)]);
    double  d_theta = theta/(25-(i_line_final+1));
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
void interp_j(double *y,int i,int j,int jmax,int ni)
{
    double  my = (y[offset2d(i,jmax,ni)]-y[offset2d(i,0,ni)])/jmax;
    y[offset2d(i,j,ni)] = my*j+y[offset2d(i,0,ni)];
}
void interp_i(double *x,int i,int j,int imax,int  ni)
{
    double m_x = (x[offset2d(imax,j,ni)]-x[offset2d(0,j,ni)])/imax;
    x[offset2d(i,j,ni)] = m_x*i+x[offset2d(0,j,ni)];
}
void dx_ksi (double *x,double *x_ksi, int ni, int i, int j){
    x_ksi[offset2d(i,j,ni)] = ((x[offset2d(i+1,j,ni)]-x[offset2d(i-1,j,ni)])/2);
}
void dx_etta (double *x,double *x_etta, int ni, int i, int j){
    x_etta[offset2d(i,j,ni)] = ((x[offset2d(i,j+1,ni)]-x[offset2d(i,j-1,ni)])/2);
}
void dy_ksi (double *y,double *y_ksi, int ni, int i, int j){
     y_ksi[offset2d(i,j,ni)] = ((y[offset2d(i+1,j,ni)]-y[offset2d(i-1,j,ni)])/2);
}
void dy_etta (double *y, double *y_etta, int ni, int i, int j){
      y_etta[offset2d(i,j,ni)] = ((y[offset2d(i,j+1,ni)]-y[offset2d(i,j-1,ni)])/2);
}
void dx_ksi_ksi (double *x,double *x_ksi_ksi, int ni, int i, int j){
    x_ksi_ksi[offset2d(i,j,ni)] = (x[offset2d(i+1,j,ni)]-2*x[offset2d(i,j,ni)]+x[offset2d(i-1,j,ni)]);
}
void dx_etta_etta (double *x, double *x_etta_etta, int ni, int i, int j){
     x_etta_etta[offset2d(i,j,ni)] = (x[offset2d(i,j+1,ni)]-2*x[offset2d(i,j,ni)]+x[offset2d(i,j-1,ni)]);
}
void dy_ksi_ksi (double *y, double *y_ksi_ksi, int ni, int i, int j){
     y_ksi_ksi[offset2d(i,j,ni)] = (y[offset2d(i+1,j,ni)]-2*y[offset2d(i,j,ni)]+y[offset2d(i-1,j,ni)]);
}
void dy_etta_etta (double *y, double *y_etta_etta, int ni, int i, int j){
     y_etta_etta[offset2d(i,j,ni)] = (y[offset2d(i,j+1,ni)]-2*y[offset2d(i,j,ni)]+y[offset2d(i,j-1,ni)]);
}
void calc_derivatives(double *x_ksi, double *x_etta, double *y_ksi, double *y_etta, double *x_ksi_ksi, double *x_etta_etta,double *y_ksi_ksi, double *y_etta_etta, double *x, double *y, int ni, int imax, int jmax){
    for (int j = 0; j <= jmax; j++){
        for (int i = 1; i < imax; i++){
            dx_ksi(x, x_ksi, ni, i, j);
            dy_ksi(y, y_ksi, ni, i, j);
            dx_ksi_ksi(x, x_ksi_ksi, ni, i, j);
            dy_ksi_ksi(y, y_ksi_ksi, ni, i, j);    
        }
    }

    for (int i = 0; i <= imax; i++){
        for (int j = 1; j < jmax; j++){
            dx_etta(x, x_etta, ni, i, j);
            dy_etta(y, y_etta, ni, i, j);
            dx_etta_etta(x, x_etta_etta, ni, i, j);
            dy_etta_etta(y, y_etta_etta, ni, i, j);       
        }
    }

}
void calc_alpha (double *x_etta, double *y_etta,double *alpha,  int ni, int i, int j){
    alpha[offset2d(i,j,ni)] = pow(y_etta[offset2d(i,j,ni)],2)+pow(x_etta[offset2d(i,j,ni)],2);}
void calc_beta (double *x_etta, double *y_etta, double *x_ksi, double *y_ksi, double *beta,  int ni, int i, int j){
    beta[offset2d(i,j,ni)] = x_etta[offset2d(i,j,ni)]*x_ksi[offset2d(i,j,ni)] + y_etta[offset2d(i,j,ni)]*y_ksi[offset2d(i,j,ni)];
}
void calc_gamma (double *x_ksi, double *y_ksi,double *gamma,  int ni, int i, int j){
    gamma[offset2d(i,j,ni)] = pow(y_ksi[offset2d(i,j,ni)],2)+pow(x_ksi[offset2d(i,j,ni)],2);
}   
void calc_greek(double *x_etta, double *y_etta, double *x_ksi, double *y_ksi,double *alpha, double *beta, double *gamma, int ni, int imax, int jmax){
    for(int i = 1; i < imax; i++){
        for (int j = 1; j < jmax; j++){
             calc_alpha (x_etta, y_etta,alpha, ni, i, j);
             calc_beta (x_etta, y_etta, x_ksi, y_ksi, beta, ni, i, j);
             calc_gamma (x_ksi, y_ksi, gamma, ni, i, j);
        }
    }
}
void Calc_A (double *alpha, double *gamma,  int ni,int sweep, int imax, double *A, int j, int jmax, int i){
    if (sweep  == 1){
        for (int i = 0; i <= imax; i++){
            A[i] = -alpha[offset2d(i,j,ni)];
        }
    }
    else{
        for (int j = 0; j <= jmax; j++)
        {
            A[j] = -1*gamma[offset2d(i,j,ni)];
        }
    }
}
void Calc_B (double *alpha, double *gamma, int ni,int sweep, int imax, double *B, int j, int jmax, int i, double r){
    if (sweep  == 1){
        for (int i = 0; i <= imax; i++){
            B[i] = r + 2*alpha[offset2d(i,j,ni)];
        }
    }
    else{
        for (int j = 0; j <= jmax; j++)
        {
             B[j] = r + 2*gamma[offset2d(i,j,ni)];
        }
    }
}
void Calc_C (double *alpha, double *gamma,  int ni, int sweep, int imax, double *C, int j, int jmax, int i){
    if (sweep  == 1){
        for (int i = 0; i <= imax; i++){
            C[i] = -1*alpha[offset2d(i,j,ni)];
        }
    }
    else{
        for (int j = 0; j <= jmax; j++)
        {
            C[j] = -1*gamma[offset2d(i,j,ni)];
        }
    }
}
void LHS (double *alpha, double *gamma, int ni, int sweep, int imax,double *A, double *B, double *C, int jmax, double r, int i, int j){
    
    Calc_A (alpha, gamma, ni, sweep, imax, A, j, jmax,i);
    Calc_B (alpha, gamma, ni, sweep, imax, B, j, jmax, i, r);
    Calc_C (alpha, gamma, ni, sweep, imax, C, j, jmax,i);
    A[0] = 0;
    B[0] = 1;
    C[0] = 0;
   
    if (sweep == 1){
        A[imax] = 0;
        B[imax] = 1;
        C[imax] = 0;
        
    }
    else {
        A[imax] = 0;
        B[imax] = 1;
        C[imax] = 0;
        
    }

}
void calc_psi_boundary (double *x, double *y,double *x_etta, double *y_etta, double *x_etta_etta, double *y_etta_etta, int ni, double *psi, int jmax, int imax){
    for (int i = 0; i<= imax; i+= imax){ // running for either i = 0 or i = imax
        for (int j = 0; j <= jmax; j++){
            if (fabs(y_etta[offset2d(i,j,ni)]) > fabs(x_etta[offset2d(i,j,ni)])){
                psi[offset2d(i,j,ni)] = -y_etta_etta[offset2d(i,j,ni)]/y_etta[offset2d(i,j,ni)];
            }
            else{
                psi[offset2d(i,j,ni)] = -x_etta_etta[offset2d(i,j,ni)]/x_etta[offset2d(i,j,ni)];
            }
        }
    }
}
void calc_phi_boundary (double *x, double *y, double *x_ksi, double *y_ksi, double *x_ksi_ksi, double *y_ksi_ksi, int ni, double *phi, int jmax, int imax){
    for (int j =0; j<=jmax; j+=jmax) //for i max after finishing i min
    {
        for (int i=0; i<=imax; i++) //for i min
        {
            if (fabs(y_ksi[offset2d(i,j,ni)])>fabs(x_ksi[offset2d(i,j,ni)]))
            {
                phi[offset2d(i,j,ni)] = -y_ksi_ksi[offset2d(i,j,ni)]/y_ksi[offset2d(i,j,ni)];
            }
            else
            {
                phi[offset2d(i,j,ni)] = -x_ksi_ksi[offset2d(i,j,ni)]/x_ksi[offset2d(i,j,ni)];
            }
        }
    }
}
void control_function_inside (double *x, double *y, int imax, int jmax, int ni, int nj,  double *psi,double *x_etta, double *y_etta, double *x_ksi, double *y_ksi, double *x_etta_etta, double *y_etta_etta, double *x_ksi_ksi, double *y_ksi_ksi, double *phi , double *A, double *B, double *C, double *D, int control_case ){
    if (control_case==1)
    {
        for (int i=1; i<imax; i++)
        {
            for (int j=1; j<jmax; j++)
            {
                psi[offset2d(i,j,ni)]=0;
                phi[offset2d(i,j,ni)]=0;
            }
        }
    }
    else if (control_case==2)
    {
        for (int j=1; j<jmax; j++)
        {
            for (int i=1; i<imax; i++)
            {
                phi[offset2d(i,j,ni)]=0;
                interp_j(psi,i,j,imax,ni);
            }
        }
    }
    else if (control_case==3)
    {
        for (int i=1; i<imax; i++)
        {
            for (int j=1; j<jmax; j++)
            {
                psi[offset2d(i,j,ni)]=0;
                interp_j(phi,i,j,jmax,ni);
            }
        }
    }
    else if (control_case==4)
    {
        for (int j=1; j<jmax; j++)
        {
            for (int i=1; i<imax; i++)
            {
                interp_i(psi,i,j,imax,ni);
            }
        }
        for (int i=1; i<imax; i++)
        {
            for (int j=1; j<jmax; j++)
            {
                interp_j(phi,i,j,jmax,ni);
            }
        }
    }
}
void control_function(double *x, double *y, int imax, int jmax, int ni, int nj,  double *psi,double *x_etta, double *y_etta, double *x_ksi, double *y_ksi, double *x_etta_etta, double *y_etta_etta, double *x_ksi_ksi, double *y_ksi_ksi, double *phi , double *A, double *B, double *C, double *D, int control_case ){
    calc_psi_boundary(x, y, x_etta, y_etta, x_etta_etta, y_etta_etta, ni, psi, jmax, imax);
    calc_phi_boundary(x, y, x_ksi, y_ksi, x_ksi_ksi, y_ksi_ksi, ni, phi, jmax, imax);
    control_function_inside (x, y, imax, jmax, ni, nj,  psi, x_etta, y_etta, x_ksi, y_ksi, x_etta_etta, y_etta_etta, x_ksi_ksi, y_ksi_ksi, phi , A, B, C, D, control_case );
}
void L_operator(double *alpha, double *beta, double *gamma, double *var, double *L, int ni, double *psi, double *phi, int imax, int jmax){
    for (int i = 1; i< imax; i++){
        for (int j = 1; j < jmax; j++){
            double x_i = var[offset2d(i,j,ni)];
            double x_i_1 = var[offset2d(i+1,j,ni)];
            double x_i_m1 = var[offset2d(i-1,j,ni)];
            double x_j_1 = var[offset2d(i,j+1,ni)];
            double x_j_m1 = var[offset2d(i,j-1,ni)];
            double x_i_1_j_1 = var[offset2d(i+1,j+1,ni)];
            double x_i_1_j_m1 = var[offset2d(i+1,j-1,ni)];
            double x_i_m1_j_1 = var[offset2d(i-1,j+1,ni)];
            double x_i_m1_j_m1 = var[offset2d(i-1,j-1,ni)];
            L[offset2d(i,j,ni)] = alpha[offset2d(i,j,ni)]*((x_i_1-2*x_i+x_i_m1)+0.5*phi[offset2d(i,j,ni)]*(x_i_1-x_i_m1))
            -0.5*beta[offset2d(i,j,ni)]*(x_i_1_j_1-x_i_1_j_m1-x_i_m1_j_1+x_i_m1_j_m1)+gamma[offset2d(i,j,ni)]
            *(x_j_1-2*x_i+x_j_m1+0.5*psi[offset2d(i,j,ni)]*(x_j_1-x_j_m1));
        }
    }
}
void RHS (double *D, int i, int j, int ni, double r, double w, double *L, int sweep, int imax, int jmax){
    if (sweep == 1){
        for (int i = 1; i < imax; i++){
        D[i] = r*w*L[offset2d(i,j,ni)]; 
        }
        D[0] = 0;
        D[imax] = 0;
    }
    else{
        for (int j = 0; j <= jmax; j++){
            D[j] = L[offset2d(i,j,ni)]; // L is the solutions vector!!!!!!!!
        }
        D[0] = 0;
        D[imax] = 0;
    }

}
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
void step (double *x, double *y,int jmax, double *alpha, double *beta, double *gamma, double *psi, double *phi, int ni, int imax, double *A , double *B, double *C, double *D, double r, double w, double *L_x, double *L_y, double *f_x, double *f_y, double *solutions_x, double *solutions_y, double *C_x, double *C_y, double *Del_x, double *Del_y, double *x_ksi, double *y_ksi, double *x_etta, double *y_etta, double *x_ksi_ksi, double *y_ksi_ksi, double *x_etta_etta, double *y_etta_etta ){
     // Sweep 1
    int sweep = 1;
    int is = 1; int ie = imax -1;
    for (int j = 1; j < jmax; j++){
        int i = 0;
        // X values
        LHS(alpha, gamma, ni, sweep, imax, A, B, C, jmax, r, i, j);
        RHS(D, i, j, ni, r, w, L_x, sweep, imax,jmax);
        tridiag(A, B, C, D, f_x, is, ie);
        for (int k  = 1; k < imax; k++ ){
            solutions_x[offset2d(k,j,ni)] = f_x[k];    
        }
        // Y values
        LHS(alpha, gamma, ni, sweep, imax, A, B, C, jmax, r, i, j);
        RHS(D, i, j, ni, r, w, L_y, sweep, imax,jmax);
        tridiag(A, B, C, D, f_y, is, ie);
        for (int k  = 1; k < imax; k++ ){
            solutions_y[offset2d(k,j,ni)] = f_y[k];
        }           
    }



    // Sweep 2
    sweep = 2;
     is = 1;  ie = jmax -1;
    for (int i = 1; i < imax; i++){
        int j = 0;
        // X values
        LHS(alpha, gamma, ni, sweep, imax, A, B, C, jmax, r, i, j);
        RHS(D, i, j, ni, r, w, solutions_x, sweep, imax,jmax);
        tridiag(A, B, C, D, C_x, is, ie);
        for (int k  = 1; k < jmax; k++ ){
            Del_x[offset2d(i,k,ni)] = C_x[k];    
        }
        // Y values
        LHS(alpha, gamma, ni, sweep, imax, A, B, C, jmax, r, i, j);
        RHS(D, i, j, ni, r, w, solutions_y, sweep, imax,jmax);
        tridiag(A, B, C, D, C_y, is, ie);
        for (int k  = 1; k < jmax; k++ ){
            Del_y[offset2d(i,k,ni)] = C_y[k];
            
        }           
    }
    // Updating both X and Y values
    for(int i = 1; i < imax; i++){
        for(int j = 1; j < jmax; j++){
            x[offset2d(i,j,ni)] += Del_x[offset2d(i,j,ni)];
            y[offset2d(i,j,ni)] += Del_y[offset2d(i,j,ni)];
        }
    }
    calc_derivatives(x_ksi, x_etta, y_ksi, y_etta, x_ksi_ksi, x_etta_etta, y_ksi_ksi, y_etta_etta, x, y, ni, imax, jmax);
    calc_greek(x_etta, y_etta, x_ksi,y_ksi, alpha, beta, gamma, ni, imax, jmax );
    L_operator(alpha, beta, gamma, x, L_x, ni, psi, phi, imax, jmax);
    L_operator(alpha, beta, gamma, y, L_y, ni, psi, phi, imax, jmax);   
}

  
int main() {
    double t;int imax;int jmax;int tel;int le;int teu;double dy;double xsf;double ysf;
    double r = 0.005;
    double w = 1.2;
    int control_case = 2;
    input(&t,&imax,&jmax,&tel,&le,&teu,&dy,&xsf,&ysf);
    double dx = (double)1/(le-tel);
    int  ni = imax+1;int  nj = jmax+1;
    double *x =malloc((ni*nj) *sizeof(double));
    double *y =malloc((ni*nj) *sizeof(double));
    double *psi =malloc((ni*nj) *sizeof(double));
    double *phi =malloc((ni*nj) *sizeof(double));
    double *x_ksi =malloc((ni*nj) *sizeof(double));
    double *x_ksi_ksi =malloc((ni*nj) *sizeof(double));
    double *x_etta =malloc((ni*nj) *sizeof(double));
    double *x_etta_etta =malloc((ni*nj) *sizeof(double));
    double *y_etta =malloc((ni*nj) *sizeof(double));
    double *y_etta_etta =malloc((ni*nj) *sizeof(double));
    double *y_ksi =malloc((ni*nj) *sizeof(double));
    double *y_ksi_ksi =malloc((ni*nj) *sizeof(double));
    double *alpha =malloc((ni*nj) *sizeof(double));
    double *beta =malloc((ni*nj) *sizeof(double));
    double *gamma =malloc((ni*nj) *sizeof(double));
    double *A =malloc((MAX(ni, nj)) *sizeof(double));
    double *B =malloc((MAX(ni, nj)) *sizeof(double));
    double *C =malloc((MAX(ni, nj)) *sizeof(double));
    double *D =malloc((MAX(ni, nj)) *sizeof(double));
    double *L_x =malloc((ni*nj) *sizeof(double));
    double *L_y =malloc((ni*nj) *sizeof(double));
    double *f_x =malloc((MAX(ni, nj)) *sizeof(double));
    double *f_y =malloc((MAX(ni, nj)) *sizeof(double));
    double *solutions_x =malloc((ni*nj) *sizeof(double));
    double *solutions_y =malloc((ni*nj) *sizeof(double)); 
    double *C_x =malloc((MAX(ni, nj)) *sizeof(double));
    double *C_y =malloc((MAX(ni, nj)) *sizeof(double)); 
    double *Del_x =malloc((ni*nj) *sizeof(double));
    double *Del_y =malloc((ni*nj) *sizeof(double)); 
    
    // Call the functions
    j_min(x, y, ni, dx, t, tel, le, teu, xsf);
    i_max(y, x, ni, ysf, le, dy);
    i_min(y, x, ni, le, dy);
    j_max(x, y, jmax, ni);
    for (int i = 1; i<=imax-1; i++) {
        for(int j = 1; j<= jmax-1; j++){
            interp_j(x, i, j, jmax,ni);
            interp_j(y, i, j, jmax,ni);
        }
    }
    
    calc_derivatives(x_ksi, x_etta, y_ksi, y_etta, x_ksi_ksi, x_etta_etta, y_ksi_ksi, y_etta_etta, x, y, ni, imax, jmax);
    calc_greek(x_etta, y_etta, x_ksi,y_ksi, alpha, beta, gamma, ni, imax, jmax );
    control_function(x, y, imax, jmax, ni, nj, psi, x_etta, y_etta, x_ksi, y_ksi, x_etta_etta, y_etta_etta, x_ksi_ksi, y_ksi_ksi, phi , A, B, C, D, control_case );
    L_operator(alpha, beta, gamma, x, L_x, ni, psi, phi, imax, jmax);
    L_operator(alpha, beta, gamma, y, L_y, ni, psi, phi, imax, jmax);
    print_output(ni, nj, phi, psi);
   
    
    int iteration = 1;
    while (log10(find_fabs_Max(L_x,ni*nj)) > -6 && log10(find_fabs_Max(L_y,ni*nj)) > -6)
    {
       step(x, y, jmax, alpha, beta, gamma, psi, phi, ni, imax, A , B, C, D, r, w, L_x, L_y, f_x, f_y, solutions_x, solutions_y, C_x, C_y, Del_x, Del_y, x_ksi, y_ksi, x_etta, y_etta, x_ksi_ksi, y_ksi_ksi, x_etta_etta, y_etta_etta); 
       iteration ++;
    }
    print_output(ni, nj, x, y);

    free(x);
    free(y);
    free(psi);
    free(phi);
    free(x_ksi);
    free(x_ksi_ksi);
    free(x_etta);
    free(x_etta_etta);
    free(y_etta);
    free(y_etta_etta);
    free(y_ksi);
    free(y_ksi_ksi);
    free(alpha);
    free(beta);
    free(gamma);
    free(A);
    free(B);
    free(C);
    free(D);
    free(L_x);
    free(L_y);
    free(solutions_x);
    free(solutions_y);
    free(C_x);
    free(C_y);
    free(Del_x);
    free(Del_y);



    return 0;
}
