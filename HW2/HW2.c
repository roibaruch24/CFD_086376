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
int input(int *ni, int *nj, int *teu, int *tel, double *M, double *alpha, double *p, double *rho) 
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

    // Read Mach number, alpha, pressure, and rho from parameters file
    if (fscanf(input_param, "%lf %lf %lf %lf", M, alpha, p, rho) != 4) 
    {
        printf("Error reading M, alpha, p, rho from the parameters file.\n");
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
void matrics(int imax, int jmax, int ni, int nj, double *x, double *y, double *jacobian, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y)
{
    
    double *y_ksi = (double *)malloc(ni*nj * sizeof(double));
    double *x_ksi = (double *)malloc(ni*nj * sizeof(double));
    double *x_etta = (double *)malloc(ni*nj * sizeof(double));
    double *y_etta = (double *)malloc(ni*nj * sizeof(double));
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
    for (int i  = 0; i <= imax; i++){
            for (int j = 0; j <= jmax; j++){
                jacobian[offset2d(i,j,ni)] = 1/(x_ksi[offset2d(i,j,ni)]*y_etta[offset2d(i,j,ni)] - y_ksi[offset2d(i,j,ni)]*x_etta[offset2d(i,j,ni)]);
            }
    }
    // inverting the derivitives
    for (int i = 0; i <= imax; i++){
        for( int j = 0; j <= jmax; j++){
            ksi_x[offset2d(i, j, ni)] = jacobian[offset2d(i, j, ni)]*y_etta[offset2d(i, j, ni)];
            ksi_y[offset2d(i, j, ni)] = -jacobian[offset2d(i, j, ni)]*x_etta[offset2d(i, j, ni)];
            etta_x[offset2d(i, j, ni)] = -jacobian[offset2d(i, j, ni)]*y_ksi[offset2d(i, j, ni)];
            etta_y[offset2d(i, j, ni)] = jacobian[offset2d(i, j, ni)]*x_ksi[offset2d(i, j, ni)];
        }
    }
free(y_ksi);
free(x_ksi);
free(x_etta);
free(y_etta);


}
void freestream (double rho, double p, double M, double alpha, double *Q, int imax, int jmax, int ni, int nj){
    double a = sqrt(gamma*rho/p);
    double u = M*a*cos(alpha/57.3);
    double v = M*a*sin(alpha/57.3);
    for(int i = 0; i<= imax; i++){
        for( int j=0; j<= jmax; j++){
            Q[offset3d(i, j, 0, ni, nj)] = rho;
            Q[offset3d(i, j, 1, ni, nj)] = rho*u;
            Q[offset3d(i, j, 2, ni, nj)] = rho*v;
            Q[offset3d(i, j, 3, ni, nj)] = p/(gamma-1) + 0.5*rho*(pow(u,2)+pow(v,2));
        }
    }

}
void bc_wall(double *x, double *y, int tel, int teu, double *Q, int ni, int nj, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y){
    for (int i = tel; i<=teu; i++){
            double u_1 = Q[offset3d(i, 1, 1,ni, nj)]/Q[offset3d(i, 0, 1, ni, nj)];
            double v_1 = Q[offset3d(i, 1, 2,ni, nj)]/Q[offset3d(i, 1, 0, ni, nj)];
            double u = ((ksi_x[offset2d(i, 1, ni)]*u_1 )+ksi_y[offset2d(i, 1, ni)]*v_1)/(ksi_x[offset2d(i, 0, ni)]-ksi_y[offset2d(i, 0, ni)]*(etta_x[offset2d(i, 0, ni)]/etta_y[offset2d(i, 0, ni)]));
            double v =(-etta_x[offset2d(i,0,ni)])/(etta_y[offset2d(i,0,ni)])*u;
            double rho = Q[offset3d(i,1,0,ni,nj)];
            double e1 = Q[offset3d(i,1,3,ni,nj)];
            double P = (gamma-1)*(e1-0.5*rho*(pow(u_1,2)+(pow(v_1,2))));
            double e = P/(gamma-1)+0.5*rho*(pow(u,2)+(pow(v,2)));
            Q[offset3d(i,0,0,ni,nj)] = rho;
            Q[offset3d(i,0,1,ni,nj)] = rho*u;
            Q[offset3d(i,0,2,ni,nj)] = rho*v;
            Q[offset3d(i,0,3,ni,nj)] = e;
    }

}
void Kutta (double *x, double *y, int tel, int teu, double *Q, int ni, int nj, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y){
    double u_1_tel = Q[offset3d(tel,1,1,ni,nj)]/Q[offset3d(tel,1,0,ni,nj)];
    double v_1_tel = Q[offset3d(tel,1,2,ni,nj)]/Q[offset3d(tel,1,0,ni,nj)];
    double u_tel = ((ksi_x[offset2d(tel, 1, ni)]*u_1_tel )+ksi_y[offset2d(tel, 1, ni)]*v_1_tel)/(ksi_x[offset2d(tel, 0, ni)]-ksi_y[offset2d(tel, 0, ni)]*(etta_x[offset2d(tel, 0, ni)]/etta_y[offset2d(tel, 0, ni)]));
    double v_tel =(-etta_x[offset2d(tel,0,ni)])/(etta_y[offset2d(tel,0,ni)])*u_tel;
    double rho_tel = Q[offset3d(tel,1,0,ni,nj)];
    double e1_tel = Q[offset3d(tel,1,3,ni,nj)];
    double P_tel = (gamma-1)*(e1_tel-0.5*rho_tel*(pow(u_1_tel,2)+(pow(v_1_tel,2))));

    double u_1_teu = Q[offset3d(teu,1,1,ni,nj)]/Q[offset3d(teu,1,0,ni,nj)];
    double v_1_teu = Q[offset3d(teu,1,2,ni,nj)]/Q[offset3d(teu,1,0,ni,nj)];
    double u_teu = ((ksi_x[offset2d(teu, 1, ni)]*u_1_teu )+ksi_y[offset2d(teu, 1, ni)]*v_1_teu)/(ksi_x[offset2d(teu, 0, ni)]-ksi_y[offset2d(teu, 0, ni)]*(etta_x[offset2d(teu, 0, ni)]/etta_y[offset2d(teu, 0, ni)]));
    double v_teu =(-etta_x[offset2d(teu,0,ni)])/(etta_y[offset2d(teu,0,ni)])*u_teu;
    double rho_teu = Q[offset3d(teu,1,0,ni,nj)];
    double e1_teu = Q[offset3d(teu,1,3,ni,nj)];
    double P_teu = (gamma-1)*(e1_teu-0.5*rho_teu*(pow(u_1_teu,2)+(pow(v_1_teu,2))));
    //averaging
    double P_te = 0.5*(P_tel+P_teu);
    double u_te = 0.5*(u_tel+u_teu);
    double v_te = 0.5*(v_tel+v_teu);
    double rho_te = 0.5*(rho_tel+rho_teu);
    double e_te = P_te/(gamma-1)+0.5*rho_te*(pow(u_te,2)+(pow(v_te,2)));
    for (int i=tel;i <= teu;i+=teu-tel)
    {
        Q[offset3d(i,0,0,ni,nj)] = rho_te;
        Q[offset3d(i,0,1,ni,nj)] = rho_te*u_te;
        Q[offset3d(i,0,2,ni,nj)] = rho_te*v_te;
        Q[offset3d(i,0,3,ni,nj)] = e_te;
    }
}
void bc_cut(double *x, double *y, int tel, int teu, double *Q, int ni, int nj, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y){
    for(int i = 1; i < tel; i++){
        for (int k = 0; k <= 3; k++){
            Q[offset3d(i, 0, k, ni, nj)] = 0.5*(Q[offset3d(i,1,k,ni,nj)] + Q[offset3d(ni-1, 1,k,ni, nj)]);
            Q[offset3d(ni - 1 - i, 0, k, ni, nj)] = Q[offset3d(i, 0 , k, ni, nj)];
        }
    }
}
void bc_outflow(double *x, double *y, int tel, int teu, double *Q, int ni, int nj, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y, int jmax){
    for(int j = 1; j < jmax; j++){
        for (int k = 0; k <= 3; k++){
            Q[offset3d(0, j, k, ni, nj)] =  Q[offset3d(1, j, k, ni, nj)];
            Q[offset3d(ni - 1, j, k, ni, nj)] = Q[offset3d(ni-2, j , k, ni, nj)];
        }
    }
}
void RHS(int ni, int nj, int imax, int jamx, double *Jacobian, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y, double *Q, double *W, double *S){   
    /* ksi direction */ 
    for (int j = 1; j < jamx; j++){
        for ( int i = 0; i <= imax; i++){
            double u = Q[offset3d(i, j, 1, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double v = Q[offset3d(i, j, 2, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double U_conv = ksi_x[offset2d(i, j, ni)]*u + ksi_x[offset2d(i, j, ni)]*v;
            double P = ((gamma-1)*Q[offset3d(i, j, 3, ni, nj)]-0.5*Q[offset3d(i, j, 0, ni, nj)]*(pow(u,2)+pow(v,2)))/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 0, ni)] = (Q[offset3d(i, j, 0, ni, nj)]*U_conv)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 1, ni)] = (Q[offset3d(i, j, 1, ni, nj)]*U_conv + ksi_x[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 2, ni)] = (Q[offset3d(i, j, 2, ni, nj)]*U_conv +ksi_y[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
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
            double V_conv = etta_x[offset2d(i, j, ni)]*u + etta_y[offset2d(i, j, ni)]*v;
            double P = ((gamma-1)*Q[offset3d(i, j, 3, ni, nj)]-0.5*Q[offset3d(i, j, 0, ni, nj)]*(pow(u,2)+pow(v,2)))/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 0, ni)] = (Q[offset3d(i, j, 0, ni, nj)]*V_conv)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 1, ni)] = (Q[offset3d(i, j, 1, ni, nj)]*V_conv + etta_x[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 2, ni)] = (Q[offset3d(i, j, 2, ni, nj)]*V_conv + etta_y[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 3, ni)] = V_conv*(P+Q[offset3d(i, j, 3, ni, nj)])/Jacobian[offset2d(i, j, ni)];
        }
        for(int j = 1; j < jamx ; j++){
            for( int k = 0; k < 4; k++ ){
                S[offset3d(i, j, k, ni, nj)] = -0.5*(W[offset2d(j+1,k, ni)]-W[offset2d(j-1,k,ni)]); 
            }
        }

    }

}
int smooth(double *q, double *s, double *jac, double *xx, double *xy, double *yx, double *yy, int id, int jd, double *s2, double *rspec, double *qv, double *dd, double epse, double fsmach, double dt){
double *rho, *u_vel, *v_vel, *t_e;

    double eratio, smool, gm1, ggm1, cx, cy, eps, ra, u, v, qq, ss, st, 
          qav, qxx, ssfs, qyy;
    int ib, ie, jb, je, i, j, offset, offsetp1, offsetm1, ip, ir, n,
        jp, jr;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    smool = 1.0;
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    ib = 1;
    ie = id - 1;
    jb = 1;
    je = jd - 1;

    cx = 2.;
    cy = 1.;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*     smoothing in xi direction */

    for (j = jb; j < je; j++) {
	for (i = 0; i < id; i++) {
	    offset = id * j + i;
	    eps = epse / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq);
            rspec[i] = eps * (fabs(xx[offset] * u + xy[offset] * v) + sqrt((xx[offset] * xx[offset] + xy[offset] * xy[offset]) * ss + 0.01));
	    qv[i] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (i = ib; i < ie; i++) {
	    ip = i + 1;
	    ir = i - 1;
	    qxx = qv[ip] - qv[i] * 2. + qv[ir];
	    qav = (qv[ip] + qv[i] * 2. + qv[ir]) * .25;
	    dd[i] = eratio * fabs(qxx / qav);
	}

	dd[0] = dd[1];
	dd[id - 1] = dd[id - 2];

	for (n = 0; n < 4; n++) {
	    for (i = ib; i < ie; i++) {
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j) * id + i + 1;
		offsetm1 = (jd * n + j) * id + i - 1;
		s2[i] = q[offsetp1] - 2.0 * q[offset] + q[offsetm1];
	    }

	    s2[0] = s2[1] * -1.;
	    s2[id - 1] = s2[id - 2] * -1.;

	    for (i = ib; i < ie; i++) {
		ip = i + 1;
		ir = i - 1;
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j) * id + i + 1;
		offsetm1 = (jd * n + j) * id + i - 1;
		st = ((dd[ip] + dd[i]) * .5 * (q[offsetp1] - q[offset]) - cx / (cx + dd[ip] + dd[i]) * (s2[ip] - s2[i])) * (rspec[ip] + rspec[i]) + ((dd[ir] + dd[i]) * .5 * (q[offsetm1] - q[offset]) - cx / (cx + dd[ir] + dd[i]) * (s2[ir] - s2[i])) * (rspec[ir] + rspec[i]);
		s[offset] += st * .5 * dt;
	    }
	}
    }

/*     smoothing in eta direction */

    ssfs = 1. / (0.001 + fsmach * fsmach);
    for (i = ib; i < ie; i++) {
	for (j = 0; j < jd; j++) {
	    offset = id * j + i;
	    eps = epse / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq) * (1.0 - smool) + smool * qq * ssfs;
            rspec[j] = eps * (fabs(yx[offset] * u + yy[offset] * v) + sqrt((yx[offset] * yx[offset] + yy[offset] * yy[offset]) * ss + 0.01));
	    qv[j] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (j = jb; j < je; j++) {
	    jp = j + 1;
	    jr = j - 1;
	    qyy = qv[jp] - qv[j] * 2. + qv[jr];
	    qav = (qv[jp] + qv[j] * 2. + qv[jr]) * .25;
	    dd[j] = eratio * fabs(qyy / qav);
	}

	dd[0] = dd[1];
	dd[jd - 1] = dd[jd - 2];

	for (n = 0; n < 4; n++) {
	    for (j = jb; j < je; j++) {
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j + 1) * id + i;
		offsetm1 = (jd * n + j - 1) * id + i;
		s2[j] = q[offsetp1] - 2.0 * q[offset] + q[offsetm1];
	    }

	    s2[0] = s2[1] * -1.;
	    s2[jd - 1] = s2[jd - 2] * -1.;

	    for (j = jb; j < je; j++) {
		jp = j + 1;
		jr = j - 1;
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j + 1) * id + i;
		offsetm1 = (jd * n + j - 1) * id + i;
		st = ((dd[jp] + dd[j]) * .5 * (q[offsetp1] - q[offset]) - cy / (cy + dd[jp] + dd[j]) * (s2[jp] - s2[j])) * (rspec[jp] + rspec[j]) + ((dd[jr] + dd[j]) * .5 * (q[offsetm1] - q[offset]) - cy / (cy + dd[jr] + dd[j]) * (s2[jr] - s2[j])) * (rspec[jr] + rspec[j]);
		s[offset] += st * .5 * dt;
	    }
	}
    }

    return 0;
} /* smooth */
void BC(double *x, double *y, int tel, int teu, double *Q, int ni, int nj, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y, int jmax){
    bc_wall(x, y, tel, teu, Q, ni, nj, ksi_x, ksi_y, etta_x, etta_y);
    Kutta(x, y, tel, teu, Q, ni, nj, ksi_x, ksi_y, etta_x, etta_y);
    bc_cut(x ,y , tel, teu, Q, ni, nj, ksi_x, ksi_y, etta_x, etta_y);
    bc_outflow (x, y, tel, teu, Q, ni, nj, ksi_x, ksi_y, etta_x, etta_y, jmax);
    
}
int print_output(int ni,int nj,double *x, double *y)
{
    //char *output_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\output.txt";
    char *output_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\output.txt";
    
    FILE *output= fopen(output_file_path, "wt"); //crating output file
    for (int i = 0; i < ni*nj; i++){
            fprintf(output, "%f %f\n", x[i], y[i]);
        }
    
    fclose(output);
    return 0;
}
int print_Q(int ni,int nj,double *Q)
{
    //char *output_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\output.txt";
    char *output_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\Q.txt";
    
    FILE *output= fopen(output_file_path, "wt"); //crating output file
    for (int i = 0; i < ni*nj*4; i++){
            fprintf(output, "%f\n", Q[i]);
        }
    
    fclose(output);
    return 0;
}

int main() {
    int ni, nj, teu, tel;
    double M_0, alpha, p_0, rho_0; // Ensure these are doubles
    // Read the parameters (ni, nj, teu, tel, Mach number, etc.)
    if (input(&ni, &nj, &teu, &tel, &M_0, &alpha, &p_0, &rho_0) != 0) {
        printf("Failed to read input parameters.\n");
        return -1; // Exit if input fails
    }
    printf("ni: %d, nj: %d, teu: %d, tel: %d\n", ni, nj, teu, tel);
    printf("Mach number: %f, Alpha: %f, Pressure: %f, rho: %f\n", M_0, alpha, p_0, rho_0);

    // Allocate memory for x and y in main
    int num_points = ni * nj;
    int jmax = nj - 1;
    int imax = ni - 1;
    double *x = (double *)malloc(num_points * sizeof(double));
    double *y = (double *)malloc(num_points * sizeof(double));
    double *jacobian = (double *)malloc(num_points * sizeof(double));
    double *p = (double *)malloc(num_points * sizeof(double));
    double *rho = (double *)malloc(num_points * sizeof(double));
    double *Q = (double *)malloc(num_points * 4 * sizeof(double));
    double *S = (double *)malloc(num_points * 4 * sizeof(double));
    double *W = (double *)malloc(MAX(ni, nj) * 4 * sizeof(double));
    //double *S = (double *)malloc(MAX(ni,nj) * 4 * sizeof(double));
    double *ksi_x = (double *)malloc(num_points * sizeof(double));
    double *ksi_y = (double *)malloc(num_points * sizeof(double));
    double *etta_x = (double *)malloc(num_points * sizeof(double));
    double *etta_y = (double *)malloc(num_points * sizeof(double));
    double *s2 = (double *)malloc(MAX(ni,nj) * sizeof(double));
    double *rspec = (double *)malloc(MAX(ni,nj) * sizeof(double));
    double *qv = (double *)malloc(MAX(ni,nj) * sizeof(double));
    double *dd = (double *)malloc(MAX(ni,nj) * sizeof(double));

    matrics(imax, jmax, ni, nj, x, y, jacobian, ksi_x, ksi_y, etta_x, etta_y );
    freestream (rho_0, p_0, M_0, alpha, Q, imax, jmax, ni, nj);
    BC(x, y, tel, teu, Q, ni, nj, ksi_x, ksi_y, etta_x, etta_y, jmax);

    //RHS(ni, nj, imax, jmax,  jacobian, ksi_x, ksi_y, etta_x, etta_y, Q, W, S);
    double epse = 0.006;
    double dt = 0.01;
    //smooth(Q, S, jacobian, ksi_x, ksi_y, etta_x, etta_y, ni, nj, s2, rspec, qv, dd, epse, M_0, dt);
    //print_output(ni, nj, Q, jacobian);
    print_Q(ni,nj,Q);



    // Free the dynamically allocated memory
    free(x);
    free(y);
    free(ksi_y);
    free(ksi_x);
    free(etta_x);
    free(etta_y);
    free(jacobian);
    free(p);
    free(rho);
    free(Q);

    printf("Memory freed. Program finished successfully.\n");
    return 0; 
}

