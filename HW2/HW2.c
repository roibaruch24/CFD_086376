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
    const char *parameters_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\parametrs.txt";
    //const char *parameters_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\parametrs.txt";
    const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\input_mesh.txt";
    //const char *input_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\Grid_File2.txt";

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
    const char *input_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\input_mesh.txt";
    //const char *input_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\Grid_File2.txt";

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
    double a = sqrt(gamma*p/rho);
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
            double u_1 = Q[offset3d(i, 1, 1,ni, nj)]/Q[offset3d(i, 1, 0, ni, nj)];
            double v_1 = Q[offset3d(i, 1, 2,ni, nj)]/Q[offset3d(i, 1, 0, ni, nj)];
            double u = ((ksi_x[offset2d(i, 1, ni)]*u_1 )+ksi_y[offset2d(i, 1, ni)]*v_1)/(ksi_x[offset2d(i, 0, ni)]-ksi_y[offset2d(i, 0, ni)]*(etta_x[offset2d(i, 0, ni)]/etta_y[offset2d(i, 0, ni)]));
            double v =(-etta_x[offset2d(i,0,ni)])/(etta_y[offset2d(i,0,ni)])*u;
            if (etta_y[offset2d(i,0,ni)]==0)
            {
            v = (ksi_x[offset2d(i,1,ni)]*u_1+(ksi_y[offset2d(i,1,ni)])*v_1-ksi_x[offset2d(i,0,ni)]*u)/ksi_y[offset2d(i,0,ni)];
            }       
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
void RHS(int ni, int nj, int imax, int jmax, double *Jacobian, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y, double *Q, double *W, double *S, double dt ){   
    // Zero S
    for (int i = 0;i <= imax;i++){
        for(int j = 0;j <= jmax;j++){
            for (int k = 0; k <=3;k++){
                S[offset3d(i,j,k,ni,nj)] = 0;
            }
        }
    }
    /* ksi direction */ 
    for (int j = 1; j < jmax; j++){
        for ( int i = 0; i <= imax; i++){
            double u = Q[offset3d(i, j, 1, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double v = Q[offset3d(i, j, 2, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double U_conv = ksi_x[offset2d(i, j, ni)]*u + ksi_y[offset2d(i, j, ni)]*v;
            double P = (gamma-1)*(Q[offset3d(i,j,3,ni,nj)]-0.5*Q[offset3d(i,j,0,ni,nj)]*(pow(u,2)+pow(v,2)));

            W[offset2d(i, 0, ni)] = Q[offset3d(i, j, 0, ni, nj)]*U_conv/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 1, ni)] = (Q[offset3d(i, j, 1, ni, nj)]*U_conv + ksi_x[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 2, ni)] = (Q[offset3d(i, j, 2, ni, nj)]*U_conv + ksi_y[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(i, 3, ni)] = (Q[offset3d(i,j,3,ni,nj)]+P)*U_conv/Jacobian[offset2d(i,j,ni)];
    }
    for(int i = 1; i< imax; i++){
        for( int k = 0; k < 4; k++ ){
            S[offset3d(i, j, k, ni, nj)] += -0.5*(W[offset2d(i+1,k, ni)]-W[offset2d(i-1,k,ni)]);
        }
    }

}

/* etta direction */ 
    for (int i = 1; i < imax; i++){
        for ( int j = 0; j <= jmax; j++){
            double u = Q[offset3d(i, j, 1, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double v = Q[offset3d(i, j, 2, ni, nj)]/Q[offset3d(i, j, 0, ni, nj)];
            double V_conv = etta_x[offset2d(i, j, ni)]*u + etta_y[offset2d(i, j, ni)]*v;
            double P = (gamma-1)*(Q[offset3d(i,j,3,ni,nj)]-0.5*Q[offset3d(i,j,0,ni,nj)]*(pow(u,2)+pow(v,2)));
            W[offset2d(j, 0, nj)] = (Q[offset3d(i, j, 0, ni, nj)]*V_conv)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 1, nj)] = (Q[offset3d(i, j, 1, ni, nj)]*V_conv + etta_x[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 2, nj)] = (Q[offset3d(i, j, 2, ni, nj)]*V_conv + etta_y[offset2d(i, j, ni)]*P)/Jacobian[offset2d(i, j, ni)];
            W[offset2d(j, 3, nj)] = (Q[offset3d(i,j,3,ni,nj)]+P)*V_conv/Jacobian[offset2d(i,j,ni)];
        }
        for(int j = 1; j < jmax ; j++){
            for( int k = 0; k < 4; k++ ){
                S[offset3d(i, j, k, ni, nj)] += -0.5*(W[offset2d(j+1,k, nj)]-W[offset2d(j-1,k,nj)]); 
            }
        }

    }
    for(int i = 0; i<= imax; i++){
        for(int j = 0; j <= jmax; j++){
            for(int k = 0; k <= 3; k++ ){
                S[offset3d(i, j, k, ni, nj)] *=dt;
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
void Calc_Jacob_LHS(int i, int ni, double k_x, double k_y, double Q_0, double Q_1, double Q_2, double Q_3, double *B){
    
    double u = Q_1/Q_0;
    double v = Q_2/Q_0;
    double psi_square = 0.5*(gamma -1)*(pow(u,2) + pow(v,2));
    double theta = k_x*u + k_y*v;
    double gamma_1 = gamma - 1;
    double gamma_2 = gamma - 2;
    double beta = gamma*Q_3/Q_0 - psi_square;

    B[offset3d(i, 0, 0, ni, 4)] = 0;
    B[offset3d(i, 0, 1, ni, 4)] = k_x;
    B[offset3d(i, 0, 2, ni, 4)] = k_y;
    B[offset3d(i, 0, 3, ni, 4)] = 0;

    B[offset3d(i, 1, 0, ni, 4)] = k_x*psi_square - u*theta;
    B[offset3d(i, 1, 1, ni, 4)] = theta -k_x*gamma_2*u;
    B[offset3d(i, 1, 2, ni, 4)] = k_y*u - gamma_1*k_x*v;
    B[offset3d(i, 1, 3, ni, 4)] = k_x*gamma_1;

    B[offset3d(i, 2, 0, ni, 4)] = k_y*psi_square - v*theta;
    B[offset3d(i, 2, 1, ni, 4)] = k_x*v - k_y*gamma_1*u;
    B[offset3d(i, 2, 2, ni, 4)] = theta - k_y*gamma_2*v;
    B[offset3d(i, 2, 3, ni, 4)] = k_y*gamma_1;

    B[offset3d(i, 3, 0, ni, 4)] = theta*(2*psi_square - gamma*Q_3/Q_0);
    B[offset3d(i, 3, 1, ni, 4)] = k_x*beta - gamma_1*u*theta;
    B[offset3d(i, 3, 2, ni, 4)] = k_y*beta - gamma_1*v*theta;
    B[offset3d(i, 3, 3, ni, 4)] = gamma*theta;


}
void LHSX(int ni, int nj, int j, double dt, double *Q,double *A, double *B, double *C, double *ksi_x, double *ksi_y){
    for (int i = 0; i < ni; i++){
        double Q_0 = Q[offset3d(i,j,0,ni,nj)];
        double Q_1 = Q[offset3d(i,j,1,ni,nj)];
        double Q_2 = Q[offset3d(i,j,2,ni,nj)];
        double Q_3 = Q[offset3d(i,j,3,ni,nj)];
        Calc_Jacob_LHS(i, ni, ksi_x[offset2d(i, j, ni)], ksi_y[offset2d(i, j, ni)], Q_0, Q_1, Q_2, Q_3,B);
    }
    for (int n = 0; n < 4; n++){
        for (int m = 0; m < 4; m++){
            for (int i = 1; i < ni - 1; i++){
                A[offset3d(i, m, n, ni, 4)] = -0.5*B[offset3d(i -1, m , n, ni, 4)]*dt;
                C[offset3d(i, m, n, ni, 4)] = 0.5*B[offset3d(i + 1, m , n, ni, 4)]*dt;
            }
        }
    }
for (int i = 0; i < ni-1; i++){
    for (int m = 0; m < 4; m++){
        for (int n = 0; n < 4; n++){
            if (m!=n){
                B[offset3d(i, m, n, ni, 4)] = 0;
            }
            if (m == n){
                B[offset3d(i, m, n, ni, 4)] = 1;
            }
        }
    }
}

}
void LHSY(int ni, int nj, int i, double dt, double *Q,double *A, double *B, double *C, double *etta_x, double *etta_y){
    for (int j = 0; j < nj; j++){
        double Q_0 = Q[offset3d(i,j,0,ni,nj)];
        double Q_1 = Q[offset3d(i,j,1,ni,nj)];
        double Q_2 = Q[offset3d(i,j,2,ni,nj)];
        double Q_3 = Q[offset3d(i,j,3,ni,nj)];
        Calc_Jacob_LHS(j, nj, etta_x[offset2d(i, j, ni)], etta_y[offset2d(i, j, ni)], Q_0, Q_1, Q_2, Q_3,B);
    }
    for (int n = 0; n < 4; n++){
        for (int m = 0; m < 4; m++){
            for (int j = 1; j < nj - 1; j++){
                A[offset3d(j, m, n, nj, 4)] = -0.5*B[offset3d(j -1, m , n, nj, 4)]*dt;
                C[offset3d(j, m, n, nj, 4)] = 0.5*B[offset3d(j + 1, m , n, nj, 4)]*dt;
            }
        }
    }
for (int j = 0; j < nj-1; j++){
    for (int m = 0; m < 4; m++){
        for (int n = 0; n < 4; n++){
            if (m!=n){
                B[offset3d(j, m, n, nj, 4)] = 0;
            }
            if (m == n){
                B[offset3d(j, m, n, nj, 4)] = 1;
            }
        }
    }
}
}

int smoothx(double *q, double *xx, double *xy, int id, int jd, double *a,
	   double *b, double *c, int j,double *jac, double *drr, double *drp, 
           double *rspec, double *qv, double *dd,
           double epsi, double fsmach, double dt)
{

    double *rho, *u_vel, *v_vel, *t_e;

    double eratio, gm1, ggm1, eps, ra, u, v, qq, ss,
          qav, qxx, rr, rp;
    int ib, ie, i, offset, offsetp1, offsetm1, ip, ir, n;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    ib = 1;
    ie = id - 1;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*     smoothing in xi direction */

        for (i = 0; i < id; i++) {
	    offset = id * j + i;
	    eps = epsi / jac[offset];
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

        for (i = ib; i < ie; i++) {
	    ip = i + 1;
	    ir = i - 1;
	    offset = j * id + i;
	    offsetp1 = j * id + i + 1;
	    offsetm1 = j * id + i - 1;
	    rp = (0.5 * (dd[ip] + dd[i]) + 2.5) * dt * 0.5 * (rspec[ip] + rspec[i]);
	    rr = (0.5 * (dd[ir] + dd[i]) + 2.5) * dt * 0.5 * (rspec[ir] + rspec[i]);
	    qv[i] = (rr + rp) * jac[offset];
	    drr[i] = rr * jac[offsetm1];
	    drp[i] = rp * jac[offsetp1];
	}

	for (n = 0; n < 4; n++) {
	    for (i = ib; i < ie; i++) {
	        offset = (n * 4 + n) * id + i;
		a[offset] -= drr[i];
		b[offset] += qv[i];
		c[offset] -= drp[i];
	    }
        }
    return 0;
} 

int smoothy(double *q, double *yx, double *yy, int id, int jd, double *a,
	   double *b, double *c, int i,double *jac, double *drr, double *drp, 
           double *rspec, double *qv, double *dd,
           double epsi, double fsmach, double dt)
{

    double *rho, *u_vel, *v_vel, *t_e;

    double eratio, smool, gm1, ggm1, eps, ra, u, v, qq, ss, 
          qav, ssfs, qyy, rp, rr;
    int jb, je, j, offset, offsetp1, offsetm1, n,
        jp, jr;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    smool = 1.0;
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    jb = 1;
    je = jd - 1;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*     smoothing in eta direction */

        ssfs = 1. / (0.001 + fsmach * fsmach);
        for (j = 0; j < jd; j++) {
	    offset = id * j + i;
	    eps = epsi / jac[offset];
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

        for (j = jb; j < je; j++) {
	    jp = j + 1;
	    jr = j - 1;
	    offset = j * id + i;
	    offsetp1 = (j + 1) * id + i;
	    offsetm1 = (j - 1) * id + i;
	    rp = (0.5 * (dd[jp] + dd[j]) + 2.5) * dt * 0.5 * (rspec[jp] + rspec[j]);
	    rr = (0.5 * (dd[jr] + dd[j]) + 2.5) * dt * 0.5 * (rspec[jr] + rspec[j]);
	    qv[j] = (rr + rp) * jac[offset];
	    drr[j] = rr * jac[offsetm1];
	    drp[j] = rp * jac[offsetp1];
	}

	for (n = 0; n < 4; n++) {
	    for (j = jb; j < je; j++) {
	        offset = (n * 4 + n) * jd + j;
		a[offset] -= drr[j];
		b[offset] += qv[j];
		c[offset] -= drp[j];
	    }
        }
    return 0;
} /* smoothy */
int btri4s(double *a, double *b, double *c, double *f, int kd, int ks, int ke)
{
  /* Local variables */
  int k, m, n, nd, md;

  double c1, d1, d2, d3, d4, c2, c3, c4, b11, b21, b22, b31, b32, b33, 
    b41, b42, b43, b44, u12, u13, u14, u23, u24, u34;
  
  
  /*   (A,B,C)F = F, F and B are overloaded, solution in F */

  md = 4;
  nd = 4;

  /*   Part 1. Forward block sweep */
  
  for (k = ks; k <= ke; k++)
    {
      
      /*      Step 1. Construct L in B */
      
      if (k != ks) 
	{
	  for (m = 0; m < md; m++) 
	    {
	      for (n = 0; n < nd; n++) 
		{
		  b[k + kd * (m + md * n)] = b[k + kd * (m + md * n)] 
		    - a[k + kd * (m + md * 0)] * b[k - 1 + kd * (0 + md * n)] 
		    - a[k + kd * (m + md * 1)] * b[k - 1 + kd * (1 + md * n)] 
		    - a[k + kd * (m + md * 2)] * b[k - 1 + kd * (2 + md * n)] 
		    - a[k + kd * (m + md * 3)] * b[k - 1 + kd * (3 + md * n)] ;
		}
	    }
	}
      
      /*      Step 2. Compute L inverse (block matrix) */
      
      /*          A. Decompose L into L and U */
      
      b11 = 1. / b[k + kd * (0 + md * 0)];
      u12 = b[k + kd * (0 + md * 1)] * b11;
      u13 = b[k + kd * (0 + md * 2)] * b11;
      u14 = b[k + kd * (0 + md * 3)] * b11;
      b21 = b[k + kd * (1 + md * 0)];
      b22 = 1. / (b[k + kd * (1 + md * 1)] - b21 * u12);
      u23 = (b[k + kd * (1 + md * 2)] - b21 * u13) * b22;
      u24 = (b[k + kd * (1 + md * 3)] - b21 * u14) * b22;
      b31 = b[k + kd * (2 + md * 0)];
      b32 = b[k + kd * (2 + md * 1)] - b31 * u12;
      b33 = 1. / (b[k + kd * (2 + md * 2)] - b31 * u13 - b32 * u23);
      u34 = (b[k + kd * (2 + md * 3)] - b31 * u14 - b32 * u24) * b33;
      b41 = b[k + kd * (3 + md * 0)];
      b42 = b[k + kd * (3 + md * 1)] - b41 * u12;
      b43 = b[k + kd * (3 + md * 2)] - b41 * u13 - b42 * u23;
      b44 = 1. / (b[k + kd * (3 + md * 3)] - b41 * u14 - b42 * u24 
		  - b43 * u34);
      
      /*      Step 3. Solve for intermediate vector */
      
      /*          A. Construct RHS */
      if (k != ks) 
	{
	  for (m = 0; m < md; m++) 
	    {
	      f[k + kd * m] = f[k + kd * m] 
		- a[k + kd * (m + md * 0)] * f[k - 1 + kd * 0] 
		- a[k + kd * (m + md * 1)] * f[k - 1 + kd * 1] 
		- a[k + kd * (m + md * 2)] * f[k - 1 + kd * 2] 
		- a[k + kd * (m + md * 3)] * f[k - 1 + kd * 3];
	    }
	}
      
      /*          B. Intermediate vector */
      
      /*          Forward substitution */
      
      d1 = f[k + kd * 0] * b11;
      d2 = (f[k + kd * 1] - b21 * d1) * b22;
      d3 = (f[k + kd * 2] - b31 * d1 - b32 * d2) * b33;
      d4 = (f[k + kd * 3] - b41 * d1 - b42 * d2 - b43 * d3) * b44;
      
      /*          Backward substitution */
      
      f[k + kd * 3] = d4;
      f[k + kd * 2] = d3 - u34 * d4;
      f[k + kd * 1] = d2 - u23 * f[k + kd * 2] - u24 * d4;
      f[k + kd * 0] = d1 - u12 * f[k + kd * 1] - u13 * f[k + kd * 2] - u14 * d4;
      
      /*      Step 4. Construct U = L ** (-1) * C */
      /*              by columns and store in B */
      
      if (k != ke) 
	{
	  for (n = 0; n < nd; n++) 
	    {

	      /*          Forward substitution */
	      
	      c1 = c[k + kd * (0 + md * n)] * b11;
	      c2 = (c[k + kd * (1 + md * n)] - b21 * c1) * b22;
	      c3 = (c[k + kd * (2 + md * n)] - b31 * c1 - b32 * c2) * 
		b33;
	      c4 = (c[k + kd * (3 + md * n)] - b41 * c1 - b42 * c2 - 
		    b43 * c3) * b44;
	      
	      /*          Backward substitution */
	      
	      b[k + kd * (3 + md * n)] = c4;
	      b[k + kd * (2 + md * n)] = c3 - u34 * c4;
	      b[k + kd * (1 + md * n)] = c2 - u23 * b[k + kd * (2 + md * n)] - u24 * c4;
	      b[k + kd * (0 + md * n)] = c1 - u12 * b[k + kd * (1 + md * n)] 
		- u13 * b[k + kd * (2 + md * n)] - u14 * c4;
	    }
	}
    }
  
  /*   Part 2. Backward block sweep */
  
  if (ke == ks) 
    {
      return 0;
    }

  for (k = ke - 1; k >= ks; --k) 
    {
      for (m = 0; m < md; m++) 
	{
	  f[k + kd * m] = f[k + kd * m] 
	    - b[k + kd * (m + md * 0)] * f[k + 1 + kd * 0] 
	    - b[k + kd * (m + md * 1)] * f[k + 1 + kd * 1] 
	    - b[k + kd * (m + md * 2)] * f[k + 1 + kd * 2] 
	    - b[k + kd * (m + md * 3)] * f[k + 1 + kd * 3];
	}
    }
  
  return 0;
  
} /* btri4s_ */
void step(int ni, int nj, double dt, int imax, int jmax, double *jacobian, double *Q, double *W, double *S,double *s2, double *rspec, double *qv, double *dd, double epse, double M_0, double *A, double *B, double *C, double *D, double *ksi_x, double *ksi_y, double *etta_x, double *etta_y, double *drr, double *drp){
    // Call RHS
    RHS(ni, nj, imax, jmax,  jacobian, ksi_x, ksi_y, etta_x, etta_y, Q, W, S, dt);
    smooth(Q, S, jacobian, ksi_x, ksi_y, etta_x, etta_y, ni, nj, s2, rspec, qv, dd, epse, M_0, dt);
    // ksi inversions
    for(int j = 1; j< nj -1; j++){
       LHSX(ni, nj, j, dt, Q, A, B, C, ksi_x, ksi_y);
       smoothx(Q, ksi_x, ksi_y, ni, nj,  A, B, C, j, jacobian, drr, drp, rspec, qv, dd, 2*epse, M_0, dt);
       for (int k = 0; k < 4; k++){
            for (int i = 0; i < ni -1; i++){
                D[offset2d(i, k, ni)] = S[offset3d(i, j, k, ni, 4)];
            }
       }
       btri4s(A, B, C, D, ni, 1, 1);
       for (int k = 0; k < 4; k++){
            for (int i = 0; i < ni -1; i++){
                S[offset3d(i, j, k, ni, 4)] = D[offset2d(i, k, ni)];
            }
       }
    }

    // etta inversions
    for(int i = 1; i< ni -1 ; i++){
       LHSY(ni, nj, i, dt, Q, A, B, C, etta_x, etta_y); 
       smoothy(Q, etta_x, etta_y, ni, nj,  A, B, C, i, jacobian, drr, drp, rspec, qv, dd, 2*epse, M_0, dt);
       for (int k = 0; k < 4; k++){
            for (int j = 0; j < nj -1; j++){
                D[offset2d(j, k, nj)] = S[offset3d(i, j, k, ni, 4)]; // Load D from S
            }
       }

       btri4s(A, B, C, D, nj, 1, 1);

       for (int k = 0; k < 4; k++){
            for (int j = 0; j < nj -1; j++){
                S[offset3d(i, j, k, ni, 4)] = D[offset2d(j, k, nj)];
            }
       }
    }
    // Upating the solution
    for (int i = 1; i < imax; i++){
            for (int j = 1; j < jmax; j++){
                for(int k = 0; k <= 3; k++){
                    Q[offset3d(i, j, k, ni, nj)] += S[offset3d(i, j, k, ni, nj)]*jacobian[offset2d(i, j, ni)];
                }
            }
        }
}
int convergence(int i, int j, int ni, int nj, double Q){


}

int print_output(int ni,int nj,double *x, double *y)
{
    char *output_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\output.txt";
    //char *output_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\output.txt";
    
    FILE *output= fopen(output_file_path, "wt"); //crating output file
    for (int i = 0; i < ni*nj; i++){
            fprintf(output, "%f %f\n", x[i], y[i]);
        }
    
    fclose(output);
    return 0;
}
int print_Q(int ni,int nj,double *Q)
{
    char *output_file_path = "C:\\Users\\roiba\\Documents\\CFD_086376\\HW2\\Q.txt";
    //char *output_file_path = "C:\\Users\\roiB\\Desktop\\CFD\\CFD_086376\\HW2\\Q.txt";
    
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
    double epse = 0.06;
    double dt = pow(10,-6);
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
    double *A = (double *)malloc(MAX(ni, nj) * 16 * sizeof(double));
    double *B = (double *)malloc(MAX(ni, nj) * 16 * sizeof(double));
    double *C = (double *)malloc(MAX(ni, nj) * 16 * sizeof(double));
    double *D = (double *)malloc(MAX(ni, nj) * 4 * sizeof(double));
    double *ksi_x = (double *)malloc(num_points * sizeof(double));
    double *ksi_y = (double *)malloc(num_points * sizeof(double));
    double *etta_x = (double *)malloc(num_points * sizeof(double));
    double *etta_y = (double *)malloc(num_points * sizeof(double));
    double *s2 = (double *)malloc(MAX(ni,nj) * sizeof(double));
    double *rspec = (double *)malloc(MAX(ni,nj) * sizeof(double));
    double *qv = (double *)malloc(MAX(ni,nj) * sizeof(double));
    double *dd = (double *)malloc(MAX(ni,nj) * sizeof(double));
    double *drr = (double *)malloc(MAX(ni,nj) * sizeof(double));
    double *drp = (double *)malloc(MAX(ni,nj) * sizeof(double));
    input_mesh(x, y, ni, nj);

    matrics(imax, jmax, ni, nj, x, y, jacobian, ksi_x, ksi_y, etta_x, etta_y );
    freestream (rho_0, p_0, M_0, alpha, Q, imax, jmax, ni, nj);
    /*for (int a = 0; a < 1000; a++)
	{
		BC(x, y, tel, teu, Q, ni, nj, ksi_x, ksi_y, etta_x, etta_y, jmax);
		RHS(ni, nj, imax, jmax,  jacobian, ksi_x, ksi_y, etta_x, etta_y, Q, W, S, dt);
		smooth(Q,S,jacobian,ksi_x,ksi_y,etta_x,etta_y,ni,nj,s2,rspec,qv,dd,epse,M_0,dt);
		for (int i = 1;i < imax;i++)
		{
			for(int j = 1;j < jmax;j++)
			{
				for (int k = 0; k <=3;k++)
				{
					Q[offset3d(i,j,k,ni,nj)] += S[offset3d(i,j,k,ni,nj)]*jacobian[offset2d(i,j,ni)];
				}
			}
		}
	}*/

   for (int iter = 0; iter < 20000; iter++){
        BC(x, y, tel, teu, Q, ni, nj, ksi_x, ksi_y, etta_x, etta_y, jmax);
        step(ni, nj, dt, imax, jmax, jacobian, Q, W, S, s2, rspec, qv, dd, epse, M_0, A, B, C, D, ksi_x, ksi_y, etta_x, etta_y, drr, drp);
    }
        
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
    free(S);
    free(W);
    free(A);
    free(B);
    free(C);
    free(D);
    free(s2);
    free(rspec);
    free(qv);
    free(dd);
    free(drr);
    free(drp);


    printf("yalla\n");
    return 0; 
}

