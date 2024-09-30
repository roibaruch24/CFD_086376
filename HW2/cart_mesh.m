clear, clc;
%% Input
Mach_0 = 0.9;
Alpha_0  = 3; % [deg]
P_0 = 101325; % [Pa]
Rho_0 = 1.225; % [kg / m^3]
dt = 1e-5;
res = 1e-3;
ID = 'C:\Users\roiba\Documents\CFD_086376\HW2\parametrs.txt';
        file = fopen(ID, 'wt');
        fprintf(file, '%f %f %f %f %f %f\n', Mach_0, Alpha_0, P_0, Rho_0, dt, res);
        fclose(file);
gamma = 1.4;
ni = 51;
nj = 26;
%% Run
system('HW2.exe')
%% Read Q mat
Q = table2array(readtable("Q.txt"));
Q0 = Q(1:ni*nj);
Q1 = Q(ni*nj+1:2*ni*nj);
Q2 = Q(2*ni*nj+1:3*ni*nj);
Q3 = Q(3*ni*nj+1:4*ni*nj);
Q0_mat = zeros(51,26);
Q1_mat = zeros(51,26);
Q2_mat = zeros(51,26);
Q3_mat = zeros(51,26);
k = 1;
for i=1:1:nj
    Q0_mat(:,i)=Q0(k:k+ni-1,1);
    Q1_mat(:,i)=Q1(k:k+ni-1,1);
    Q2_mat(:,i)=Q2(k:k+ni-1,1);
    Q3_mat(:,i)=Q3(k:k+ni-1,1);
    k = k+ni;
end
u = Q1_mat./Q0_mat;
v = Q2_mat./Q0_mat;
u_mag = sqrt(u.^2+v.^2);
p = (gamma-1)*(Q3_mat-0.5*Q0_mat.*u_mag.^2);
a = sqrt(gamma*p./Q0_mat);
mach = u_mag./a;
load Grid_file_mat.txt
x = Grid_file_mat(:,1);
y = Grid_file_mat(:,2);
x_mat = zeros(51,26);
y_mat = zeros(51,26);
k = 1;
for i=1:1:26
    x_mat(:,i)=x(k:k+50,1);
    y_mat(:,i)=y(k:k+50,1);
    k = k+51;
end
%% Plotting
figure
hold on
grid on
plot(x_mat(:,1), y_mat(:,1), 'k',LineWidth=2); %j min
%quiver(x_mat,y_mat,u,v,0.05)
xlim([-0.5 1.5])
ylim([-1 1])
colormap("turbo");
contourf(x_mat,y_mat,mach,30,LineStyle="none")
colorbar
cb = colorbar; % Create the colorbar and get its handle
cb.Label.String = 'Mach'; % Set the label for the colorbar
title(['NACA 0012 Mach = ' num2str(Mach_0) ', \alpha = ' num2str(Alpha_0)])


