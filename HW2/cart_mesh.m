clear, clc;
gamma = 1.4;
ni = 51;
nj = 26;
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
%%
figure
hold on
grid on
plot(x_mat(:,1), y_mat(:,1), 'k',LineWidth=2); %j min
%%quiver(x_mat,y_mat,u,v,0.05)
%%xlim([-0.2 1.2])
%%ylim([-0.2 0.2])
colormap("turbo");
contourf(x_mat,y_mat,mach,30,LineStyle="none")
colorbar

