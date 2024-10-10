clear, clc, close all;
% ================================= HW2 ==============================
% This script is used for running the CFD solver written in HW2. The input
% for the script is:
%
% 1. main_path               --> directory of the main folder for HW2 where the .exe file is
% located.
%
% 2. mesh_path               --> directory of the folder where HW1 is located
% including the .exe file.
%
% 3. Input_file_name         --> the name of the file with the input parameters.
%  Must be located in the main_path !!! 
%
% 4. Input parametrs         --> free stream values:
% M_0     - free stream Mach number
% Alpha_0 - free stream angle of attack in degrees
% P_0     - free stream pressure in pascal
% Rho_0   - free stream density in kg/m^3
% dt      - time step for the calculation
% res     - orders of magnitude for the residue to reduce for convergences
% ni, nj  - mesh size values for the grid from HW1
%
% 5. Output_file_name        --> name of the output file with the Q matrix
%
% 6. Residue_file_name       --> name of the residue file
%
% 7. Mesh_file_name          --> name of the mesh file
%
% 8. Mesh_file_name_original --> name of the mesh file in HW1 directory 
% (can be output.txt etc...)
%
% 9. Flags                   --> Flags for the type of run wanted:
% Run_CFD_flag    - 1 for running the CFD program and 0 for just post
% processing  exiting data.
% Run_Mesh_flag   - either 1, 2, 3, 4 for different control cases or 0 for
% not creating a new mesh.
% Mach_plot_flag - 1 for plotting the Mach field
% Cp_plot_flag   - 1 for plotting the pressure coefficient field
% Mesh_plot      - 1 for plotting the Mesh
% Residue_plot   - 1 for plotting the residue
% ================================= Paths =================================

main_path = 'C:\Users\roiba\Documents\CFD_086376\HW2';
mesh_path = 'C:\Users\roiba\Documents\CFD_086376\HW1';

% ================================= Filenames =============================

Input_file_name = 'parametrs.txt';
Output_file_name = 'Q.txt';
Residue_file_name = 'residue.txt';
Mesh_file_name = 'input_mesh.txt';
Mesh_file_name_original = 'output.txt';
% ================================= Flags =================================
Flag.Run_CFD_flag   = 1;
Flag.Run_Mesh_flag  = 0;
Flag.Mach_plot_flag = 1;
Flag.Cp_plot_flag   = 1;
Flag.Mesh_plot      = 0;
Flag.Residue_plot   = 0;

% ============================= Manual Inputs =============================

Input.Mach_0 = 0.3;
Input.Alpha_0  = 0; % [deg]
Input.P_0 = 101325; % [Pa]
Input.Rho_0 = 1.225; % [kg / m^3]
Input.dt = 1e-5;
Input.res = 1e-3;
Input.Gamma_0 = 1.4;
ni = 51;
nj = 26;

% ================================= Main ==================================
if Flag.Run_Mesh_flag
    Run_Mesh (main_path, mesh_path, Mesh_file_name, Mesh_file_name_original, Flag)
end

if Flag.Run_CFD_flag 
    Input_parameters(Input,Input_file_name, main_path);
    system('HW2.exe')
end
Output = Post_process(Output_file_name, main_path, Mesh_file_name, Residue_file_name, ni, nj,Input);
if Flag.Mach_plot_flag
    plot_mach(Input, Output)
end
if Flag.Cp_plot_flag
    plot_cp(Input, Output)
end
if Flag.Mesh_plot
    plot_mesh(Output)
end
if Flag.Residue_plot
    plot_residue(Input,Output)
end

Output_forces =  Aero_forces (Output, Input)

% =============================== Functions ===============================

function Input_parameters(Input,Input_file_name, main_path)
cd(main_path)
file = fopen(Input_file_name, 'wt');
fprintf(file, '%f %f %f %f %f %f\n', Input.Mach_0, Input.Alpha_0, Input.P_0, Input.Rho_0, Input.dt, Input.res);
fclose(file);
end

function Run_Mesh (main_path, mesh_path, Mesh_file_name, Mesh_file_name_original, Flag)
cd(mesh_path)
Mesh_input(Flag.Run_Mesh_flag)
system('HW1.exe');
copyfile([mesh_path '\' Mesh_file_name_original], [main_path '\' Mesh_file_name]);
end

function plot_mesh(Output)
figure
hold on
grid minor
plot(Output.x_mat(:,1), Output.y_mat(:,1), 'k', 'DisplayName', 'j min');
plot(Output.x_mat(end,:), Output.y_mat(end,:), 'k', 'DisplayName', 'i max');
plot(Output.x_mat(1,:), Output.y_mat(1,:), 'k', 'DisplayName', 'i min');
plot(Output.x_mat(:,26), Output.y_mat(:,26), 'k', 'DisplayName', 'j max');
for j = 2:1:25
    plot(Output.x_mat(:,j), Output.y_mat(:,j),'b','LineWidth',0.5)
end
for i = 2:1:50
    plot(Output.x_mat(i,:), Output.y_mat(i,:),'b','LineWidth',0.5)
end
end

function Output = Post_process(Output_file_name, main_path, Mesh_file_name, Residue_file_name, ni, nj, Input)
cd(main_path)
Q = table2array(readtable(Output_file_name));
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
% extracting paramters from the Q matrix 
u = Q1_mat./Q0_mat;
v = Q2_mat./Q0_mat;
u_mag = sqrt(u.^2+v.^2);
Output.p = (Input.Gamma_0-1)*(Q3_mat-0.5*Q0_mat.*u_mag.^2);
a = sqrt(Input.Gamma_0*Output.p./Q0_mat);
Output.Mach = u_mag./a;
fid = fopen(Mesh_file_name, 'r');
fgetl(fid);
Grid_file_mat = textscan(fid, '%f %f');
fclose(fid);


x = cell2mat(Grid_file_mat(:,1));
y = cell2mat(Grid_file_mat(:,2));
Output.x_mat = zeros(ni,nj);
Output.y_mat = zeros(ni,nj);
k = 1;
for i=1:1:nj
    Output.x_mat(:,i)=x(k:k+ni-1,1);
    Output.y_mat(:,i)=y(k:k+ni-1,1);
    k = k+ni;
end
q_dyn_0 = 0.5*1.4*Input.P_0*(Input.Mach_0)^2;
Output.Cp = (Output.p - Input.P_0)./q_dyn_0;
Output.Residue = load (Residue_file_name);
end

function plot_mach(Input, Output)
figure
hold on
grid on
plot(Output.x_mat(:,1), Output.y_mat(:,1), 'k',LineWidth=2); %j min
xlim([-0.5 1.5])
ylim([-1 1])
colormap("turbo");
contourf(Output.x_mat,Output.y_mat,Output.Mach,390,LineStyle="none")
colorbar
cb = colorbar; 
cb.Label.String = 'Mach'; % Set the label for the colorbar
title(['Mach plot NACA 0012 Mach = ' num2str(Input.Mach_0) ', \alpha = ' num2str(Input.Alpha_0)])
end

function plot_cp(Input, Output)
figure
hold on
grid on
plot(Output.x_mat(:,1), Output.y_mat(:,1), 'k',LineWidth=2); %j min
xlim([-0.5 1.5])
ylim([-1 1])
colormap("turbo");
contourf(Output.x_mat,Output.y_mat,Output.Cp,300,LineStyle="none")
colorbar
cb = colorbar; 
cb.Label.String = 'C_p'; % Set the label for the colorbar
title(['C_p plot NACA 0012 Mach = ' num2str(Input.Mach_0) ', \alpha = ' num2str(Input.Alpha_0)])
end

function plot_residue(Input,Output)
figure
semilogy(Output.Residue(:,1),Output.Residue(:,2))
title(['Residue plot for M = ' num2str(Input.Mach_0) ', \alpha = ' num2str(Input.Alpha_0)])
xlabel('Iteration')
ylabel('Residue')
grid on
end

function Output_forces =  Aero_forces (Output, Input)
%%
% Define airfoil coordinates (x, y) and corresponding Cp values
x_airfoil = Output.x_mat(12:40,1);  % x-coordinates of surface points
y_airfoil = Output.y_mat(12:40,1);  % y-coordinates of surface points
Cp = Output.Cp(12:40,1);
% Number of points on the surface
n_points = length(x_airfoil);

% Preallocate force component arrays
Cx = zeros(n_points-1, 1);  % Force in x-direction
Cz = zeros(n_points-1, 1);  % Force in z-direction

% Loop over each surface segment to calculate forces
for i = 1:(n_points-1)
    % Coordinates of the current surface segment
    x1 = x_airfoil(i);
    y1 = y_airfoil(i);
    x2 = x_airfoil(i+1);
    y2 = y_airfoil(i+1);

    % Calculate the length of the surface segment
    dS = sqrt((x2 - x1)^2 + (y2 - y1)^2);

    % Calculate the normal vector (nx, nz)
    tangent_x = x2 - x1;
    tangent_y = y2 - y1;

    % Calculate the normal vector (perpendicular to the tangent)
    normal_x = -tangent_y;  % Outward normal in x-direction
    normal_z = tangent_x;   % Outward normal in z-direction
    magnitude = sqrt(normal_x^2 + normal_z^2);
    nx = normal_x / magnitude;  % Normal vector in x-direction
    nz = normal_z / magnitude;  % Normal vector in z-direction

    % Make sure normals point outward:
    % If the airfoil's surface points are ordered clockwise, this is OK.
    % If they are ordered counterclockwise, reverse the normal direction:
    nx = -nx;  % Uncomment this if your normals point inward
    nz = -nz;  % Uncomment this if your normals point inward

    % Average pressure coefficient for this surface segment
    Cp_avg = (Cp(i) + Cp(i+1)) / 2;


    % Calculate force components on this segment
   Cx(i) = Cp_avg * dS * nx;  % Force in x-direction
   Cz(i) = Cp_avg * dS * nz;  % Force in z-direction
end
figure
plot(x_airfoil, Cp)
Total_Fx = trapz(y_airfoil(1:end-1), Cx);  
Total_Fy = trapz(x_airfoil(1:end-1), Cz);  

Output_forces.Lift = Total_Fy * cosd(Input.Alpha_0) - Total_Fx * sind(Input.Alpha_0);
Output_forces.Drag = Total_Fx * cosd(Input.Alpha_0) + Total_Fy * sind(Input.Alpha_0);

end