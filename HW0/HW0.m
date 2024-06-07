clc, clear;
%% Input parameters
N = 20;
Boundary_condition = 'Dirichlet'; % can be either Dirichlet or Neumann
if strcmp(Boundary_condition, 'Neumann')
    bc = 'N';
elseif strcmp(Boundary_condition, 'Dirichlet')
    bc = 'D';
end
 a = 0; % starting x value for the calculation 
 b = 1; % ending x value for the calculation 
Boundary_condition_initial = 0; % for Dirichlet it's Y_0 and for  Neumann its Y'_0
Boundary_condition_final = 1; % for Dirichlet it's Y_N and for  Neumann its Y'_N
ID = 'C:\Users\roiba\Documents\CFD_086376\HW0\input.txt';
file = fopen(ID, 'wt');
fprintf(file, '%d %d %d %c %d %d\n', N, a, b, bc, Boundary_condition_initial, Boundary_condition_final);
fclose(file);
x = linspace(a,b,N+1);
%% Running the C program
% For new users please run mex -setup in the command window and follow the
% instructions in the Readme file on Github: https://github.com/roibaruch24/CFD_086376 
mex cfd1.c
[u, A_dig, B_dig, C_dig, d] = cfd1();
%% Ploting the function
figure
plot (x',u');
title("y(x) for N = "+num2str(N))
xlabel("x")
ylabel("y")
grid on
%% error calculation 
switch Boundary_condition
    case 'Dirichlet'
        low_dig = diag(A_dig(3:N-1),-1);
        mid_dig = diag(B_dig(2:N-1));
        high_dig = diag(C_dig(2:N-2),1);
        LHS  = low_dig + mid_dig + high_dig;
        figure;
        plot(LHS*u(2:N-1)'-d(2:N-1)');
        title("Error for N = "+num2str(N))
        xlabel("N")
        ylabel("Error")
        grid on;
    case 'Neumann'
        low_dig = diag(A_dig(2:N),-1);
        mid_dig = diag(B_dig(1:N));
        high_dig = diag(C_dig(1:N-1),1);
        LHS  = low_dig + mid_dig + high_dig;figure
        plot(LHS*u(1:N)'-d(1:N)')
        title("Error for N = "+num2str(N))
        xlabel("N")
        ylabel("Error")
        grid on;
end






