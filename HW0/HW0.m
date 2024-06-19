clc,clear, close all;
%% Input parameters
Boundary_condition = 'Dirichlet';         % can be either Dirichlet or Neumann
a = 0;                                  % starting x value for the calculation 
b = 1;                                  % ending x value for the calculation 
Boundary_condition_initial = 0;         % for Dirichlet it's Y_0 and for  Neumann its Y'_0
Boundary_condition_final   = 1;        % for Dirichlet it's Y_N and for  Neumann its Y'_N
N_values = [5 10];             % enter mesh sizes
precision_types = {'single', 'double'}; % enter precision types
% compile the C programs
mex cfd1.c;
mex cfd1_double_precision.c;
time = zeros(length(N_values), length(precision_types));
for i = 1:length(N_values)
    N = N_values(i);
    for j = 1:length(precision_types)
        precision = precision_types{j};
        %% Writing the input.txt file
        if strcmp(Boundary_condition, 'Neumann')
            bc = 'N';
        elseif strcmp(Boundary_condition, 'Dirichlet')
            bc = 'D';
        end
        ID = 'C:\Users\roiba\Documents\CFD_086376\HW0\input.txt';
        file = fopen(ID, 'wt');
        fprintf(file, '%d %f %f %c %f %f\n', N, a, b, bc, Boundary_condition_initial, Boundary_condition_final);
        fclose(file);
        x = linspace(a,b,N+1);
        %% Running the C pro    gram
        % For new users please run mex -setup in the command window and follow the
        % instructions in the Readme file on Github: https://github.com/roibaruch24/CFD_086376 
        if strcmp(precision, 'single')
            tic;
            [u, A_dig, B_dig, C_dig, d] = cfd1();
            time(i,j) = toc;
        elseif strcmp(precision, 'double')
            tic;
            [u, A_dig, B_dig, C_dig, d] = cfd1_double_precision();
            time(i,j) = toc;
        end
        %% Ploting the function
        figure(1)
        plot(x', u', 'DisplayName', ['N = ' num2str(N) ', Precision: ' char(precision)]);
        hold on;
        legend ('show','Interpreter', 'latex')
        title("y(x)")
        subtitle(['Boundary condition are: ' Boundary_condition])
        xlabel("x")
        ylabel("y")
        grid on
        %% Error calculation 
%         N_vec=linspace(1,N,N+1);
%         h=(b-a)/N;
%         switch Boundary_condition
%             case 'Dirichlet'
%                 low_dig = diag(A_dig(3:N-1),-1);
%                 mid_dig = diag(B_dig(2:N-1));
%                 high_dig = diag(C_dig(2:N-2),1);
%                 LHS  = low_dig + mid_dig + high_dig;
%                 figure(2);
%                 error =  LHS*u(2:N-1)'-d(2:N-1)';
%                 semilogy(x, abs(error(1:end-1)), 'DisplayName', ['N = ' num2str(N) ', Precision: ' char(precision)], 'LineWidth', 1);
%                 legend show
%                 hold on;
%                 title("Error on a logarithmic scale")
%                 subtitle(['Boundary condition are: ' Boundary_condition])
%                 xlabel("X")
%                 ylabel("Error")
%                 grid on;
%             case 'Neumann'
%                 low_dig = diag(A_dig(2:N+1),-1);
%                 mid_dig = diag(B_dig(1:N+1));
%                 high_dig = diag(C_dig(1:N),1);
%                 LHS  = low_dig + mid_dig + high_dig;
%                 figure(3);
%                 error = LHS*u(1:N+1)'-d(1:N+1)';
%                 semilogy(x, abs(error(1:end)), 'DisplayName', ['N = ' num2str(N) ', Precision: ' char(precision)], 'LineWidth', 1);
%                 legend show
%                 hold on;
%                 title("Error on a logarithmic scale")
%                 subtitle(['Boundary condition are: ' Boundary_condition])
%                 xlabel("X")
%                 ylabel("Error")
%                 grid on;
%         end
    end
end
% %% boundary conditions 
% % switch Boundary_condition
% %     case 'Neumann'
% %         figure(1);
% %         bc_start = Boundary_condition_initial*x(1:N/3)+(u(1)-Boundary_condition_initial*x(1));
% %         bc_end = Boundary_condition_final*x(end-N/3:end)+(u(N)-Boundary_condition_final*x(N));
% %         plot(x(1:N/3),bc_start,'--','DisplayName','$y''(0) = 1$')
% %         plot(x(end-N/3:end),bc_end,'--','DisplayName','$y''(1) = -1$')
% %         legend ('show','Interpreter', 'latex')
% % end
% %% Time plots
% figure(4);
% hold on;
% grid on;
% plot(N_values,time(:,1))
% plot(N_values,time(:,2))
% legend('Single precision','Double precision');
% xlabel('N value');
% ylabel('Run time (s)');
% title('C program run time for different mesh sizes');
% 
% figure(5);
% hold on;
% grid on;
% plot(N_values,time(:,2)-time(:,1))
% yline(0,'k--')
% xlabel('N value');
% ylabel('$\Delta T$', 'Interpreter', 'latex');
% title('$\Delta T = T_{double} - T_{single}$ for different mesh sizes', 'Interpreter', 'latex');
