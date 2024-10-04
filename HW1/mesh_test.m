
%% scattering
load output.txt
x = output(:,1);
y = output(:,2);

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
grid minor
plot(x_mat(:,1), y_mat(:,1), 'k', 'DisplayName', 'j min');
plot(x_mat(end,:), y_mat(end,:), 'k', 'DisplayName', 'i max');
plot(x_mat(1,:), y_mat(1,:), 'k', 'DisplayName', 'i min');
plot(x_mat(:,26), y_mat(:,26), 'k', 'DisplayName', 'j max');
for j = 2:1:25
    plot(x_mat(:,j), y_mat(:,j),'b','LineWidth',0.5)
end
for i = 2:1:50
    plot(x_mat(i,:), y_mat(i,:),'b','LineWidth',0.5)
end
title('Mesh for control case 2');
%%
load residue.txt
n = residue(:,3);
L_x = residue(:,1);
L_y = residue(:,2);
figure
semilogy(n,L_x)
hold on
grid on
semilogy(n,L_y)
xlabel("Iteration")
ylabel('$L_{max}$', 'Interpreter','latex')
title("Convergences")
legend('L_x','L_y')


%%
w_vec = [0.5 0.8 1 1.3 1.6 1.9];
n_r_01_vec = [6413 4005 3202 2461 1998 1681];
n_r_005_vec = [3745 2337 1868 1435 1164 979];
n_r_001_vec = [7727 4822 3852 2952 2376 2091];

figure 
loglog(w_vec, n_r_01_vec)
hold on
grid on
loglog(w_vec, n_r_005_vec)
loglog(w_vec, n_r_001_vec)
title('Iteration number for vs $\omega$ values for varing r values','Interpreter','latex')
xlabel('$\omega$', 'Interpreter','latex')
ylabel('Iterations')
legend('r = 0.01', 'r = 0.005', 'r = 0.001')

figure
loglog([0.01 0.005 0.001], [n_r_01_vec(3) n_r_005_vec(3) n_r_001_vec(3)])
hold on;
grid on;
loglog([0.01 0.005 0.001], [n_r_01_vec(4) n_r_005_vec(4) n_r_001_vec(4)])
loglog([0.01 0.005 0.001], [n_r_01_vec(5) n_r_005_vec(5) n_r_001_vec(5)])
title('Iteration number for vs r values for varing $\omega$ values','Interpreter','latex')
xlabel('r')
ylabel('Iterations')
legend('$\omega$ = 1','$\omega$ = 1.3','$\omega$ = 1.6', 'Interpreter','latex')


