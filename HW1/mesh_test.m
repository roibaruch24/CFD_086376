close,clc;

x = zeros(51,26);
y = zeros(51,26);

dx = 1/(26-12);
dy = 0.02;
xsf = 1.15;
ysf = 1.15;


y(51,2) = 0.02;

for j = 3:26
    y(51,j) = y(51,j-1) + (y(51,j-1) - y(51,j-2))*ysf;
end
t = 0.12;
x_int = 1.008930411365;
for i = 12:26
    x(i,1) = 1-cos(0.5*pi*(26-i)*dx);
    y(i,1) = 5 * t * (0.2969 * sqrt(x_int * x(i,1)) - ...
                          0.1260 * (x_int * x(i,1)) - ...
                          0.3516 * (x_int * x(i,1))^2 + ...
                          0.2843 * (x_int * x(i,1))^3 - ...
                          0.1015 * (x_int * x(i,1))^4);
end
x(26:40,1) = flipud(x(12:26,1));
for i = 26:40
    y(i,1) = -5 * t * (0.2969 * sqrt(x_int * x(i,1)) - ...
                           0.1260 * (x_int * x(i,1)) - ...
                           0.3516 * (x_int * x(i,1))^2 + ...
                           0.2843 * (x_int * x(i,1))^3 - ...
                           0.1015 * (x_int * x(i,1))^4);
end


for i = 41:51
    x(i,1) = x(i-1,1)+(x(i-1,1)-x(i-2,1))*xsf;
end
for i = 1:11
    x(i,1) = x(51-i,1);
end



x(1,2:end)  = x(1,2:end)  + x(end,1);
x(51,2:end) = x(51,2:end) + x(51,1);
% x(2:end-1,end) = x(2:end-1,end) + x(end,end);
y(1,:)  = -y(51,:);

figure()
hold on;
% theta = linspace( pi/2, 3*pi/2, 100);
% R_max = y_i_max(end);
% x_j_max_circle = R_max*cos(theta);
% y_j_max_circle = R_max*sin(theta);
% 
% R_max_vec = zeros(26,1)+R_max;
% plot(x_j_min(26:51),R_max_vec,'k-','DisplayName', 'j max line');
% plot(x_j_min(26:51),-R_max_vec,'k-','DisplayName', 'j max line');
% plot(x_j_max_circle,y_j_max_circle,'k-','DisplayName', 'j max circle');
scatter(x(:,1), y(:,1), 'k', 'DisplayName', 'j min');
scatter(x(end,:), y(end,:), 'm', 'DisplayName', 'i max');
scatter(x(1,:), y(1,:), 'm', 'DisplayName', 'i min');
legend;
xlabel('Chordwise Position (x)');
ylabel('Surface Position (y)');
title('Airfoil Surface and Points');
grid on;
hold off;

