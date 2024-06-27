% clear,close,clc;
% 
% x = zeros(51,26);
% y = zeros(51,26);
% 
% dx = 1/(26-12);
% dy = 0.02;
% xsf = 1.15;
% ysf = 1.15;
% i_max = 51;
% j_max = 26;
% 
% 
% %% j min
% t = 0.12;
% x_int = 1.008930411365;
% for i = 12:26
%     x(i,1) = 1-cos(0.5*pi*(26-i)*dx);
%     y(i,1) = 5 * t * (0.2969 * sqrt(x_int * x(i,1)) - ...
%                           0.1260 * (x_int * x(i,1)) - ...
%                           0.3516 * (x_int * x(i,1))^2 + ...
%                           0.2843 * (x_int * x(i,1))^3 - ...
%                           0.1015 * (x_int * x(i,1))^4);
% end
% x(26:40,1) = flipud(x(12:26,1));
% for i = 26:40
%     y(i,1) = -5 * t * (0.2969 * sqrt(x_int * x(i,1)) - ...
%                            0.1260 * (x_int * x(i,1)) - ...
%                            0.3516 * (x_int * x(i,1))^2 + ...
%                            0.2843 * (x_int * x(i,1))^3 - ...
%                            0.1015 * (x_int * x(i,1))^4);
% end
% 
% for i = 41:51
%     x(i,1) = x(i-1,1)+(x(i-1,1)-x(i-2,1))*xsf;
% end
% for i = 1:11
%     x(i,1) = x(52-i,1);
% end
% %% i max
% y(51,2) = 0.02;  
% for j = 3:26    
%     y(51,j) = y(51,j-1) + (y(51,j-1) - y(51,j-2))*ysf;
% end
% x(51,2:end) = x(51,2:end) + x(51,1);
% %% i min
% x(1,2:end)  = x(1,2:end)  + x(end,1);
% y(1,:)  = -y(51,:);
% 
% 
% %% j max
% 
% R_max = y(end,end);
% L = x(1,1)*2+ pi*R_max;
% Delta_L = L/51;
% n_x_line = ceil(x(1,1)/Delta_L);
% x(1,end) = x(1,1);
% 
% for i = 1:(n_x_line)
%     x(i+1,end) = x(1,1) - Delta_L * i;
%     y(i+1,end) = -y(end,end);
% end
% n_arc = 51-n_x_line*2;
% d_theta = pi/n_arc;
% j=1;
% for i = n_x_line+2:(n_arc+n_x_line+1)
%     x(i,end) = R_max*sin(pi+d_theta*j);
%     y(i,end) = R_max*cos(pi+d_theta*j);
%     j=j+1;
% end
% 
% j = 1;
% for i = (n_arc+n_x_line+2):length(x(:,1))
%     x(i,end) = Delta_L * j;
%     y(i,end) = y(end,end);
%     j = j+1;
% end
% x(end,end) = x(1,end);
% 
% %% interp
% for i = 2:(i_max-1)
%     for j = 2:(j_max-1)
%         m_x = (x(i,j_max)-x(i,1))/25;
%         m_y = (y(i,j_max)-y(i,1))/25;
%         x(i,j) = x(i,1)+m_x*(j-1);
%         y(i,j) = y(i,1)+m_y*(j-1);
%     end
% end

%% scattering
load output.txt
x = output(:,1);
y = output(:,2);
x_mat_r = zeros(51,26);
y_mat_r = zeros(51,26);
k = 1;
for i=1:1:26
    x_mat_r(:,i)=x(k:k+50,1);
    y_mat_r(:,i)=y(k:k+50,1);
    k = k+51;
end
figure
hold on
grid on
plot(x_mat(:,1), y_mat(:,1), 'k', 'DisplayName', 'j min');
plot(x_mat(end,:), y_mat(end,:), 'm', 'DisplayName', 'i max');
plot(x_mat(1,:), y_mat(1,:), 'g', 'DisplayName', 'i min');
plot(x_mat(:,26), y_mat(:,26), 'b', 'DisplayName', 'j max');
for j = 2:1:25
    plot(x_mat(:,j), y_mat(:,j),'r')
end
for i = 2:1:50
    plot(x_mat(i,:), y_mat(i,:),'r')
end

