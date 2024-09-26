clear, clc;
ni = 51;
nj = 26;
k = 1;
% %%
% output = zeros(ni*nj,2);
% for j = 1:nj
%     for i = 1:ni
%         output(k,1) = i;
%         output(k,2) = j;
%         k = k+1;
%     end
% end
% output(:,1) = output(:,1) - 1; 
% output(:,2) = output(:,2) - 1; 
% %fileID = fopen('E:\Technion\CFD\hw2\cmake-build-debug\Grid_File2.txt', 'w');
% fprintf(fileID, '%d %d %d %d\n', ni, nj, 4, 7 );
% for row = 1:size(output, 1)
%     fprintf(fileID, '%d %d\n', output(row, 1), output(row, 2));
% end
% fclose(fileID);
%system('hw2.exe');
%%
Q = table2array(readtable("Q.txt"));
Q1 = Q(1:ni*nj);
Q2 = Q(ni*nj+1:2*ni*nj);
Q3 = Q(2*ni*nj+1:3*ni*nj);
Q4 = Q(3*ni*nj+1:4*ni*nj);
k = 1;
for i=1:1:nj
    Q1_mat(:,i)=Q1(k:k+ni-1,1);
    Q2_mat(:,i)=Q2(k:k+ni-1,1);
    Q3_mat(:,i)=Q3(k:k+ni-1,1);
    Q4_mat(:,i)=Q4(k:k+ni-1,1);
    k = k+ni;
end
%%
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


