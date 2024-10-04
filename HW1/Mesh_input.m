function Mesh_input(control_case)
% Open a file for writing
fileID = fopen('input.txt', 'w');

thickness = 0.12;
imax = 50;
jmax = 25;
tel  = 11;
le = 25;
teu = 39;
dy = 0.02;
xsf = 1.15;
ysf = 1.15;
r = 0.01;
w = 1.5;
% Write the data to the file
fprintf(fileID, '%.2f %.0f %.0f %.0f %.0f %.0f %.2f %.2f %.2f %.2f %.1f %.0f\n', thickness, imax, jmax, tel, le, teu, dy, xsf, ysf, r, w, control_case);

% Close the file
fclose(fileID);
end