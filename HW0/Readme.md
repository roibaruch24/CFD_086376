# HW 0
The HW0 assignment is solving a second-order differential equation using numerical methods.

The program should only run the Matlab HW0.m script. The script gets the initial values: 

N, the boundary conditions, and the interval [a b]. Then, the script writes the input.txt, compiles, and runs the C program that solves the equation.

Then, all the relevant graphs are plotted.

The code not using the Matlab script can be found here: https://github.com/roibaruch24/CFD_086376/commit/96ea68fa531ecec515f5b97705841f9b3c92e146?diff=unified&w=0# 

Press expand all to see the entire code:![image](https://github.com/roibaruch24/CFD_086376/assets/171028386/2861fe3c-6229-4d78-80ff-038450814d03)


For new users, either use the non-Matlab (Recommended for easy use) version linked above or follow these steps to be able to run the C program through Matlab:

## 1. Download the Repository:
If you don't have Git installed, you can download the repository as a ZIP file from GitHub:

Go to the CFD_086376 GitHub repository.
Click on the "Code" button.
Select "Download ZIP".
Extract the ZIP file's contents to a directory on your local machine, e.g., C:\Users\YourUsername\Documents\CFD_086376.
## 2. Set Up MATLAB
Ensure MATLAB is set up with a C compiler. You can configure this using:
 ```bash
mex -setup
```
## 3. Edit JSON file
Locate the JSON file for the compiler you want to configure, now edit it, and add in the includePath as follows:
```bash
{
    "configurations": [
        {
            "name": "windows-gcc-x86",
            "includePath": [
                "${workspaceFolder}/**",
                "C:/Program Files/MATLAB/R2023a/extern/include"
            ],
            "compilerPath": "C:/MinGW/bin/gcc.exe",
            "cStandard": "${default}",
            "cppStandard": "${default}",
            "intelliSenseMode": "windows-gcc-x86",
            "compilerArgs": [
                ""
            ]
        }
    ],
    "version": 4
```
## 4. Run HW0.m
Ensure all your files are in the same directory, add the input.txt file to the folder and change the ID in HW0.m to the correct directory.



