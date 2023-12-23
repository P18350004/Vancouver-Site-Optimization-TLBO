This manuscript presents the applicabilty of TLBO-KRSM optimization algorithm to estimate site-specific Vs profile along with the Vs profile from other optimization algorithms (PSO, ABCO, GA, DE) using MASW techniques.

Important: The matlab code of TLBO algorithm is given here.
           The explanation of each matlab code are written infront of respective lines.
           Repository "Vancouver-Site-Optimization-TLBO" consist of supporitng files.

           
Please adopt the following steps to run the code.

Guideline for quick test of Vancouver Site data
Step-1:  Save Matlab files "c_curve0.mat" and "lambda_curve0" with the Phase velocity and waveleght vector of Vancouver site from the experimental fundamental mode.
Step-2: Open folers "main_TLBO.m" in which lower and upper search spaces are defined as given in the manuscript.
Step-3: Input the parameters of TLBO algorithm as per given in manuscript. 
Step-4: Run the "main_TLBO.m" code and wait to complete 100 iterations.
Step-5: At 100th iteration, Best solutions (thickness and Vs) are obtained. Further, these best solutions are used again in "MASWaves.m" code, which is given in the Repository "Vancouver-Site-Optimization-TLBO".
