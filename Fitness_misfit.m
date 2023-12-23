function e = Fitness_misfit(x)


% Import sample data from fundamental mode dispersion curve
load('lambda_curve0.mat')  %Wavelengths of fundamantal mode
load('c_curve0.mat')      % Rayleigh wave phase velocities of fundamental mode

%% inversion

c_test_min = 0; % minimum Rayleigh wave phase velocity (m/s)
c_test_max = 700;  %maximum Rayleigh wave phase velocity (m/s)
delta_c_test = 5; % smallest interval of Rayleigh wave phase velocity (m/s)
c_test = c_test_min:delta_c_test:c_test_max; % m/s

%inputs
n = 5;  % number of soil layers
alpha = [ 1000 1000 1000 1000 1000 1000]; % P-wave velocity (m/s)
%this population will be simulated
h(1,:)=[x(1) x(2) x(3) x(4) x(5) x(6)];   %soil layer thickness
h(2,:)=[ x(7) x(8) x(9) x(10) x(11) x(12)];  %shear wave velocity of soil layers
rho = [2000 2000 2000 2000 2000 2000]; % density of soil layers  (kg/m^3)


% Compute wave numbers that correspond to wavelengths lambda
k = (2*pi)./lambda_curve0;

D = zeros(length(c_test),length(k));
c_t = zeros(length(k),1);

% For each wave number, recompute the system stiffness matrix using
% different values of c_test until its determinant has a sign change. 

for l = 1:length(k)
    for m = 1:length(c_test)
        %stiffness matrix
                        % System stiffness matrix
                K = zeros(2*(n+1),2*(n+1));

                % Check to see if the trial phase velocity is equal to the shear wave velocity
                % or compression wave velocity of one of the layers
                epsilon = 0.0001;
                while any(abs(c_test(m)-h(2,n))<epsilon) || any(abs(c_test(m)-alpha(n))<epsilon)
                    c_test(m) = c_test(m)*(1-epsilon);
                end

                % Finite thickness layers j = 1,...,n
                for j = 1:n
                    % Compute element stiffness matrix for layer j
                    
                    %ke layer function
                        r(j) = sqrt(1-c_test(m)^2/alpha(j)^2);
                        s(j) = sqrt(1-c_test(m)^2/h(2,j)^2);

                        Cr(j) = cosh(k(l)*r(j)*h(1,j));
                        Sr(j) = sinh(k(l)*r(j)*h(1,j));
                        Cs(j) = cosh(k(l)*s(j)*h(1,j));
                        Ss(j) = sinh(k(l)*s(j)*h(1,j));

                        D(j) = 2*(1-Cr(j)*Cs(j)) + (1/(r(j)*s(j)) + r(j)*s(j))*Sr(j)*Ss(j);

                        k11_e(j) = (k(l)*rho(j)*c_test(m)^2)/D(j) * (s(j)^(-1)*Cr(j)*Ss(j) - r(j)*Sr(j)*Cs(j));
                        k12_e(j) = (k(l)*rho(j)*c_test(m)^2)/D(j) * (Cr(j)*Cs(j) - r(j)*s(j)*Sr(j)*Ss(j) - 1) - k(l)*rho(j)*h(2,j)^2*(1+s(j)^2);
                        k13_e(j) = (k(l)*rho(j)*c_test(m)^2)/D(j) * (r(j)*Sr(j) - s(j)^(-1)*Ss(j));
                        k14_e(j) = (k(l)*rho(j)*c_test(m)^2)/D(j) * (-Cr(j) + Cs(j));
                        k21_e(j) = k12_e(j);
                        k22_e(j) = (k(l)*rho(j)*c_test(m)^2)/D(j) * (r(j)^(-1)*Sr(j)*Cs(j) - s(j)*Cr(j)*Ss(j));
                        k23_e(j) = -k14_e(j);
                        k24_e(j) = (k(l)*rho(j)*c_test(m)^2)/D(j) * (-r(j)^(-1)*Sr(j) + s(j)*Ss(j));
                        k31_e(j) = k13_e(j);
                        k32_e(j) = k23_e(j);
                        k33_e(j) = k11_e(j);
                        k34_e(j) = -k12_e(j);
                        k41_e(j) = k14_e(j);
                        k42_e(j) = k24_e(j);
                        k43_e(j) = -k21_e(j);
                        k44_e(j) = k22_e(j);

                        Ke = [k11_e(j) k12_e(j) k13_e(j) k14_e(j);
                            k21_e(j) k22_e(j) k23_e(j) k24_e(j);
                            k31_e(j) k32_e(j) k33_e(j) k34_e(j);
                            k41_e(j) k42_e(j) k43_e(j) k44_e(j)];
                         % Add to the system stiffness matrix
                            DOFS = [(2*j-1):(2*j+2)];
                            K(DOFS,DOFS) = K(DOFS,DOFS)+Ke;
                end
             %half space
                            r(j+1) = sqrt(1-c_test(m)^2/alpha(j+1)^2);
                            s(j+1) = sqrt(1-c_test(m)^2/h(2,j+1)^2);

                            k_11(j+1) = k(l)*rho(j+1)*h(2,j+1)^2*(r(j+1)*(1-s(j+1)^2))/(1-r(j+1)*s(j+1));
                            k_12(j+1) = k(l)*rho(j+1)*h(2,j+1)^2*(1-s(j+1)^2)/(1-r(j+1)*s(j+1)) - 2*k(l)*rho(j+1)*h(2,j+1)^2;
                            k_21(j+1) = k_12(j+1);
                            k_22(j+1) = k(l)*rho(j+1)*h(2,j+1)^2*(s(j+1)*(1-s(j+1)^2))/(1-r(j+1)*s(j+1));

                            Ke_halfspace = [k_11(j+1) k_12(j+1) ; k_21(j+1) k_22(j+1)];
                    % Add to the system stiffness matrix
                    DOFS = [2*(n+1)-1:2*(n+1)];
                    K(DOFS,DOFS) = K(DOFS,DOFS)+Ke_halfspace;

                    % Evaluate determinant of system stiffness matrix
                    D(l,m) = real(det(K));
        if m==1;
            sign_old = sign(D(l,m));
        else
            sign_old = signD;
        end
        signD = sign(D(l,m));
        if sign_old*signD == -1
            c_t(l)=c_test(m);
            lambda_t(l)=2*pi/k(l);
            break
        end
    end 
end

%%misfit
Q = length(c_t); % Number of data points in theoretical and experimental dispersion curves.
c_curve0;
temp = 0;

for q = 1:Q
    temp = temp + sqrt((c_curve0(q) - c_t(q))^2)/c_curve0(q);
end

e = 1/Q * temp * 100;

end
