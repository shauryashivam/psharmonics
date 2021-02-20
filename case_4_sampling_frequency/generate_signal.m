function s = generate_signal(A1,phi_1,A3,phi_3,A5,phi_5,A9,phi_9,A11,phi_11,A13,phi_13,t_vector,f)

%This function gives out the signal of known harmonic content when the
%parameters are passed. The input arguments should be amplitudes and phases
%of the harmonics. time vector and frequency also should be passed as input.

phi_1 = (phi_1*pi)/180;
phi_3 = (phi_3*pi)/180;
phi_5 = (phi_5*pi)/180;
phi_9 = (phi_9*pi)/180;
phi_11 = (phi_11*pi)/180;
phi_13 = (phi_13*pi)/180;


omega=2*pi*f;

s=A1*cos(omega*t_vector + phi_1)+ A3*cos(3*omega*t_vector + phi_3)+A5*cos(5*omega*t_vector + phi_5)+ ...
     A9*cos(9*omega*t_vector + phi_9)+A11*cos(11*omega*t_vector + phi_11)+A13*cos(13*omega*t_vector + phi_13);
