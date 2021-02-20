clear;
close all;

%--------INITIALIZATION OF VARIABLES-----------------%

%--Normal variables---%
f=50;                       % frequency
fs=2*16*50;                 % sampling frequency should be atleast twice the highest frequency of signal to be analysed
T_cycle = 1/f;             % time period for one cycle
T_max=2*T_cycle;            % max time period
T_vector=0:1/fs:T_max;      % time vector
T_plot = 1000*T_vector;     % Time vector in milli seconds in order to use during plotting


%----Amplitudes and Phases of the Artificially simulated cosine signal-----%
% First harmonic phase and Amplitude
A1    = 1.0;
phi_1 = 30;
% Third harmonic phase and Amplitude
A3    = 0.1;
phi_3 = 210;
%Fifth harmonic phase and Amplitude
A5    = 0.06;
phi_5 = 180;
%Ninth harmonic phase and Amplitude
A9    = 0.009;
phi_9 = -145;
%Eleventh harmonic phase and Amplitude
A11    = 0.005;
phi_11 = 30;
%Thirteen harmonic phase and Amplitude
A13    = 0.003;
phi_13 = 0;

%---Matrices of state space equation and kalman filter algorithm----%
F = eye(12,12);             %state feed back matrix. In this case it is identity matrix
Q = 0.05*eye(12,12);      %Input white noise covariance matrix
R       = 0.05;           %Measurement white noise covariance matrix

X_0_minus = zeros(12,1);    %Initial state vector for the algorithm
P_0_minus = 10*eye(12,12);  %Initial error cavariance matrix for the algorithm

%------- Computation of H matrix----------------------------------%
H_k=calculate_h(f,fs);          % circular H_k is computed

%-----Computation of Gain K_k for each state----------------------%
[K_k] = calculate_k(F,H_k,Q,R,P_0_minus,fs,T_vector);

K_1 = K_k(:,1);
K_2 = K_k(:,2);
K_3 = K_k(:,3);
K_4 = K_k(:,4);
K_5 = K_k(:,5);
K_6 = K_k(:,6);
K_7 = K_k(:,7);
K_8 = K_k(:,8);
K_9 = K_k(:,9);
K_10 = K_k(:,10);
K_11 = K_k(:,11);
K_12 = K_k(:,12);

%---Generate output vector i.e frabricate a cosine wave of known harmonic content in it--%
Z = generate_signal(A1,phi_1,A3,phi_3,A5,phi_5,A9,phi_9,A11,phi_11,A13,phi_13,T_vector,f);

%----Estimate the states using the kalman filtering algorithm-----%
[X_k] = calculate_states(F,H_k,Z,K_k,T_vector,X_0_minus);

X_1 = X_k(:,1);
X_2 = X_k(:,2);
X_3 = X_k(:,3);
X_4 = X_k(:,4);
X_5 = X_k(:,5);
X_6 = X_k(:,6);
X_7 = X_k(:,7);
X_8 = X_k(:,8);
X_9 = X_k(:,9);
X_10 = X_k(:,10);
X_11 = X_k(:,11);
X_12 = X_k(:,12);

%------Harmonic content calculation from the estimated states------%
Harmonic_1 = (X_1.^2 + X_2.^2).^0.5;
Harmonic_3 = (X_3.^2 + X_4.^2).^0.5;
Harmonic_5 = (X_5.^2 + X_6.^2).^0.5;
Harmonic_9 = (X_7.^2 + X_8.^2).^0.5;
Harmonic_11 = (X_9.^2 + X_10.^2).^0.5;
Harmonic_13 = (X_11.^2 + X_12.^2).^0.5;

%-----Total Harmonic Distortion Calculation-----%
THD_theory = (A3^2+A5^2+A9^2+A11^2+A13^2)^0.5/A1;
THD_estimated = (Harmonic_3(end)^2+Harmonic_5(end)^2+Harmonic_9(end)^2+...
                    Harmonic_11(end)^2+Harmonic_13(end)^2)^0.5/Harmonic_1(end);
                
disp('Theoritical THD =');disp(THD_theory);
disp('Estimated THD =');disp(THD_estimated);

%------Plotting of Gains and Harmonics----------------%
%Plotting of fabricated cosine wave
figure(1)
plot(T_plot,Z,'r','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('Time----> in milli seconds');
ylabel('Magnitude in P.U');
title('Artificially created cosine wave');
grid on;

%Plotting K_ks of each state
figure(2);
subplot(3,2,1)
plot(T_plot,K_1,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_1');
grid on;

subplot(3,2,2)
plot(T_plot,K_2,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_2')
grid on;

subplot(3,2,3)
plot(T_plot,K_3,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_3')
grid on;

subplot(3,2,4)
plot(T_plot,K_4,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_4')
grid on;

subplot(3,2,5)
plot(T_plot,K_5,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_5')
grid on;

subplot(3,2,6)
plot(T_plot,K_6,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_6')
grid on;

figure(3);
subplot(3,2,1)
plot(T_plot,K_7,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_7');
grid on;

subplot(3,2,2)
plot(T_plot,K_8,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_8')
grid on;

subplot(3,2,3)
plot(T_plot,K_9,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_9')
grid on;

subplot(3,2,4)
plot(T_plot,K_10,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_1_0')
grid on;

subplot(3,2,5)
plot(T_plot,K_11,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_1_1')
grid on;

subplot(3,2,6)
plot(T_plot,K_12,'-s','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
xlabel('T----> milli seconds');
ylabel('Gain');
title('Kalman Gain of X_1_2')
grid on;

%--------- Plotting of harmonics-------------------
figure(4);
subplot(3,2,1)
plot(T_plot,Harmonic_1,'-r','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
hold on;
plot(T_plot,A1,'ks','LineWidth',1,'MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('Time----> in milli seconds');
ylabel('Magnitude in P.U');
title('1^s^t harmonic');
grid on;

subplot(3,2,2)
plot(T_plot,Harmonic_3,'-r','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
hold on;
plot(T_plot,A3,'ks','LineWidth',1,'MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('Time----> in milli seconds');
ylabel('Magnitude in P.U');
title('3^r^d harmonic');
grid on;

subplot(3,2,3)
plot(T_plot,Harmonic_5,'-r','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
hold on;
plot(T_plot,A5,'ks','LineWidth',1,'MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('Time----> in milli seconds');
ylabel('Magnitude in P.U');
title('5^t^h harmonic');
grid on;

subplot(3,2,4)
plot(T_plot,Harmonic_9,'-r','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
hold on;
plot(T_plot,A9,'ks','LineWidth',1,'MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('Time----> in milli seconds');
ylabel('Magnitude in P.U');
title('9^t^h harmonic');
grid on;

subplot(3,2,5)
plot(T_plot,Harmonic_11,'-r','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
hold on;
plot(T_plot,A11,'ks','LineWidth',1,'MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('Time----> in milli seconds');
ylabel('Magnitude in P.U');
title('11^t^h harmonic');
grid on;

subplot(3,2,6)
plot(T_plot,Harmonic_13,'-r','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','w');
hold on;
plot(T_plot,A13,'ks','LineWidth',1,'MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('Time----> in milli seconds');
ylabel('Magnitude in P.U');
title('13^t^h harmonic');
grid on;