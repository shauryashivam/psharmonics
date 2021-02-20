function [K_k] = calculate_k(F,H_k,Q,R,P_0_minus,fs,T_k_vec)
% This function calculates the cyclic K_k

K_k = [];
delta_t = 1/fs;
P_k_minus = P_0_minus;

for n=0:length(T_k_vec)-1
    cyc_n = rem(n,32);          %Inorder to use the cyclic H_k
    dummy_h = H_k(cyc_n+1,:);
    dummy_k = P_k_minus*dummy_h'*inv(dummy_h*P_k_minus*dummy_h'+R);
    K_k = [K_k;dummy_k'];
    P_k_plus = (eye(12,12)-dummy_k*dummy_h)*P_k_minus;
    P_k_minus = F*P_k_plus*F'+Q;
end