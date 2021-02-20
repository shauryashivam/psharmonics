function [X_k] = calculate_states(F,H_k,Z,K_k,T_vector,X_0_minus,fs,f)
% This function calculates the cyclic K_k

X_k = [];
X_k_minus = X_0_minus;
for n=0:length(T_vector)-1
    cyc_n = rem(n,fs/f);
    dummy_h = H_k(cyc_n+1,:);
    dummy_k = K_k(n+1,:)';
    dummy_z = Z(n+1);
    X_k_plus = X_k_minus + dummy_k*(dummy_z - dummy_h*X_k_minus);
    X_k =[X_k;X_k_plus'];
    X_k_minus = F*X_k_plus;
end