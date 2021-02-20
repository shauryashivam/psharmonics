function H_k = calculate_h(f,fs)
% This is the function to calculate circular H_k values.
omega   = 2*pi*f;
delta_t = 1/fs;

H_k =[];
for k=0:(fs/f)-1
    h = [cos(omega*k*delta_t) -1*sin(omega*k*delta_t) cos(3*omega*k*delta_t) -sin(3*omega*k*delta_t) cos(5*omega*k*delta_t) ...
         -sin(5*omega*k*delta_t) cos(9*omega*k*delta_t) -sin(9*omega*k*delta_t) cos(11*omega*k*delta_t) -sin(11*omega*k*delta_t) ...
         cos(13*omega*k*delta_t) -sin(13*omega*k*delta_t)];
             
     H_k = [H_k; h];
end