function [v_div,SIC,V] = div_velocity(big_T_alpha_S,big_K,big_TAS,big_AFC,parm)
%DIV_VELOCITY Summary of this function goes here
%   Detailed explanation goes here
SIC=big_T_alpha_S*(big_K\big_TAS');
[V,D]=eig(big_AFC*SIC);
ev=real(diag(D));            % drop numerical imaginary noise
lamda=max(ev);               % largest signed eigenvalue
if lamda <= 0
    v_div=Inf;               % no divergence
else
    q_d=1/lamda;
    v_div=sqrt((2*q_d)/parm.ro);
end

end

