function [delta_y,theta_z,theta_x,deform,v_div] = deformation(N,parm)
%DEFORMATION Summary of this function goes here
%   Detailed explanation goes here
parm.q=(1/2)*parm.ro*parm.v^2; %[m]
delta_span=parm.span/N;
[big_AFC,big_K,big_TAS,big_T_alpha_S] = matrix_builder(parm,delta_span,N);
[v_div,SIC,~] = div_velocity(big_T_alpha_S,big_K,big_TAS,big_AFC,parm);


deform.alpha_rigid=ones(2*N,1)*parm.alpha;
deform.F_top_rigid=parm.q*big_AFC*deform.alpha_rigid;
deform.F_top=inv((eye(2*N)-parm.q*big_AFC*SIC))*deform.F_top_rigid;
deform.F=big_TAS'*deform.F_top;
deform.u=inv(big_K)*deform.F;

j=1;
for i=1:3:3*N
    delta_y(j)=deform.u(i);
    theta_z(j)=deform.u(i+1);
    theta_x(j)=deform.u(i+2);
    j=j+1;
end
end

