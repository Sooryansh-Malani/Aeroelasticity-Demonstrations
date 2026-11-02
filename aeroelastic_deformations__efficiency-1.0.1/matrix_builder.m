function [big_AFC,big_K,big_TAS,big_T_alpha_S] = matrix_builder(parm,delta_span,N)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


%% AFC
delta=zeros(N,1);
alpha=ones(N,1)*deg2rad(5);
S=ones(N,1)*parm.chord_length*delta_span;
C = cell(2);
AFC=cell(N,1);

for i=1:1:N
    C{1, 1}  =S(i)*parm.C_l_alpha;
    C{1, 2}  =S(i)*parm.C_l_delta*delta(i);
    C{2, 1}  =0*alpha(i);
    C{2, 2}  =(-1)*S(i)*parm.C_l_delta*(abs(parm.ac_chord-parm.chord_length))*delta(i);
    
    AFC{i,1}=C;
end


big_AFC=zeros(2*N);
j=1;
for i=1:2:(2*N)
    big_AFC(i,i)=AFC{j, 1}{1, 1};
    big_AFC(i,i+1)=AFC{j, 1}{1, 2};
    big_AFC(i+1,i)=AFC{j, 1}{2, 1};
    big_AFC(i+1,i+1)=AFC{j, 1}{2, 2};
    j=j+1;
end

%% SIC
% KE
KE=[(12*parm.EI)/(delta_span^3) (6*parm.EI)/(delta_span^2) 0 (-12*parm.EI)/(delta_span^3)  (6*parm.EI)/(delta_span^2) 0;
    (6*parm.EI)/(delta_span^2) (4*parm.EI)/(delta_span) 0 (-6*parm.EI)/(delta_span^2) (2*parm.EI)/(delta_span) 0;
    0 0 parm.GJ/delta_span 0 0 -parm.GJ/delta_span;
    (-12*parm.EI)/(delta_span^3) (-6*parm.EI)/(delta_span^2) 0 (12*parm.EI)/(delta_span^3)  (-6*parm.EI)/(delta_span^2) 0;
    (6*parm.EI)/(delta_span^2) (2*parm.EI)/(delta_span) 0 (-6*parm.EI)/(delta_span^2) (4*parm.EI)/(delta_span) 0;
    0 0 -parm.GJ/delta_span 0 0 parm.GJ/delta_span;
    ];

big_K=zeros(3*N);
big_K(1:3,1:3)=KE(4:6,4:6);
for i=1:3:3*(N-1)
    big_K(i:i+5,i:i+5)=big_K(i:i+5,i:i+5)+KE(:,:);
end

% TAS

C = cell(2,6);
TAS=cell(N,0);
big_TAS=zeros(2*N,3*N);
for i=1:N
    e1=parm.e-(delta_span/2)*sin(parm.swept);
    e2=parm.e+(delta_span/2)*sin(parm.swept);
    C{1,1}=1;
    C{1,2}=-e1*sin(parm.swept);
    C{1,3}=e1*cos(parm.swept);
    C{1,4}=1;
    C{1,5}=-e2*sin(parm.swept);
    C{1,6}=e2*cos(parm.swept);
    C{2,1}=0;
    C{2,2}=-sin(parm.swept);
    C{2,3}=cos(parm.swept);
    C{2,4}=0;
    C{2,5}=-sin(parm.swept);
    C{2,6}=cos(parm.swept);
    
    TAS{i,1}=C;
end

big_TAS(1,1)=TAS{1, 1}{1, 4};
big_TAS(1,2)=TAS{1, 1}{1, 5};
big_TAS(1,3)=TAS{1, 1}{1, 6};
big_TAS(2,1)=TAS{1, 1}{2, 4};
big_TAS(2,2)=TAS{1, 1}{2, 5};
big_TAS(2,3)=TAS{1, 1}{2, 6};

k=2;
j=1;
t=1;
for i=3:1:(2*N)
    big_TAS(i,j)=TAS{k, 1}{t, 1};
    big_TAS(i,j+1)=TAS{k, 1}{t, 2};
    big_TAS(i,j+2)=TAS{k, 1}{t, 3};
    big_TAS(i,j+3)=TAS{k, 1}{t, 4};
    big_TAS(i,j+4)=TAS{k, 1}{t, 5};
    big_TAS(i,j+5)=TAS{k, 1}{t, 6};
    if(mod(i,2)==0) % if i is odd;
        k=k+1;
        j=j+3;
        t=1;
    else % next i is even
        t=2;
    end
end

% T_alpha_S


C = cell(2,6);
T_alpha_S=cell(N,0);
big_T_alpha_S=zeros(2*N,3*N);
for i=1:N

    C{1,1}=0;
    C{1,2}=-sin(parm.swept);
    C{1,3}=cos(parm.swept);
    C{1,4}=0;
    C{1,5}=-sin(parm.swept);
    C{1,6}=cos(parm.swept);
    C{2,1}=0;
    C{2,2}=0;
    C{2,3}=0;
    C{2,4}=0;
    C{2,5}=0;
    C{2,6}=0;
    
    T_alpha_S{i,1}=C;
end

big_T_alpha_S(1,1)=T_alpha_S{1, 1}{1, 4};
big_T_alpha_S(1,2)=T_alpha_S{1, 1}{1, 5};
big_T_alpha_S(1,3)=T_alpha_S{1, 1}{1, 6};
big_T_alpha_S(2,1)=T_alpha_S{1, 1}{2, 4};
big_T_alpha_S(2,2)=T_alpha_S{1, 1}{2, 5};
big_T_alpha_S(2,3)=T_alpha_S{1, 1}{2, 6};

k=2;
j=1;
t=1;
for i=3:1:(2*N)
    big_T_alpha_S(i,j)=T_alpha_S{k, 1}{t, 1};
    big_T_alpha_S(i,j+1)=T_alpha_S{k, 1}{t, 2};
    big_T_alpha_S(i,j+2)=T_alpha_S{k, 1}{t, 3};
    big_T_alpha_S(i,j+3)=T_alpha_S{k, 1}{t, 4};
    big_T_alpha_S(i,j+4)=T_alpha_S{k, 1}{t, 5};
    big_T_alpha_S(i,j+5)=T_alpha_S{k, 1}{t, 6};
    if(mod(i,2)==0) % if i is odd;
        k=k+1;
        j=j+3;
        t=1;
    else % next i is even
        t=2;
    end
end
big_TAS = 0.5*big_TAS;
big_T_alpha_S = 0.5*big_T_alpha_S;
end

