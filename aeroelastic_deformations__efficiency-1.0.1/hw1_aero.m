%HW1 - Guy Gordon Biberstein 208010439 & Tamir gazit 208002931
%% Initalization
parm.EI=7.3; %[Nm^2]
parm.GJ=3.7; %[Nm^2]
parm.C_l_alpha=2*pi; %[1/rad]
parm.C_l_delta=0; %[1/rad]
parm.chord_length=100*(10^-3); %[m]
parm.span=600*(10^-3); %[m]
parm.ro=1.225; %[kg/m^3]
parm.v=10; %[m/sec]
parm.q=(1/2)*parm.ro*parm.v^2; %[m]
parm.ac_chord=0.25*parm.chord_length; %[m]
parm.elastic_chord=0.3*parm.chord_length; %[m]
parm.e=abs(parm.elastic_chord-parm.ac_chord); %[m]


% y=0:delta_span:parm.span;

%%
prompt = 'Insert swept angle in deg';
x = input(prompt);
parm.swept=deg2rad(x);


% [big_AFC,big_K,big_TAS,big_T_alpha_S] = matrix_builder(parm,delta_span,N);


%% Div velocity Vs. N

N_arr=[50 75 100 250 300];
v_div_for_N=zeros(length(N_arr),1);
for i=1:length(N_arr)
    
    N=N_arr(i);
    delta_span=parm.span/N;
    [big_AFC,big_K,big_TAS,big_T_alpha_S] = matrix_builder(parm,delta_span,N);
    [v_div_for_N(i),  ] = div_velocity(big_T_alpha_S,big_K,big_TAS,big_AFC,parm);
    
    
end
figure();
plot(N_arr,v_div_for_N,'-o')
xlabel('N - number of strips')
ylabel('Velocity [m/sec]')
grid minor
title(['Divergence velocity Vs. number of strips - \Lambda=',num2str(rad2deg(parm.swept)),''])

figure();
semilogy(N_arr,v_div_for_N-v_div_for_N(end),'-o')
xlabel('N - number of strips')
ylabel('Velocity error [m/sec]')
grid minor
title(['Divergence velocity error Vs. number of strips - \Lambda=',num2str(rad2deg(parm.swept)),''])



%% divergence speed vary with the location of the elastic axis
N=100;
delta_span=parm.span/N;
e_arr=0:0.05:1;
v_div_for_e_location=zeros(length(e_arr),1);
for i=1:length(e_arr)
    parm.elastic_chord=e_arr(i)*parm.chord_length; %[m]
    parm.e=abs(parm.elastic_chord-parm.ac_chord); %[m]
    
    [big_AFC,big_K,big_TAS,big_T_alpha_S] = matrix_builder(parm,delta_span,N);
    [v_div_for_e_location(i), ] = div_velocity(big_T_alpha_S,big_K,big_TAS,big_AFC,parm);
    
    
end
figure();
subplot(1,2,1)
plot(e_arr,v_div_for_e_location,'-o')
xlabel('elastic axis location [% of chord]')
ylabel('Velocity [m/sec]')
grid minor
title(['Divergence velocity Vs. elastic axis location - \Lambda=',num2str(rad2deg(parm.swept))])
xline(0.25,'r--');

subplot(1,2,2)
semilogy(e_arr,v_div_for_e_location,'-o')
xlabel('elastic axis location [% of chord]')
ylabel('Velocity [m/sec] (log scale)')
grid minor
title(['Divergence velocity (semilog) - \Lambda=',num2str(rad2deg(parm.swept))])
xline(0.25,'r--');

% restore the default elastic axis for the sections below
parm.elastic_chord = 0.3*parm.chord_length;
parm.e = abs(parm.elastic_chord - parm.ac_chord);


%% wing tip deformation and torsion angle at an angle of attack of 5 deg in a range
%of wind-tunnel speeds between 0 to 30m/s.
parm.elastic_chord = 0.3*parm.chord_length;
parm.e = abs(parm.elastic_chord - parm.ac_chord);
N=100;
parm.alpha=deg2rad(5); %[rad]

v_arr=0:2:30;
delta_y_tip=zeros(length(v_arr),1);
theta_z_tip=zeros(length(v_arr),1);
theta_x_tip=zeros(length(v_arr),1);
for i=1:1:length(v_arr)
    parm.v=v_arr(i);
    
    [delta_y,theta_z,theta_x,~,~] = deformation(N,parm);
    delta_y_tip(i)=delta_y(end);
    theta_z_tip(i)=theta_z(end);
    theta_x_tip(i)=theta_x(end);
    
end

figure();
subplot(1,2,1)
plot(v_arr,delta_y_tip);
xlabel('v - [m/sec]');
ylabel('deformation - [m]')
title(['Wing tip deformation [m] Vs. Velocity [m/sec] - \alpha=5 - \Lambda=',num2str(rad2deg(parm.swept)),'']);
grid minor


subplot(1,2,2)
plot(v_arr,rad2deg(theta_x_tip));
xlabel('v - [m/sec]');
ylabel('Torsion angle - [deg]')
title(['Wing tip torsion angle [deg] Vs. Velocity [m/sec] - \alpha=5 - \Lambda=',num2str(rad2deg(parm.swept)),''])
grid minor

%%  aeroelastic lift efficiency in a range of wind-tunnel speeds between 0 to 30m/s.
N=100;
parm.alpha=deg2rad(5); %[rad]


v_arr=0:2:30;
eta=zeros(length(v_arr),0);

s=parm.chord_length*delta_span;
for i=1:1:length(v_arr)
    parm.v=v_arr(i);
    [~,~,~,deform,v_div] = deformation(N,parm);
    
    q=(1/2)*parm.ro*v_arr(i)^2;
    q_div=(1/2)*parm.ro*v_div^2;
    k_thetha=q_div*parm.e*parm.C_l_alpha*s;
    eta(i) = 1/(1 - q/q_div);   % same formula as k/(k - q*s*Cla*e); equals 1 if q_div = Inf
end

figure()
plot(v_arr,eta)
xlabel('v - [m/sec]')
ylabel('\eta ')
title(['Aeroelastic lift efficiency Vs. velocity [m/sec] - \Lambda=',num2str(rad2deg(parm.swept)),''])
grid minor

% diffrent swept angle 
figure()
swept=[-20 0 20]*pi/180;

for i=1:1:3
    parm.swept=swept(i);
    for i=1:1:length(v_arr)
        parm.v=v_arr(i);
        [~,~,~,deform,v_div] = deformation(N,parm);
        
        q=(1/2)*parm.ro*v_arr(i)^2;
        q_div=(1/2)*parm.ro*v_div^2;
        k_thetha=q_div*parm.e*parm.C_l_alpha*s;
        eta(i) =k_thetha/(k_thetha-q*s*parm.C_l_alpha*parm.e);
    end
    
    hold on;
    plot(v_arr,eta)
end
    xlabel('v - [m/sec]')
    ylabel('\eta ')
    title('Aeroelastic lift efficiency for diffrent swept angles Vs. velocity [m/sec]')
    grid minor
    legend('-20','0','20');



    %%
res = plot_wing_deformations(parm, 30, [-10 0 20], 0.40);   % elastic axis at 40% chord   % speed 30 m/s
[s,dy,thx,thz,v_div] = wing_state(parm,100,30,15);    % any sweep, any speed