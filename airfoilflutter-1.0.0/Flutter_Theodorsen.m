close all;
clear all;

mu = 500;  % airfoil density (kg/m3)
c = 1;  % chord length (m)
x_f = 0.5 * c; % Location of elastical axis measured from leading edge of airfoil
a = -0.0 * c;  % Location of elastic axis measured from center of airfoil
t = 0.01;   % Airfoil thickness (m)
b = c / 2;  % Airfoil mid chord (m)



omega_h = 5 * 2 * pi;    % Plunge natural frequency (rad/s)
omega_alpha = 17* 2 * pi;   % Pitch natural frequency (rad/s)

m = mu * c * t; % Mass per unit span of airfoil (kg/m)

S = -m * a;  

K_h = m * (omega_h ^ 2);    % Plunge stiffness

I_c = (1 / 12) * m * (c ^ 2);   % Moment of inertia of airfoil about center of airfoil
I_alpha = I_c + (m * (a ^ 2));  % Parallel axis theorem to obtain moment of inertia about elastic axis

K_alpha = I_alpha * (omega_alpha ^ 2);  % Pitch stiffness

rho = 1.225; % Density of air (kg/m3)

V1 = [];    %plunge eigenvalue
d1 = [];    % plunge damping
V2 = [];    % Pitch eigenvalue
d2 = [];    % Pitch damping
d_all = [];

U_min = 0.1;  % (m/s)
U_max = 80; % (m/s)
dU = 0.1; % (m/s)

tolerance = 1e-3;


U_range = U_min:dU:U_max;   % Range of velocity values from 0 to 80 m/s in steps of 0.1 m/s


for U = U_range

    ii = 3; % Eigenvalue index to extract

    for omega = [sqrt(K_h / m) sqrt(K_alpha / I_alpha)] % Initial estimate of flutter frequency

        error = 10; % Arbitrary number to determine convergence to stop p-k method iteration

        while error > tolerance     % While loop for eigenvalue iteration, stops when eigenvalue converges

            k = (omega * b) / U;
            C_k = 1 - (0.165 ./ (1 - 1i * (0.0455 ./ k))) - (0.335 ./ (1 - 1i * (0.30 ./ k)));  % Reduced frequency k using Theodorsen Function approximation

            M_matrix = [m, S; S, I_alpha];
            K_matrix = [K_h, 0; 0, K_alpha];
            D_matrix = rho * pi * (b ^ 2) * [-1, a; a, -((a ^ 2) + ((b ^ 2) / 8))];
            E_k_matrix = rho * pi * b * U * [-2 * C_k, -b + (2 * C_k * (a - (b / 2))); b - (b - (2 * a + b) * C_k), (b - (2 * a + b) * C_k) * (a - (b / 2))];
            F_k_matrix = rho * pi * b * (U ^ 2) * [0, -2 * C_k; 0, b - (b - (2 * a + b) * C_k)];
       
            A_k = [zeros(2), eye(2); (M_matrix - D_matrix) \ (F_k_matrix - K_matrix), (M_matrix - D_matrix) \ E_k_matrix];

            Eigenvalues = eig(A_k); % Eigenvalue analysis

            [Eigenvalues_Imag, Index] = sort(imag(Eigenvalues), 'ascend');  % Arranging the eigenvalues in ascending order by their imaginary parts

            Eigenvalue_Imag_1 = Eigenvalues_Imag(ii);   % Extracting imaginary part of largest/2nd largest eigenvalue for checking against initial assumed frequency (omega)
            error = abs(Eigenvalue_Imag_1 - omega);     % Calculate error between imaginary part of eigenvalue against initial assumed frequency (omega)
            omega = Eigenvalue_Imag_1;  % Initial assumed frequency is updated to the new eigenvalue imaginary part

        end

        Eigenvalues_sorted = Eigenvalues(Index);

        % Section to obtain and store frequencies and damping

        % plunge values
        if ii == 3

            V1 = [V1 abs(Eigenvalues_sorted(ii))];
            d1 = [d1 -1 * real(Eigenvalues_sorted(ii)) / abs(Eigenvalues_sorted(ii))];

        end

        % pitch values
        if ii == 4

            V2 = [V2 abs(Eigenvalues_sorted(ii))];

            d2 = [d2 -1 * real(Eigenvalues_sorted(ii)) / abs(Eigenvalues_sorted(ii))];


        end

        d_all = [d_all Eigenvalues_sorted];

        ii = ii + 1;

    end

end

% Loop to find the flutter speed and frequency at the point where damping
% goes negative
for i = 1:1:length(U_range)
    
    if d1(i)<0 
        Flutter_speed = U_range(i);
        Flutter_frequency = V1(i)/2/pi;
        Flutter_damping = d1(i);
        break;
    end
    
    if d2(i)<0 
        Flutter_speed = U_range(i);
        Flutter_frequency = V2(i)/2/pi;
        Flutter_damping = d2(i);
        break;
    end
    
end

fprintf("Flutter Speed: %.2f m/s \n", Flutter_speed);
fprintf("Flutter Frequency: %.3f Hz \n", Flutter_frequency);

figure(1)
title("Natural Frequency vs Airspeed")
xlabel("U (m/s)")
ylabel("Natural Frequency (Hz)")
hold on
plot(U_range,V1/2/pi,"DisplayName","P-K Plunge");
plot(U_range,V2/2/pi,"DisplayName","P-K Pitch");
plot(Flutter_speed, Flutter_frequency,"x","DisplayName", "Flutter Point", "Color", "red", "MarkerSize", 15)
legend show
hold off
grid on

figure(2)
title("Damping Ratio vs Airspeed")
xlabel("U (m/s)")
ylabel("Damping Ratio")
hold on
plot(U_range,d1,"DisplayName","P-K Plunge");
plot(U_range,d2,"DisplayName","P-K Pitch");
plot(Flutter_speed, Flutter_damping,"x","DisplayName", "Flutter Point", "Color", "red", "MarkerSize", 15)
legend show
hold off
grid on
%%



























%% Mode visualization at a chosen airspeed
U_vis      = 75;     % airspeed to visualize (m/s)
show_decay = true;   % true: include e^(sigma*t) decay/growth, false: steady oscillation
T_total    = 1;      % time window (s)
n_frames   = 600;    % samples / animation frames

p = struct('m',m,'S',S,'I_alpha',I_alpha,'K_h',K_h,'K_alpha',K_alpha, ...
           'a',a,'b',b,'rho',rho);

% Run the p-k iteration for each mode and keep eigenvalue + eigenvector
omega_start = [sqrt(K_h/m), sqrt(K_alpha/I_alpha)];   % plunge guess, pitch guess
lam = zeros(1,2);
vec = zeros(4,2);
for j = 1:2
    [lam(j), vec(:,j)] = pk_mode(omega_start(j), U_vis, j+2, p, tolerance);

    % Normalise: dominant DOF real and positive, overall size set for visibility
    ref = vec(1,j);
    if b*abs(vec(2,j)) > abs(vec(1,j)), ref = b*vec(2,j); end
    vec(:,j) = vec(:,j) * conj(ref)/abs(ref);
    vec(:,j) = vec(:,j) * (0.25*b) / max(abs(vec(1,j)), b*abs(vec(2,j)));
end

fprintf("\nAt U = %.1f m/s:\n", U_vis);
for j = 1:2
    fprintf("Mode %d: f = %.3f Hz, damping ratio = %.4f, |h|/(b|alpha|) = %.3f\n", ...
        j, imag(lam(j))/2/pi, -real(lam(j))/abs(lam(j)), abs(vec(1,j))/(b*abs(vec(2,j))));
end

% Time histories of each mode: x(t) = Re(v * e^(lambda t))
t = linspace(0, T_total, n_frames);
hh = zeros(3, n_frames);   % rows: mode 1, mode 2, sum of both
aa = zeros(3, n_frames);
for j = 1:2
    if show_decay
        x = real(vec(:,j) * exp(lam(j)*t));
    else
        x = real(vec(:,j) * exp(1i*imag(lam(j))*t));
    end
    hh(j,:) = x(1,:);
    aa(j,:) = x(2,:);
end
hh(3,:) = hh(1,:) + hh(2,:);
aa(3,:) = aa(1,:) + aa(2,:);

labels = {sprintf('Mode 1 (%.2f Hz)', imag(lam(1))/2/pi), ...
          sprintf('Mode 2 (%.2f Hz)', imag(lam(2))/2/pi), ...
          'Mode 1 + Mode 2'};

% Figure 3: time histories
figure(3)
subplot(2,1,1); hold on; grid on
for r = 1:3, plot(t, hh(r,:)/b, "DisplayName", labels{r}); end
ylabel("h / b"); title(sprintf("Plunge and pitch response at U = %.1f m/s", U_vis)); legend show
subplot(2,1,2); hold on; grid on
for r = 1:3, plot(t, aa(r,:), "DisplayName", labels{r}); end
ylabel("\alpha (rad)"); xlabel("Time (s)"); legend show

% Figure 4: animation (slow motion; the real time window is T_total seconds)
figure(4); clf
xi = [-b - a, b - a];            % leading and trailing edge positions relative to elastic axis
Hd = gobjects(1,3);
for r = 1:3
    subplot(1,3,r); hold on; grid on; axis equal
    xlim([-0.8 0.8]*c); ylim([-0.5 0.5]*c)
    plot(0, 0, 'k+', 'MarkerSize', 10)   % elastic axis at rest
    yline(0, '--')
    Hd(r) = plot(nan, nan, 'LineWidth', 5, 'Color', [0 0.45 0.74]);
    title(labels{r})
end
for n = 1:n_frames
    for r = 1:3
        al = aa(r,n);
        xs = xi * cos(al);
        ys = -hh(r,n) - xi * sin(al);    % h positive downward, alpha positive nose-up
        set(Hd(r), 'XData', xs, 'YData', ys);
    end
    sgtitle(sprintf("U = %.1f m/s,  t = %.3f s", U_vis, t(n)))
    drawnow
end














%% ===== Local functions =====
function [lam, v] = pk_mode(omega, U, ii, p, tol)
    % p-k iteration for one mode; returns converged eigenvalue and eigenvector
    err = 10; iter = 0;
    while err > tol && iter < 200
        iter = iter + 1;
        A_k = build_Ak(omega, U, p);
        [V, D] = eig(A_k);
        lam_all = diag(D);
        [im_sorted, idx] = sort(imag(lam_all), 'ascend');
        err = abs(im_sorted(ii) - omega);
        omega = im_sorted(ii);
    end
    lam = lam_all(idx(ii));
    v   = V(:, idx(ii));
end

function A_k = build_Ak(omega, U, p)
    % Same matrices as your main loop
    m = p.m; S = p.S; I_alpha = p.I_alpha; a = p.a; b = p.b; rho = p.rho;
    k   = (omega * b) / U;
    C_k = 1 - (0.165 ./ (1 - 1i * (0.0455 ./ k))) - (0.335 ./ (1 - 1i * (0.30 ./ k)));

    M_matrix = [m, S; S, I_alpha];
    K_matrix = [p.K_h, 0; 0, p.K_alpha];
    D_matrix = rho * pi * (b ^ 2) * [-1, a; a, -((a ^ 2) + ((b ^ 2) / 8))];
    E_k_matrix = rho * pi * b * U * [-2 * C_k, -b + (2 * C_k * (a - (b / 2))); ...
                 b - (b - (2 * a + b) * C_k), (b - (2 * a + b) * C_k) * (a - (b / 2))];
    F_k_matrix = rho * pi * b * (U ^ 2) * [0, -2 * C_k; 0, b - (b - (2 * a + b) * C_k)];

    A_k = [zeros(2), eye(2); (M_matrix - D_matrix) \ (F_k_matrix - K_matrix), ...
                             (M_matrix - D_matrix) \ E_k_matrix];
end

function [U_f, f_f] = flutter_point(a, P)
% p-k flutter speed for one value of a. Returns NaN if no flutter in U_range.
b = P.c/2;
m = P.mu*P.c*P.t;
S = -m*a;
K_h = m*P.omega_h^2;
I_al = m*P.c^2/12 + m*a^2;
K_al = I_al*P.omega_alpha^2;
rho = P.rho;

M_matrix = [m, S; S, I_al];
K_matrix = [K_h, 0; 0, K_al];
D_matrix = rho*pi*b^2*[-1, a; a, -(a^2 + b^2/8)];
MD = M_matrix - D_matrix;

U_range = P.U_range;
nU = numel(U_range);
zeta = nan(nU,2);
wabs = nan(nU,2);
om = [P.omega_h, P.omega_alpha];     % start guesses, then tracked from the previous U

for i = 1:nU
    U = U_range(i);
    for j = 1:2
        lam = NaN;
        om_j = om(j);
        for it = 1:100
            k = om_j*b/U;
            C_k = 1 - 0.165/(1 - 1i*(0.0455/k)) - 0.335/(1 - 1i*(0.30/k));
            E = rho*pi*b*U*[-2*C_k, -b + 2*C_k*(a - b/2); ...
                b - (b - (2*a+b)*C_k), (b - (2*a+b)*C_k)*(a - b/2)];
            F = rho*pi*b*U^2*[0, -2*C_k; 0, b - (b - (2*a+b)*C_k)];
            A = [zeros(2), eye(2); MD\(F - K_matrix), MD\E];

            ev = eig(A);
            ev = ev(imag(ev) > 1e-9);          % oscillatory roots only (skips real/divergent roots)
            if isempty(ev), lam = NaN; break, end
            [~, idx] = min(abs(imag(ev) - om_j));   % track the root closest to the assumed frequency
            lam = ev(idx);
            err = abs(imag(lam) - om_j);
            om_j = imag(lam);
            if err < P.tol, break, end
        end
        if ~isnan(lam)
            zeta(i,j) = -real(lam)/abs(lam);
            wabs(i,j) = abs(lam);
            om(j) = om_j;
        end
    end
end

% First downward crossing of zero damping (threshold avoids round-off triggers)
U_f = NaN; f_f = NaN;
for j = 1:2
    idx = find(zeta(:,j) < -1e-6, 1);
    if isempty(idx), continue, end
    if idx > 1 && ~isnan(zeta(idx-1,j))
        z0 = zeta(idx-1,j); z1 = zeta(idx,j);
        r = -z0/(z1 - z0);                       % linear interpolation to zeta = 0
        Uc = U_range(idx-1) + r*(U_range(idx) - U_range(idx-1));
        wc = wabs(idx-1,j) + r*(wabs(idx,j) - wabs(idx-1,j));
    else
        Uc = U_range(idx); wc = wabs(idx,j);
    end
    if isnan(U_f) || Uc < U_f
        U_f = Uc; f_f = wc/2/pi;
    end
end
end