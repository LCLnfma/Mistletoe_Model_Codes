%% ============================================================
% Tree-mistletoe-bird model
% Numerical simulations and equilibrium stability
% Author: Laura Cruzado Lima
% ============================================================

clear
clc
close all

%% Model parameters

r     = 0.25;
h     = 1;
gamma = 0.1;
m2    = 0.5;
m3    = 3;
m4    = 0.1;
v     = 1;

%% %% Parameter d

case_d = "d2";     % Options: "d1", "d2", "d3"

switch case_d

    case "d1"
        d = 0.6582;      % Stable focus

    case "d2"
        d = 0.6555;      % Stable limit cycle

    case "d3"
        d = 0.6654;      % Stable node

    otherwise
        error('Unknown value. Select "d1", "d2", or "d3".')
end

fprintf('\nSelected case: %s\n',case_d);
fprintf('d = %.4f\n\n',d);
%% Initial conditions

IC = [
    0.5    0.3     0.3;
    0.23   1.507   0.1929;
    0      0.33    1.74;
    0.9    0.01    1.2;
    0.45   0.01    1.35;
    1      0.7     0.7;
    0.936  1.1405  0.0158;
    0.9    0.6     0.1;
    0.9    0.1     2.1;
    0.92   0.2     0.01;
    0.4    1.3     0.4;
    0.2    0.6     1;
    0      0.2     0.5
];
%% Numerical integration

Tsim = 10000;

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

nTraj = size(IC,1);

time = cell(nTraj,1);
sol  = cell(nTraj,1);

for j = 1:nTraj

    [time{j},sol{j}] = ode45( ...
        @(t,u) rhs_model(t,u,r,h,gamma,m2,m3,m4,v,d), ...
        [0 Tsim],IC(j,:),opts);

end
%% Equilibria and local stability

syms x y z

F = [
    r*x*(1-x) - h*x*y;
    gamma*x*y + (m2*y*z)/(v+y) - d*y;
    z*(1-z) + (m3*y*z)/(v+y) + m4*x*z
];

J = jacobian(F,[x y z]);

% Solve equilibrium equations
[solx,soly,solz] = vpasolve(F == 0,[x y z]);

C = double([solx soly solz]);

% Keep only real, nonnegative equilibria
tol = 1e-8;

C_pos = [];

for i = 1:size(C,1)

    Ci = C(i,:);

    if max(abs(imag(Ci))) < tol

        Ci = real(Ci);

        if all(Ci >= -tol)

            Ci(abs(Ci) < tol) = 0;

            C_pos = [C_pos; Ci];

        end
    end
end
%% Stability classification

nEq = size(C_pos,1);

stability = strings(nEq,1);
eigValues = cell(nEq,1);

tolEig = 1e-7;

fprintf('\nEquilibria and stability\n');
fprintf('====================================================\n');

for i = 1:nEq

    Ei = C_pos(i,:);

    Ji = double(subs(J,[x y z],Ei));

    lambda = eig(Ji);

    eigValues{i} = lambda;

    if all(real(lambda) < -tolEig)

        stability(i) = "Stable";

    elseif any(real(lambda) > tolEig)

        stability(i) = "Unstable";

    else

        stability(i) = "Nonhyperbolic";

    end

    fprintf('\nEquilibrium %d\n',i);
    fprintf('(x,y,z) = (%.6f, %.6f, %.6f)\n', ...
        Ei(1),Ei(2),Ei(3));

    fprintf('Eigenvalues:\n');

    for q = 1:3
        fprintf(' lambda_%d = % .8f %+.8fi\n', ...
            q,real(lambda(q)),imag(lambda(q)));
    end

    fprintf('Classification: %s\n',stability(i));

end
%% Three-dimensional phase portrait

figure
hold on

for j = 1:nTraj

    plot3(sol{j}(:,1), ...
          sol{j}(:,2), ...
          sol{j}(:,3), ...
          'b','LineWidth',1.2);

end

% Plot equilibria according to their stability

for i = 1:nEq

    Ei = C_pos(i,:);

    switch stability(i)

        case "Stable"

            plot3(Ei(1),Ei(2),Ei(3),'r*', ...
                'MarkerSize',20, ...
                'MarkerFaceColor','r');

        case "Unstable"

            plot3(Ei(1),Ei(2),Ei(3),'g*', ...
                'MarkerSize',20, ...
                'MarkerFaceColor','g');

        case "Nonhyperbolic"

            plot3(Ei(1),Ei(2),Ei(3),'ko', ...
                'MarkerSize',10, ...
                'MarkerFaceColor','k');

    end

end

xlabel('x - Trees')
ylabel('y - Mistletoes')
zlabel('z - Birds')

grid on
box on

switch case_d
    case "d1"
        d_label = 'd_1';
    case "d2"
        d_label = 'd_2';
    case "d3"
        d_label = 'd_3';
end

title(sprintf('$%s = %.4f$',d_label,d), ...
      'Interpreter','latex')
view(3)

%% ============================================================
% Local function
% ============================================================

function du = rhs_model(~,u,r,h,gamma,m2,m3,m4,v,d)

x = u(1);
y = u(2);
z = u(3);

dx = r*x*(1-x) - h*x*y;

dy = gamma*x*y ...
   + (m2*y*z)/(v+y) ...
   - d*y;

dz = z*(1-z) ...
   + (m3*y*z)/(v+y) ...
   + m4*x*z;

du = [dx;dy;dz];

end