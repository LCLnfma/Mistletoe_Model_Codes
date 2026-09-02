%%%%%%%%%%%%%%%%%%%%%
% Monodromy matrix for the three-equation system
%%%%%%%%%%%%%%%%%%%%%
clc;
clearvars;
format long;

% Plot text parameters
x_fs               = 28;                % x-axis label font size  
y_fs               = 28;                % y-axis label font size
z_fs               = 28;
xt_fs              = 18;                % x-axis tick font size
yt_fs              = xt_fs;             % y-axis tick font size
zt_fs              = xt_fs;
t_fs               = 12;                % title font size
l_fs               = 16;                % legend font size
l_width            = 2.0;               % line width
grey               = [0.5 0.5 0.5];     % gray color

% Model parameters
 r      = 0.25;
 h      = 1;
 c      = 0.1;
 m2     = 0.5;
 k      = 1;
 %d      = 0.655;
 d      = 0.6555;
 m3     = 3;
 m4     = 0.1; 
 
init_cond          = [0.15556 , 0.2267, 1.5701];      % initial conditions
time_num           = 4000;              % number of intervals used to divide the period
pts_per_interval   = 20;                % number of integration points per interval
num_periods        = 1;                 % number of periods
T                  = 143.16681491;       % period
T0                 = 0 * T;             % initial time
T1                 = num_periods * T;   % total integration time; it may be equal to or greater than one period
T_threshold        = 1 * T;             % time used to locate X_end, the monodromy matrix
steps              = 10;  

% Parameters used for derivative calculations and their approximations
a_arr              = linspace(-1, 1, steps);
b_arr              = linspace(-1, 1, steps);
time_arr           = [0 T1];
time_lin           = linspace(time_arr(1), time_arr(2), time_num); % ms
time_arr_inv       = [0 T];
time_lin_inv       = linspace(time_arr_inv(1), time_arr_inv(2), time_num); % ms
dt                 = time_lin(2) - time_lin(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Start..\n');
tic();

% x(1) = v %
% x(2) = n % 
% Integration of the nonlinear model; insert your model equations here
f = @(t,x) [r*x(1)*(1 - x(1)) - h*x(1)* x(2); 
            c*x(1)* x(2) + (m2*x(2)*x(3))/(k+x(2)) - d*x(2);
            x(3)*(1 - x(3)) + (m3*x(2)*x(3))/(k+x(2)) + m4*x(1)*x(3)];

        [t,xa] = ode45(f, time_lin, init_cond);

% Computation of the Jacobian matrix
syms x y z

% Insert the first model equation in f = 
f       = r*x*(1 - x) - h*x*y; 
a11     = diff(f, x);
a12     = diff(f, y);
a13     = diff(f, z);
f_a11   = symfun(a11, [x y z]);
f_a12   = symfun(a12, [x y z]);
f_a13   = symfun(a13, [x y z]);
num_a11 = eval(f_a11(xa(:,1), xa(:,2), xa(:,3)));
num_a12 = eval(f_a12(xa(:,1), xa(:,2), xa(:,3)));
num_a13 = eval(f_a13(xa(:,1), xa(:,2), xa(:,3)));
% Insert the second model equation in g = 
g       = c*x*y + (m2*y*z)/(k+y) - d*y;
a21     = diff(g, x);
a22     = diff(g, y);
a23     = diff(g, z);
f_a21   = symfun(a21, [x y z]);
f_a22   = symfun(a22, [x y z]);
f_a23   = symfun(a23, [x y z]);
num_a21 = eval(f_a21(xa(:,1), xa(:,2), xa(:,3)));
num_a22 = eval(f_a22(xa(:,1), xa(:,2), xa(:,3)));
num_a23 = eval(f_a23(xa(:,1), xa(:,2), xa(:,3)));
% Insert the third model equation in hh = 
hh       = z*(1 - z) + (m3*y*z)/(k+y) + m4*x*z;
a31     = diff(hh, x);
a32     = diff(hh, y);
a33     = diff(hh, z);
f_a31   = symfun(a31, [x y z]);
f_a32   = symfun(a32, [x y z]);
f_a33   = symfun(a33, [x y z]);
num_a31 = eval(f_a31(xa(:,1), xa(:,2), xa(:,3)));
num_a32 = eval(f_a32(xa(:,1), xa(:,2), xa(:,3)));
num_a33 = eval(f_a33(xa(:,1), xa(:,2), xa(:,3)));
% Spline approximation
% From this point onward, no changes are needed except for figure formatting
c11     = spline(time_lin, num_a11);
c12     = spline(time_lin, num_a12);
c13     = spline(time_lin, num_a13);
c21     = spline(time_lin, num_a21);
c22     = spline(time_lin, num_a22);
c23     = spline(time_lin, num_a23);
c31     = spline(time_lin, num_a31);
c32     = spline(time_lin, num_a32);
c33     = spline(time_lin, num_a33);


[br11, a11c, L, order, dim] = unmkpp(c11);
[br12, a12c, L, order, dim] = unmkpp(c12);
[br13, a13c, L, order, dim] = unmkpp(c13);
[br21, a21c, L, order, dim] = unmkpp(c21);
[br22, a22c, L, order, dim] = unmkpp(c22);
[br23, a23c, L, order, dim] = unmkpp(c23);
[br31, a31c, L, order, dim] = unmkpp(c31);
[br32, a32c, L, order, dim] = unmkpp(c32);
[br33, a33c, L, order, dim] = unmkpp(c33);

tzz     = [];
ca      = [];
tww     = [];
wa      = [];
ts      = [];
sr      = [];

ic_0_1  = [1 0 0];
time    = [];

for i=1:time_num-1
  f_i = @(ti, zz)[(a11c(i,1)*(ti - time_lin(i))^3 + a11c(i,2)*(ti - time_lin(i))^2 + a11c(i,3)*(ti - time_lin(i)) + a11c(i,4))*zz(1)+(a12c(i,1)*(ti - time_lin(i))^3 + a12c(i,2)*(ti - time_lin(i))^2 + a12c(i,3)*(ti - time_lin(i)) + a12c(i,4))*zz(2)+(a13c(i,1)*(ti - time_lin(i))^3 + a13c(i,2)*(ti - time_lin(i))^2 + a13c(i,3)*(ti - time_lin(i)) + a13c(i,4))*zz(3);
                  (a21c(i,1)*(ti - time_lin(i))^3 + a21c(i,2)*(ti - time_lin(i))^2 + a21c(i,3)*(ti - time_lin(i)) + a21c(i,4))*zz(1)+(a22c(i,1)*(ti - time_lin(i))^3 + a22c(i,2)*(ti - time_lin(i))^2 + a22c(i,3)*(ti - time_lin(i)) + a22c(i,4))*zz(2)+(a23c(i,1)*(ti - time_lin(i))^3 + a23c(i,2)*(ti - time_lin(i))^2 + a23c(i,3)*(ti - time_lin(i)) + a23c(i,4))*zz(3);
                  (a31c(i,1)*(ti - time_lin(i))^3 + a31c(i,2)*(ti - time_lin(i))^2 + a31c(i,3)*(ti - time_lin(i)) + a31c(i,4))*zz(1)+(a32c(i,1)*(ti - time_lin(i))^3 + a32c(i,2)*(ti - time_lin(i))^2 + a32c(i,3)*(ti - time_lin(i)) + a32c(i,4))*zz(2)+(a33c(i,1)*(ti - time_lin(i))^3 + a33c(i,2)*(ti - time_lin(i))^2 + a33c(i,3)*(ti - time_lin(i)) + a33c(i,4))*zz(3)];          
  time_int = linspace(time_lin(i), time_lin(i+1), pts_per_interval);
  [t2, za] = ode23(f_i, time_int, ic_0_1);
  tzz       = cat(1, tzz, t2);
  ca       = cat(1, ca, za);
  ic_0_1   = [ca(end,1) ca(end,2) ca(end,3)];
  time     = [time time_int];
end

ic_0_2  = [0 1 0];
for i=1:time_num-1
  f_i = @(ti, ww) [(a11c(i,1)*(ti - time_lin(i))^3 + a11c(i,2)*(ti - time_lin(i))^2 + a11c(i,3)*(ti - time_lin(i)) + a11c(i,4))*ww(1)+(a12c(i,1)*(ti - time_lin(i))^3 + a12c(i,2)*(ti - time_lin(i))^2 + a12c(i,3)*(ti - time_lin(i)) + a12c(i,4))*ww(2)+(a13c(i,1)*(ti - time_lin(i))^3 + a13c(i,2)*(ti - time_lin(i))^2 + a13c(i,3)*(ti - time_lin(i)) + a13c(i,4))*ww(3);
                  (a21c(i,1)*(ti - time_lin(i))^3 + a21c(i,2)*(ti - time_lin(i))^2 + a21c(i,3)*(ti - time_lin(i)) + a21c(i,4))*ww(1)+(a22c(i,1)*(ti - time_lin(i))^3 + a22c(i,2)*(ti - time_lin(i))^2 + a22c(i,3)*(ti - time_lin(i)) + a22c(i,4))*ww(2)+(a23c(i,1)*(ti - time_lin(i))^3 + a23c(i,2)*(ti - time_lin(i))^2 + a23c(i,3)*(ti - time_lin(i)) + a23c(i,4))*ww(3);
                  (a31c(i,1)*(ti - time_lin(i))^3 + a31c(i,2)*(ti - time_lin(i))^2 + a31c(i,3)*(ti - time_lin(i)) + a31c(i,4))*ww(1)+(a32c(i,1)*(ti - time_lin(i))^3 + a32c(i,2)*(ti - time_lin(i))^2 + a32c(i,3)*(ti - time_lin(i)) + a32c(i,4))*ww(2)+(a33c(i,1)*(ti - time_lin(i))^3 + a33c(i,2)*(ti - time_lin(i))^2 + a33c(i,3)*(ti - time_lin(i)) + a33c(i,4))*ww(3)];
  
  time_int = linspace(time_lin(i), time_lin(i+1), pts_per_interval);
  [t2, ba] = ode23(f_i, time_int, ic_0_2);
  tww       = cat(1, tww, t2);
  wa       = cat(1, wa, ba);
  ic_0_2   = [wa(end,1) wa(end,2) wa(end,3)];
end

ic_0_3  = [0 0 1];
for i=1:time_num-1
  f_i = @(ti, s) [(a11c(i,1)*(ti - time_lin(i))^3 + a11c(i,2)*(ti - time_lin(i))^2 + a11c(i,3)*(ti - time_lin(i)) + a11c(i,4))*s(1)+(a12c(i,1)*(ti - time_lin(i))^3 + a12c(i,2)*(ti - time_lin(i))^2 + a12c(i,3)*(ti - time_lin(i)) + a12c(i,4))*s(2)+(a13c(i,1)*(ti - time_lin(i))^3 + a13c(i,2)*(ti - time_lin(i))^2 + a13c(i,3)*(ti - time_lin(i)) + a13c(i,4))*s(3);
                  (a21c(i,1)*(ti - time_lin(i))^3 + a21c(i,2)*(ti - time_lin(i))^2 + a21c(i,3)*(ti - time_lin(i)) + a21c(i,4))*s(1)+(a22c(i,1)*(ti - time_lin(i))^3 + a22c(i,2)*(ti - time_lin(i))^2 + a22c(i,3)*(ti - time_lin(i)) + a22c(i,4))*s(2)+(a23c(i,1)*(ti - time_lin(i))^3 + a23c(i,2)*(ti - time_lin(i))^2 + a23c(i,3)*(ti - time_lin(i)) + a23c(i,4))*s(3);
                  (a31c(i,1)*(ti - time_lin(i))^3 + a31c(i,2)*(ti - time_lin(i))^2 + a31c(i,3)*(ti - time_lin(i)) + a31c(i,4))*s(1)+(a32c(i,1)*(ti - time_lin(i))^3 + a32c(i,2)*(ti - time_lin(i))^2 + a32c(i,3)*(ti - time_lin(i)) + a32c(i,4))*s(2)+(a33c(i,1)*(ti - time_lin(i))^3 + a33c(i,2)*(ti - time_lin(i))^2 + a33c(i,3)*(ti - time_lin(i)) + a33c(i,4))*s(3)];
  
  time_int = linspace(time_lin(i), time_lin(i+1), pts_per_interval);
  [t2, sa] = ode23(f_i, time_int, ic_0_3);
  ts       = cat(1, ts, t2);
  sr       = cat(1, sr, sa);
  ic_0_3   = [sr(end,1) sr(end,2) sr(end,3)];
end


% Computation of X_end (monodromy matrix)
x_tmp = find(time > T_threshold);
if ( ~isempty(x_tmp) )
  X_ind_tmp = x_tmp(1);
  X_end     = [ca(X_ind_tmp,:); wa(X_ind_tmp,:); sr(X_ind_tmp,:)];
else
  X_end     = [ca(end,:); wa(end,:); sr(end,:)];
end
disp('Monodromy matrix X(T)='); disp(X_end);

[V, D] = eig(X_end); % X * V = V * D
disp('Eigenvalues (Floquet multipliers)'); disp(D);
disp('Eigenvectors');  disp(V);

figure(1);
set(gcf, 'Color', 'w');
plot3(xa(:,1),xa(:,2),xa(:,3),'k','Linewidth',3); hold on;
plot3(xa(1,1),xa(1,2),xa(1,3),'r.','MarkerSize', 30);
plot3(xa(end,1),xa(end,2),xa(end,3),'b.','MarkerSize', 30);
grid on;
xt = get(gca, 'XTick'); set(gca, 'FontSize', xt_fs);
yt = get(gca, 'YTick'); set(gca, 'FontSize', yt_fs);
zt = get(gca, 'YTick'); set(gca, 'FontSize', zt_fs);
xlabel('x-trees', 'Interpreter', 'latex', 'FontSize', x_fs);
ylabel('y-mistletoe', 'Interpreter', 'latex', 'FontSize', y_fs);
zlabel('z-birds', 'Interpreter', 'latex', 'FontSize', z_fs);

figure(2)
plot(t, xa(: ,1),'r.', t , xa(: ,3),'k', t , xa(: ,2),'b.')
% title('Time');
xlabel('Time t');
ylabel('Solution y');

% % figure();
% % set(gcf, 'Color', 'w');
% % plot3(xa(:,1),xa(:,3),xa(:,2),'k','Linewidth',3); hold on;
% % plot3(xa(1,1),xa(1,3),xa(1,2),'r.','MarkerSize', 30);
% % plot3(xa(end,1), xa(end,3),xa(end,2),'b.','MarkerSize', 30);
% % grid on;
% % xt = get(gca, 'XTick'); set(gca, 'FontSize', xt_fs);
% % yt = get(gca, 'YTick'); set(gca, 'FontSize', yt_fs);
% % zt = get(gca, 'YTick'); set(gca, 'FontSize', zt_fs);
% % xlabel('x-trees', 'Interpreter', 'latex', 'FontSize', x_fs);
% % ylabel('z-birds', 'Interpreter', 'latex', 'FontSize', y_fs);
% % zlabel('y-mistletoe', 'Interpreter', 'latex', 'FontSize', z_fs);
%ylim([0 1]);

