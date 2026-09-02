%% Mistletoe_PRCC_v2.m
% Global sensitivity analysis of the tree-mistletoe-bird model
% LHS + regime classification + PRCC conditioned on the oscillatory regime
%
% Dimensionless model:
% x' = r*x*(1-x) - h*y*x
% y' = gamma*x*y + (m2*y*z)/(v+y) - d*y
% z' = z*(1-z) + (m3*y*z)/(v+y) + m4*x*z
%
% Methodological approach:
% 1) Sample parameters using LHS.
% 2) Integrate each realization for a sufficiently long time.
% 3) Classify the final solution as:
%       1 = equilibrium
%       2 = cycle/persistent oscillation
%       3 = extinction of y
%       0 = other/transient/unclassified
% 4) Compute amplitude, mean, and period ONLY for Reg == 2.
% 5) Compute PRCC for the oscillatory characteristics using only those samples.
%
% NOTE:
% Wide sampling (for example +/-30%) may cross bifurcations.
% Therefore, the PRCC obtained here is interpreted as sensitivity
% CONDITIONED on persistence of the oscillatory regime.

clear; clc; close all;

%% -------------------- Configuration --------------------

% Baseline parameters
base.r     = 0.25;
base.h     = 1.0;
base.gamma = 0.1;
base.m2    = 0.5;
base.m3    = 3.0;
base.m4    = 0.1;
base.v     = 1.0;
base.d     = 0.6555;   % fixed, near the Hopf point and with a stable cycle

% Parameters to be sampled
paramNames  = {'r','h','gamma','m2','m3','m4','v'};
paramLabels = {'r','h','\gamma','m_2','m_3','m_4','v'};
k = numel(paramNames);

% Sampling range relative to the baseline value
rangeFrac = 0.0005;

bounds = zeros(k,2);
for j = 1:k
    p0 = base.(paramNames{j});
    bounds(j,:) = [(1-rangeFrac)*p0, (1+rangeFrac)*p0];
end

% Latin Hypercube Sampling
N = 3500;
rng(1);
U = lhsdesign(N,k,'criterion','maximin','iterations',50);
P = scale_unit_to_bounds(U,bounds);

% Integration
T  = 5000;
dt = 0.5;
t_eval = 0:dt:T;

% The last half of the simulation is used to diagnose the regime.
% With T=2000, the final window contains 1000 time units,
% which is sufficient to observe several periods of the nominal cycle (~153).
transientFrac = 0.50;

% Initial condition
u0 = [0.9; 0.1; 2.1];

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Classification tolerances
tol.eq_abs       = 1e-4;   % small absolute amplitude -> equilibrium
tol.extinction   = 1e-6;   % y ~ 0
tol.amp_ratio_lo = 0.80;   % amplitude persistence between halves
tol.amp_ratio_hi = 1.20;
tol.period_cv    = 0.10;   % maximum relative variation among periods
tol.min_peaks    = 5;      % at least 5 peaks -> 4 measurable periods

%% -------------------- Preallocation --------------------

Ay    = nan(N,1);   % amplitude of y, cycles only
Ymean = nan(N,1);   % mean of y, cycles only
Per   = nan(N,1);   % period, cycles only

% Auxiliary diagnostics
Ay_all    = nan(N,1);
Ymean_all = nan(N,1);
Per_raw   = nan(N,1);

% 0=other/transient, 1=equilibrium, 2=cycle, 3=extinction of y
Reg = nan(N,1);

AmpRatio = nan(N,1);
PeriodCV = nan(N,1);
Npeaks   = zeros(N,1);

fprintf('Running %d simulations...\n',N);

%% -------------------- Simulations --------------------

for i = 1:N

    % Parameters for sample i
    r     = P(i,1);
    h     = P(i,2);
    gamma = P(i,3);
    m2    = P(i,4);
    m3    = P(i,5);
    m4    = P(i,6);
    v     = P(i,7);
    d     = base.d;

    % Integration
    sol = ode45(@(t,u) rhs_model(t,u,r,h,gamma,m2,m3,m4,v,d), ...
                [t_eval(1) t_eval(end)],u0,opts);

    Usol = deval(sol,t_eval);

    % Post-transient window
    cut = max(2,floor(transientFrac*numel(t_eval)));
    t2  = t_eval(cut:end);
    U2  = Usol(:,cut:end);

    x2 = U2(1,:);
    y2 = U2(2,:);
    z2 = U2(3,:);

    % Unfiltered outcomes, useful for diagnostics
    Ay_all(i)    = max(y2)-min(y2);
    Ymean_all(i) = mean(y2);

    % Robust classification using all three variables
    [Reg(i),diagOut] = classify_regime(t2,x2,y2,z2,tol);

    AmpRatio(i) = diagOut.amp_ratio;
    PeriodCV(i) = diagOut.period_cv;
    Npeaks(i)   = diagOut.n_peaks;
    Per_raw(i)  = diagOut.period;

    % Store outcomes ONLY if the solution is persistently oscillatory
    if Reg(i) == 2
        Ay(i)    = Ay_all(i);
        Ymean(i) = Ymean_all(i);
        Per(i)   = diagOut.period;
    end

    if mod(i,100)==0
        fprintf('  %d / %d\n',i,N);
    end
end

fprintf('Simulations completed.\n\n');

%% -------------------- Regime summary --------------------

Nother = sum(Reg==0);
Neq    = sum(Reg==1);
Nosc   = sum(Reg==2);
Next   = sum(Reg==3);

fprintf('Regime summary:\n');
fprintf('  Equilibrium              : %4d  (%6.2f%%)\n',Neq,  100*Neq/N);
fprintf('  Cycle/oscillation        : %4d  (%6.2f%%)\n',Nosc, 100*Nosc/N);
fprintf('  Extinction of y          : %4d  (%6.2f%%)\n',Next, 100*Next/N);
fprintf('  Other/transient          : %4d  (%6.2f%%)\n\n',Nother,100*Nother/N);

% Frequency figure
figure;
counts = [Neq Nosc Next Nother];
bar(counts);
set(gca,'XTick',1:4, ...
    'XTickLabel',{'Equilibrium','Cycle','Extinction y\to0','Other'}, ...
    'FontSize',12);
ylabel('Number of simulations');
title('Frequency of dynamical regimes');
grid on;
exportgraphics(gcf,'regimes_bar_v2.pdf','ContentType','vector');

%% -------------------- PRCC: oscillatory regime ONLY --------------------

osc = (Reg==2) & isfinite(Ay) & isfinite(Ymean) & isfinite(Per);

P_osc   = P(osc,:);
Ay_osc  = Ay(osc);
Ym_osc  = Ymean(osc);
Per_osc = Per(osc);

Nosc_valid = size(P_osc,1);

fprintf('Valid oscillatory samples for PRCC: %d of %d (%.2f%%)\n', ...
    Nosc_valid,N,100*Nosc_valid/N);

if Nosc_valid <= k+10
    error(['There are too few oscillatory samples to compute PRCC reliably. ' ...
           'Consider reducing the sampling range or increasing N.']);
elseif Nosc_valid < 0.30*N
    warning(['Less than 30%% of the samples remain in the oscillatory regime. ' ...
             'The PRCC should be interpreted with great caution; consider a second LHS ' ...
             'with a narrower range around the nominal parameter set.']);
end

prcc_Ay  = zeros(k,1); p_Ay  = zeros(k,1);
prcc_Ym  = zeros(k,1); p_Ym  = zeros(k,1);
prcc_Per = zeros(k,1); p_Per = zeros(k,1);

for j = 1:k

    idx = setdiff(1:k,j);

    [rho,pval] = partialcorr(P_osc(:,j),Ay_osc,P_osc(:,idx), ...
        "Type","Spearman","Rows","complete");
    prcc_Ay(j) = rho;
    p_Ay(j) = pval;

    [rho,pval] = partialcorr(P_osc(:,j),Ym_osc,P_osc(:,idx), ...
        "Type","Spearman","Rows","complete");
    prcc_Ym(j) = rho;
    p_Ym(j) = pval;

    [rho,pval] = partialcorr(P_osc(:,j),Per_osc,P_osc(:,idx), ...
        "Type","Spearman","Rows","complete");
    prcc_Per(j) = rho;
    p_Per(j) = pval;
end

%% -------------------- Tables --------------------

T_Ay = table(paramLabels',prcc_Ay,p_Ay, ...
    'VariableNames',{'Parameter','PRCC','p_value'});

T_Ym = table(paramLabels',prcc_Ym,p_Ym, ...
    'VariableNames',{'Parameter','PRCC','p_value'});

T_Per = table(paramLabels',prcc_Per,p_Per, ...
    'VariableNames',{'Parameter','PRCC','p_value'});

disp(' ');
disp('PRCC conditioned on the oscillatory regime');
disp('------------------------------------------');

disp('PRCC for amplitude of y:');
disp(T_Ay);

disp('PRCC for mean of y:');
disp(T_Ym);

disp('PRCC for period:');
disp(T_Per);

%% -------------------- PRCC figures --------------------

make_prcc_bar(prcc_Ay,paramLabels, ...
    'PRCC: Amplitude of y (oscillatory regime)', ...
    'prcc_Ay_v2.pdf');

make_prcc_bar(prcc_Ym,paramLabels, ...
    'PRCC: Mean of y (oscillatory regime)', ...
    'prcc_Ymean_v2.pdf');

make_prcc_bar(prcc_Per,paramLabels, ...
    'PRCC: Period (oscillatory regime)', ...
    'prcc_Period_v2.pdf');

%% -------------------- Monotonicity diagnostics --------------------
% PRCC is appropriate when the parameter-response relationship is
% approximately monotonic. These figures allow this assumption to be inspected.

for j = 1:k
    figure;
    scatter(P_osc(:,j),Per_osc,18,'filled');
    xlabel(paramLabels{j},'Interpreter','tex');
    ylabel('Period');
    title(['Oscillatory samples: ',paramLabels{j},' vs period'], ...
        'Interpreter','tex');
    grid on;
    exportgraphics(gcf,sprintf('scatter_%s_period.pdf',paramNames{j}), ...
        'ContentType','vector');
end

%% -------------------- Save results --------------------

Results = table((1:N)',Reg,Ay_all,Ymean_all,Per_raw,AmpRatio,PeriodCV,Npeaks, ...
    'VariableNames',{'Simulation','Regime','AmplitudeY','MeanY', ...
    'Period','AmplitudeRatio','PeriodCV','Npeaks'});

ParamTable = array2table(P,'VariableNames',paramNames);

writetable([ParamTable Results],'LHS_regime_results_v2.csv');

save('Mistletoe_PRCC_v2_results.mat', ...
    'base','bounds','P','Reg','Ay','Ymean','Per', ...
    'Ay_all','Ymean_all','Per_raw','AmpRatio','PeriodCV','Npeaks', ...
    'prcc_Ay','prcc_Ym','prcc_Per','p_Ay','p_Ym','p_Per', ...
    'paramNames','paramLabels');

fprintf('\nResults and figures exported successfully.\n');

%% ==================== LOCAL FUNCTIONS ====================

function P = scale_unit_to_bounds(U,bounds)
    [N,k] = size(U);
    P = zeros(N,k);
    for j = 1:k
        lo = bounds(j,1);
        hi = bounds(j,2);
        P(:,j) = lo + (hi-lo).*U(:,j);
    end
end

function dudt = rhs_model(~,u,r,h,gamma,m2,m3,m4,v,d)
    x = u(1);
    y = u(2);
    z = u(3);

    dx = r*x*(1-x) - h*y*x;
    dy = gamma*x*y + (m2*y*z)/(v+y) - d*y;
    dz = z*(1-z) + (m3*y*z)/(v+y) + m4*x*z;

    dudt = [dx;dy;dz];
end

function [reg,out] = classify_regime(t,x,y,z,tol)
% Classification:
% 3 = extinction of y
% 1 = equilibrium
% 2 = cycle/persistent oscillation
% 0 = other/transient/unclassified

    n = numel(t);
    mid = floor(n/2);

    idx1 = 1:mid;
    idx2 = (mid+1):n;

    % Amplitudes in the two halves of the final window
    amp1_y = max(y(idx1))-min(y(idx1));
    amp2_y = max(y(idx2))-min(y(idx2));

    amp2_x = max(x(idx2))-min(x(idx2));
    amp2_z = max(z(idx2))-min(z(idx2));

    amp_ratio = amp2_y/max(amp1_y,eps);

    % Final level of y
    nTail = min(50,n);
    y_end = mean(y(end-nTail+1:end));

    % Period and peak regularity over the entire final window
    [Tper,period_cv,n_peaks] = estimate_period_stable(t,y,tol.min_peaks);

    % Diagnostic outputs
    out.period    = Tper;
    out.period_cv = period_cv;
    out.n_peaks   = n_peaks;
    out.amp_ratio = amp_ratio;

    % 1) Extinction
    if y_end < tol.extinction && max(y(idx2)) < 10*tol.extinction
        reg = 3;
        return
    end

    % 2) Equilibrium: all three variables have very small amplitude
    if max([amp2_x amp2_y amp2_z]) < tol.eq_abs
        reg = 1;
        return
    end

    % 3) Persistent cycle:
    %    - amplitude of y does not collapse between halves
    %    - enough peaks
    %    - approximately constant periods
    persistent_amp = amp_ratio >= tol.amp_ratio_lo && ...
                     amp_ratio <= tol.amp_ratio_hi;

    regular_period = isfinite(Tper) && ...
                     n_peaks >= tol.min_peaks && ...
                     period_cv <= tol.period_cv;

    if persistent_amp && regular_period
        reg = 2;
    else
        reg = 0;
    end
end

function [Tper,cv_period,n_peaks] = estimate_period_stable(t,y,min_peaks)
% Period estimation from local maxima.
% A minimum number of peaks and consistency among periods are required.

    dy = diff(y);
    peaks = find(dy(1:end-1)>0 & dy(2:end)<=0)+1;

    n_peaks = numel(peaks);

    if n_peaks >= min_peaks
        periods = diff(t(peaks));
        Tper = median(periods);

        if mean(periods)>0
            cv_period = std(periods)/mean(periods);
        else
            cv_period = Inf;
        end
    else
        Tper = NaN;
        cv_period = NaN;
    end
end

function make_prcc_bar(prcc_vals,names,ttl,outpdf)
    figure;
    b = bar(prcc_vals);
    b.FaceColor = [0.75 0.88 0.75];
    b.EdgeColor = 'k';

    set(gca,'XTick',1:numel(names), ...
        'XTickLabel',names, ...
        'TickLabelInterpreter','tex', ...
        'FontSize',12);

    ylabel('PRCC','Interpreter','tex');
    xlabel('Parameters');
    title(ttl,'Interpreter','tex');
    grid on;
    yline(0,'-');
    ylim([-1.1 1.1]);

    x = 1:numel(prcc_vals);

    for i = 1:numel(prcc_vals)
        yy = prcc_vals(i);

        if yy >= 0
            text(x(i),yy+0.03,sprintf('%.2f',yy), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',11);
        else
            text(x(i),yy-0.03,sprintf('%.2f',yy), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','top', ...
                'FontSize',11);
        end
    end

    exportgraphics(gcf,outpdf,'ContentType','vector');
end
