% Fig2.m

% Used to make Figure 2 from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.


% This code solves a simple DDE modelling negative autoregulation of a single
% molecule. After, we plot both time series Pre and Post-Hopf bifurcation to demenstrate the different dynamics at steady state. 
% In addition, we plot bifurcation diagram as delay value changes
% tau. This code is set up for specific parameter values alphaX=10 and n=2.
% **Changing parameter values will not make bifurcation diagram accurate!!

% X -| X inhibition of X on itself
% Delay Equation:
% dx/dt = alphaX / (1 + xlag^n) - x;

% Author: Christopher Ryzowicz
%-------------------------------------------------------------------------
clear all
clc
% clf

%%% Parameters
p.alphaX = 10;
p.n = 2;

%% Figure 2B ------------------------------------------------------------------ 
%%% Timespan and delay
tau1 = [1 3];% Delay value

for i = 1:length(tau1)

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol1 = dde23(@ddefunc, tau1(i), @yhist, [0 80],opts,p);

    %%% Plot individual timeseries
    figure(1),
    plot(sol1.x,sol1.y(1,:),'LineWidth',2), hold on
    box off
    xlabel('t')
    ylabel('x')
end
set(gca,'FontSize',20);

%% Figure 2C ----------------------------------------------------------------------
%%% Timespan and delay
tau2 = 1:0.01:5;% Delay value
Tstart = 0; % Start time of solution
Tend = 200; % End time of solution

for i = 1:length(tau2)

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol2 = dde23(@ddefunc, tau2(i), @yhist, [Tstart Tend] ,opts,p);

    %%% Capture max and min of oscillation, needed to know where
    %%% oscillations started in terms of tau2 value!!
    if tau2(i)>=1.68
        max_points(i) = max(sol2.y(1,end-500:end));
        min_points(i) = min(sol2.y(1,end-500:end));
    end

end

%%% Plot steady state x vs. tau2 bifurcation diagram. I needed to know in
%%% advance where the hopf bifurcation was to get this to be accurate
figure(2),
plot(tau2(find(tau2>=1.69)),max_points(find(tau2>=1.69)),'r','LineWidth',2), hold on
plot(tau2(find(tau2>=1.69)),min_points(find(tau2>=1.69)),'r','LineWidth',2), hold on
xlabel('\tau')
ylabel('x')
plot(tau2(find(tau2<1.69)),2*ones(1,size(tau2(find(tau2<1.69)),2)),'k','LineWidth',2), hold on
plot(tau2(find(tau2>=1.69)),2*ones(1,size(tau2(find(tau2>=1.69)),2)),'k--','LineWidth',2), hold on
plot(1.7,2,'o','MarkerSize',12,'MarkerFaceColor','g'), hold on
yticks([0 2 4 6 8 10])
xticks([1 2 3 4 5])
set(gca,'FontSize',20);

%%% Save data into a .mat file
% save('Figure2B_BifurationDiagram_1DNegativeAutoregulation_tau1to5.mat')

%------------------------------------------------------------------------------

function dout = ddefunc(t, out, yl, p)
%%% Define the variables
x = out(1);
xlag = yl(1);

%%% Differential Equation
dx = p.alphaX / (1 + xlag^p.n) - x;

dout = dx;
end

%------------------------------------------------------------------------

function x_0 = yhist(t, p)
%%% History function
x_0 = 3;
end
