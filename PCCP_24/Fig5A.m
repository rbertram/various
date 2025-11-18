% Fig5A.m

% Used to make Figure 5A from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% This code simulates the positive-positive DDE models and calculates how long the oscillations
% last as a function of the delay, tau. The method of calculating how long these oscillations last is given in the text of the
% manuscript.

% Parameters_PosPos and SteadyState_Solutions_PosPos.m must be in current directory!!

% System of equations for positive-positive DDE model:
% dx/dt = (alphaX*ylag^n) / (1 + ylag^n) - x;
% dy/dt = (alphaY*xlag^n) / (1 + xlag^n) - y;

% Author: Christopher Ryzowicz
%-------------------------------------------------------------------------
clear all
clc
%  clf
%-------------------------------------------------------------------------

% Solve for steady state values for positive-positive model
[ss_x2, ss_y2] = SteadyState_Solutions_PosPos();
% Caculate midpoint of bistable steady state values
ss_midpoint2 = double((max(ss_x2)+min(ss_x2))/2);

%-------------------------------------------------------------------------

%%% Initialize list of delay values
tau_list = 2:0.5:8;

for i=1:length(tau_list)
i
    %%% Timespan and delay
    tau = tau_list(i); % Delay values
    Tstart = 0; % Start time of solution
    Tend = 4500; % End time of solution

    %%% Parameters
    Parameters_PosPos
    
    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol2 = dde23(@ddefunc2, tau, @yhist, [Tstart Tend],opts,q); % Solution for positive-positive model

    %%% Assign beta weight value for convex combination depending on if x reaches high or low steady state
    %%% Calculate time required for last 'x' value of oscillation to reach convex combination value with weight beta
    if sol2.y(1,end) < ss_midpoint2
        beta = 0.1;
        convex_combo = (1-beta)*min(ss_x2) + beta*max(ss_x2);
        oscillate_timeX2(i) = sol2.x(find(sol2.y(1,:)>convex_combo,1,'last'));
    else
        beta = 0.9;
        convex_combo = (1-beta)*min(ss_x2) + beta*max(ss_x2);
        oscillate_timeX2(i) = sol2.x(find(sol2.y(1,:)<convex_combo,1,'last'));
    end
end

%------------------------------------------------------------------------

%%% Plot scatter of oscillation duration vs. tau value
figure(2),scatter(tau_list,oscillate_timeX2./(2*tau_list),'filled','SizeData',100)
xlabel('\tau')
ylabel('Number of Oscillations')
box off
set(gca,'FontSize',24)

%%% Save data to a .mat file
% save('CodeXX9_Data.mat')

%----------------------------------------------------------------------

function dout = ddefunc2(t, out, yl, q)

%%% Define the variables
x = out(1);
y = out(2);
xlag = yl(1,:);
ylag = yl(2,:);

%%% Differential Equation
dx = (q.alphaX*ylag^q.n) / (1 + ylag^q.n) - x;
dy = (q.alphaY*xlag^q.n) / (1 + xlag^q.n) - y;

dout = [dx dy]';
end

%------------------------------------------------------------------------

function y = yhist(t,q)
%%% History function

x_0 = 1.2;
y_0 = 0.5; 

y = [x_0 y_0];
end
