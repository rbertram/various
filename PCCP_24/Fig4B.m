% Fig4B.m

% Used to make Figure 4B  from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% This code simulates the negative-negative DDE and calculates how long the oscillations
% last as a function of the distance from the stable manifold (seperatrix). The method of calculating how long these oscillations last is given in the text of the
% manuscript.

% Parameters_NegNeg.m and SteadyState_Solutions_NegNeg.m must be in current directory!!

% System of equations for negative-negative DDE model:
% dx/dt = alphaX / (1 + ylag^n) - x;
% dy/dt = alphaY / (1 + xlag^n) - y;

% Author: Christopher Ryzowicz 
%-------------------------------------------------------------------------
clear all
clc
%  clf
%-------------------------------------------------------------------------

% Solve for steady state values for negative-negative model
[ss_x1, ss_y1] = SteadyState_Solutions_NegNeg();
% Caculate midpoint of bistable steady state values
ss_midpoint1 = double((max(ss_x1)+min(ss_x1))/2);

%-------------------------------------------------------------------------

% Initialize history function pairings
x0_list = 4.2:0.1:6.2;
y0_list = 8 - x0_list;
hist_list = {x0_list, y0_list};

% Calculate distance list to stable manifold using predetermined orthoganal line y=8-x and
% history values in reference to point (4,4)
dist_to_manifold = sqrt((x0_list-4).^2 + (y0_list-4).^2);

for i=1:length(hist_list{1})
i
    %%% Timespan and delay
    tau = 6; % Delay values
    Tstart = 0; % Start time of solution
    Tend = 2000; % End time of solution

    %%% Parameters
    Parameters_NegNeg

    %%% Solving DDE
    opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
    sol1 = dde23(@ddefunc1, tau, [hist_list{1}(i), hist_list{2}(i)], [Tstart Tend],opts,p); % Solution for negative-negative model

    %%% Assign beta weight value for convex combination depending on if x reaches high or low steady state
    %%% Calculate time required for last 'x' value of oscillation to reach
    %%% convex combination value with weight beta
    if sol1.y(1,end) < ss_midpoint1
        beta = 0.1;
        convex_combo = (1-beta)*min(ss_x1) + beta*max(ss_x1);
        oscillate_timeX1(i) = sol1.x(find(sol1.y(1,:)>convex_combo,1,'last'));
    else
        beta = 0.9;
        convex_combo = (1-beta)*min(ss_x1) + beta*max(ss_x1);
        oscillate_timeX1(i) = sol1.x(find(sol1.y(1,:)<convex_combo,1,'last'));
    end
end

%------------------------------------------------------------------------

%%% Plot scatter of oscillation duration vs. distance to stable manifold
figure(2),scatter(dist_to_manifold,oscillate_timeX1./(2*tau),'filled','SizeData',100)
xlabel('\Delta')
ylabel('Number of Oscillations')
box off
set(gca,'FontSize',20)

%%% Save data in a .mat file
% save('CodeXX5_Data.mat')

%----------------------------------------------------------------------

function dout = ddefunc1(t, out, yl, p)

%%% Define the variables
x = out(1);
y = out(2);
xlag = yl(1,:);
ylag = yl(2,:);

%%% Differential Equation
dx = p.alphaX / (1 + ylag^p.n) - x;
dy = p.alphaY / (1 + xlag^p.n) - y;

dout = [dx dy]';
end
