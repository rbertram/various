% Fig8.m

% Used to make Figure 8B from "Oscillations in Delayed Positive Feedback Systems", by
% C. Ryzowicz, R. Bertram, and B. Karamched, Physical Chemistry Chemical Physics,
% 26:24861, 2024.

% This code solves simple 1D DDE systems modelling delayed positive autoregulation
% Time series is plotted.

% X -> X activation of X by X
% DDE Equation:
% dx/dt = alpha*x^n / (1 + x^n) - x;

% Author: Christopher Ryzowicz

%---------------------------------------------------------------------- ---

clear all
clc
%clf

%%% Timespan and delay
Tstart = 0; 
Tend = 1200; 

tau = 14; % Delay value

%%% Parameters
p.alpha = 2.15;
p.n = -2;

%%% Solving DDE
opts = ddeset('RelTol',1e-6,'AbsTol',1e-6);
sol = dde23(@ddefunc, tau, @xhist, [Tstart Tend],opts,p);

%%% Plot y(t) versus t
figure(1),
plot(sol.x,sol.y)
xlabel('time')
ylabel('x')
box off
axis([0 1200 -1 1.5])
set(gca,'FontSize',20)

%------------------------------------------------------------------------

function dx = ddefunc(t, x, xl, p)
%%% Differential equation function

% Differential Equation
dx = p.alpha / (1 + (xl)^p.n) - x;

end

%------------------------------------------------------------------------

function x = xhist(t, p)
%%% History function
x = -0.75;
end
