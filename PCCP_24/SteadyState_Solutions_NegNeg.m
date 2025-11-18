% SteadyState_Solutions_NegNeg.m

function [solx, soly] = SteadyState_Solutions_NegNeg()

syms x y
% import parameters from Parameters_NegNeg.m file
Parameters_NegNeg

[solx,soly] = solve(p.alphaX/(1+y^p.n)-x==0, p.alphaY/(1+x^p.n)-y==0);

end
