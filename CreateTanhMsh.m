function [msh] = CreateTanhMsh(np, x0, L, delta, alpha)
%CREATETANHMSH Summary of this function goes here
%   Detailed explanation goes here

coordinates = "cylindrical";


gamma = -0.5*log(delta/4);%delta=2.4577e-05
tau = linspace(0, L, np+1)'; % x_int coordinates in a domain from 0 to L
eta = tanh(gamma*(2*tau/L -1+alpha))./tanh(gamma); % map in the x axis of tanh and then compute tanh
a = L/(eta(end)-eta(1));
b = -a*eta(1);
%x_int = L * (eta - eta(1)) / (eta(end) - eta(1)); % map back to the x domain
x_int = a*eta+b;
msh = CreateMsh(x_int, x0, coordinates);

end
