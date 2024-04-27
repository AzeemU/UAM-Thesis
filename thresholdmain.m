
close all; clear; clc

H = readmatrix('./network_data/Matrix_H.csv'); % community-route incidence
F = readmatrix('./network_data/Matrix_F.csv'); % link-route incidence
E = readmatrix('./network_data/Matrix_E.csv'); % link-node incidence
K = E>0; % unsingned link-node incidence
c = readmatrix('./network_data/Node_Capacities.csv'); % node capacity

d = readmatrix('./network_data/Edge_Capacities.csv'); % link capacity

[nn, nl] = size(E); % nn: # of nodes, nl: # of links
[~, nr] = size(F); % nr: # of routes
[nc, ~] = size(H); % nc: # of communities
% c = max(c)*ones(nn, 1);


C = [c, 0.8*c, 0.6*c]; % possible node cap
D = [d, 0.8*d, 0.6*d]; % possible link cap
p = [0.5; 0.3; 0.2]; % nominal distribution of capacity

N = size(C, 2);
del = 0.9; % risk sensitivity parameter
%alp = 5; % fairness parameter
delthresh = 1000;

viotol = 0.1; %tolerance of risk constraint violation
FWeps = 1e-5; % stopping criterion for FW
FWmaxiter = 100; % max iter num of FW
valhist = zeros(FWmaxiter+1, 1);

crmname = 'tv'; % choose the type of coherent risk measures

% compute the initial linearization point used in Frank-Wolfe (must be positive)
yalmip('clear')

x = sdpvar(nc, 1, 'full'); % community flow
y = sdpvar(nl, 1, 'full'); % link flow
z = sdpvar(nr, 1, 'full'); % route flow
w = sdpvar(N, 1, 'full'); % dual variable for risk measures
nu = sdpvar(1, 1, 'full'); % dual variable for risk measures
lam = sdpvar(1, 1, 'full'); % dual variable for risk measures
h = sdpvar(N, 1, 'full'); % cap violation
u = sdpvar(N, 1, 'full'); % auxiliary variable for expcone
tau = sdpvar(1, 1, 'full'); % auxiliary variable for expcone
xlb = sdpvar(1, 1, 'full');

%threshold fairness variables, Chen p614
v = sdpvar(nc, 1, 'full');
wthresh = sdpvar(1,1,'full');


constr = [E*y == zeros(nn, 1), ...
    x == H*z, ...
    F*z <= y, ...
    x >= xlb*ones(nc, 1), ...
    z >= zeros(nr, 1), ...
    y >= zeros(nl, 1), ...
    h >= zeros(N, 1), ...
    w - nu*ones(N, 1) >= h, ...
    lam >= 0,...
    wthresh >= 0, ... %threshold fairness constraints
    v <= x-delthresh, ...
    v <= ones(nc, 1).*wthresh, ...
    ones(nc, 1).*wthresh <= x];

for k = 1:N
    constr = [constr, ...
        K*y <= (1 + h(k))*C(:, k), ...
        y <= (1 + h(k))*D(:, k)];
end

if strcmp(crmname, 'cvar') % conditional VaR
    constr = [constr, -nu + (1/(1-del))*lam <= viotol, ...
        norm(diag(p)*w, 1) <= lam];
elseif strcmp(crmname, 'evar') % entropic VaR
    constr = [constr, -nu + tau -log(1-del)*lam  <= viotol, ...
        p'*u <= lam, ...
        expcone([w'-tau*ones(1, N); lam*ones(1, N); u'])];
elseif strcmp(crmname, 'tv') % total variation
    constr = [constr, -nu + 2*del*lam+p'*w <= viotol, ...
        norm(w, Inf) <= lam];
else
    error('Unknown type of risk measure!')
end

options = sdpsettings('verbose',1,'solver','mosek');

objectivefunc = nc*delthresh + sum(v);

%solution = optimize(constr, -xlb, options);
solution = optimize(constr, -objectivefunc, options);

%%
xopt = value(x); % optimal community flow
yopt = value(y); % optimal link flow
zopt = value(z); % optimal route flow

% [grad, fval] = threshfunc(xopt, delthresh); % linearize threshold utility function at the current xopt
% valhist(1) = fval;


% for iter = 1:FWmaxiter % Frank-Wolfe iteration
% 
%     yalmip('clear')
% 
%     x = sdpvar(nc, 1, 'full'); % community flow
%     y = sdpvar(nl, 1, 'full'); % link flow
%     z = sdpvar(nr, 1, 'full'); % route flow
%     w = sdpvar(N, 1, 'full'); % dual variable for risk measures
%     nu = sdpvar(1, 1, 'full'); % dual variable for risk measures
%     lam = sdpvar(1, 1, 'full'); % dual variable for risk measures
%     h = sdpvar(N, 1, 'full'); % cap violation
%     u = sdpvar(N, 1, 'full'); % auxiliary variable for expcone
%     tau = sdpvar(1, 1, 'full'); % auxiliary variable for expcone
% 
%     constr = [E*y == zeros(nn, 1), ...
%         x == H*z, ...
%         F*z <= y, ...
%         z >= zeros(nr, 1), ...
%         y >= zeros(nl, 1), ...
%         h >= zeros(N, 1), ...
%         w - nu*ones(N, 1) >= h, ...
%         lam >= 0];
% 
%     for k = 1:N
%         constr = [constr, ...
%             K*y <= (1+h(k))*C(:, k), ...
%             y <= (1+h(k))*D(:, k)];
%     end
% 
%     if strcmp(crmname, 'cvar') % conditional VaR
%         constr = [constr, -nu + (1/(1-del))*lam <= viotol, ...
%             norm(diag(p)*w, 1) <= lam];
%     elseif strcmp(crmname, 'evar') % entropic VaR
%         constr = [constr, -nu + tau -log(1-del)*lam <= viotol, ...
%             p'*u <= lam, ...
%             expcone([w'-tau*ones(1, N); lam*ones(1, N); u'])];
%     elseif strcmp(crmname, 'tv') % total variation
%         constr = [constr, -nu + 2*del*lam+p'*w <= viotol, ...
%             norm(w, Inf) <= lam];
%     else
%         error('Unknown type of risk measure!')
%     end
% 
%     options = sdpsettings('verbose',0,'solver','mosek');
% 
%     grad1 = grad./sum(abs(grad));
% 
%     solution = optimize(constr, -grad1'*x, options);
% 
%     xdir = value(x); % cond grad for x, y, z
%     ydir = value(y);
%     zdir = value(z);
% 
%     xopt = iter/(iter+2)*xopt + 2/(iter+2)*xdir; % Frank-Wolfe update for x, y, z
%     yopt = iter/(iter+2)*yopt + 2/(iter+2)*ydir;
%     zopt = iter/(iter+2)*zopt + 2/(iter+2)*zdir;
% 
%     [grad, fval] = threshfunc(xopt, delthresh);
% 
%     valhist(iter+1) = fval;
%     fprintf('iteration %d, function value %d, rel change %d \n', iter, fval, abs(fval-valhist(iter))/abs(fval))
%     % screen output
% 
% end
% 
% figure('Position',[0,0,800,800]) % plot convergence history of Frank-WOlfe
% plot(1:length(valhist)-1, valhist(2:end))

writematrix(xopt,strcat(crmname,'_community_thresh_', num2str(delthresh), '_delta_', num2str(del), '.csv')) % save variables to .csv files
writematrix(yopt,strcat(crmname, '_link_thresh_', num2str(delthresh), '_delta_', num2str(del), '.csv'))

ginicoeff(xopt)



