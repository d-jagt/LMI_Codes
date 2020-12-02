%% Delay-Independent Robust Stability Condition for Continuous-Time TDS with Slowly-Varying Delay
% Example implementation of LMI from page 68 of:
%   Fridman E. 2014. Introduction to Time-Delay Systems, Analysis and Control. Springer. ISBN: 978-3-319-09392-5.
% Concerns a simple example of a system \dot{x}(t) = (A+H Delta(t) E) x(k) + (A1+H Delta(t) E1) x(t-tau(t))
% tau(t) can vary in time, assuming integer values in [0,h] for some nonnegative integer h
% Delay must be slowly-varying: \dot{tau}(t) <= d < 1
% Uncertainty must be norm-bounded: Delta(t)^T Delta(t) <= I

% Define example system
beta = 0.67;
A = [-2, 0 ; 0, -0.9];
A1 = beta*[-1, 0; -1, -1];
H = 0.1*eye(2,2);
E = ones(2,2);
E1 = beta*[1, 0; 0, -1];
d = 0;

nx = size(A,1);
r1 = size(H,2);
r2 = size(E,1);

% Define the variables
P = sdpvar(nx,nx);
Q = sdpvar(nx,nx);
eps = sdpvar(1,1);

% Construct the LMI
c_mat = [A'*P+P*A+Q , P*A1 , P*H , eps*E';
         A1'*P , -(1-d)*Q , zeros(nx,r1) , eps*E1';
         H'*P , zeros(r1,nx) , -eps*eye(r1) , zeros(r1,r2);
         eps*E , eps*E1 , zeros(r2,r1) , -eps*eye(r2)];
    
eta = 10^-5;
cts = [eps >= eta, P >= eta*eye(nx), Q >= eta*eye(nx)];
cts = [cts, c_mat <= -eta*eye(size(c_mat))];

% Seek a feasible point
obj = [];
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(cts,obj,ops);

% Display results
if sol.problem == 0
    disp(['The TDS is robustly uniformly asymptotically stable for \dot{tau} <= ',num2str(d)])
else
    disp('Stability cannot be shown');
    return
end
