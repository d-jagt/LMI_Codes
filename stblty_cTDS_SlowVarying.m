%% Delay-Independent Robust Stability Condition for Continuous-Time TDS with Slowly-Varying Delay
% Example implementation of LMI from page 64 of:
%   Fridman E. 2014. Introduction to Time-Delay Systems, Analysis and Control. Springer. ISBN: 978-3-319-09392-5.
% Concerns a simple example of a system \dot{x}(t) = A x(k) + A1 x(t-tau(t))
% tau(t) can vary in time, assuming integer values in [0,h] for some nonnegative integer h
% Delay must be slowly-varying: \dot{tau}(t) <= d < 1

% Define example system
beta = 0.35;
A = [-3, -2.5 ; 1, 0.5];
A1 = beta*[1.5, 2.5; -0.5, -1.5];
d = 0.5;

nx = size(A,1);

% Define the variables
P = sdpvar(nx,nx);
Q = sdpvar(nx,nx);

% Construct the LMI
c_mat = [A'*P+P*A+Q , P*A1;
         A1'*P , -(1-d)*Q];
    
eta = 10^-5;
cts = [P >= eta*eye(nx), Q >= eta*eye(nx)];
cts = [cts, c_mat <= -eta*eye(size(c_mat))];

% Seek a feasible point
obj = [];
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(cts,obj,ops);

% Display results
if sol.problem == 0
    disp(['The TDS is uniformly asymptotically stable for \dot{tau} <= ',num2str(d)])
else
    disp('Stability cannot be shown');
    return
end
