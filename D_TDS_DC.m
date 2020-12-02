%% Delay-Dependent Stability Condition for Discrete-Time TDS
% Example implementation of LMI from page 250 of:
%   Fridman E. 2014. Introduction to Time-Delay Systems, Analysis and Control. Springer. ISBN: 978-3-319-09392-5.
% Concerns a simple example of a system x(k+1) = A x(k) + A1 x(k-tau_k)
% tau_k can vary in time, assuming integer values in [0,h] for some nonnegative integer h

% Define example system
A = [0.8, 0 ; 0, 0.97];
A1 = [-0.1, 0; -0.1, -0.1];
h = 15;

nx = size(A,1);

% Define the variables
P = sdpvar(nx,nx);
R = sdpvar(nx,nx);
S = sdpvar(nx,nx);

P2 = sdpvar(nx,nx,'full');
P3 = sdpvar(nx,nx,'full');
S12 = sdpvar(nx,nx,'full');

% Construct the LMI
Phi11 = (A'-eye(nx))*P2 + P2'*(A-eye(nx)) + S - R;
Phi12 = P - P2' + (A'-eye(nx))*P3;
Phi22 = -P3 - P3' + P + (h^2)*R;

c_mat = [Phi11 , Phi12 , S12 , R-S12+P2'*A1;
        Phi12' , Phi22 , zeros(nx) , P3'*A1;
        S12' , zeros(nx) , -(S+R) , R-S12';
        R'-S12'+A1'*P2 , A1'*P3 , R'-S12 , -2*R+S12+S12'];
    
eta = 10^-5;
cts = [P >= eta*eye(nx), R >= eta*eye(nx), S >= eta*eye(nx)];
cts = [cts, c_mat <= -eta*eye(size(c_mat))];

% Seek a feasible point
obj = [];
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(cts,obj,ops);

% Display results
if sol.problem == 0
    disp(['The TDS is asymptotically stable for h <= ',num2str(h)])
else
    disp('Stability cannot be shown');
    return
end
