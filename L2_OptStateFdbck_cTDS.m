%% L2-Optimal Full-State-Feedback for Continuous-Time TDS (with Slowly-Varying Delay)
% Example implementation of LMI from page 211 of:
%   Fridman E. 2014. Introduction to Time-Delay Systems, Analysis and Control. Springer. ISBN: 978-3-319-09392-5.
% Concerns a simple example of a system:
%   \dot{x}(t) = A x(t) + B2 u(t-tau(t)) +B1 w(t)
%   z(t) = C1 x(t) + D12 u(t-tau(t))
% tau_k can vary in time, assuming integer values in [0,h] for some nonnegative integer h
% A bound d on \dot{\tau}<=1 can be implemented, or set greater than or equal to 1 if unknown
% Feedback control u(t)=Kx(t) is sought to minimize L2-gain of closed-loop system

% Define example system
A = [1, 0.5 ; 0.1, -1];
B2 = [1+0.3; -1];
B1 = [0; 1];
C1 = [0, 1];
D12 = 0.1;
h = 0.3;
d = 1;
eps = 0.7;

nx = size(A,1);
nw = size(B1,2);
nu = size(B2,2);
nz = size(C1,1);

% Define the variables
P = sdpvar(nx,nx);
R = sdpvar(nx,nx);
S = sdpvar(nx,nx);
if d>=1
    Q = zeros(nx,nx);
    cts = [];
else
    eta = 10^-5;
    Q = sdpvar(nx,nx);
    cts = [Q >= eta*eye(nx)];
end
P2 = sdpvar(nx,nx,'full');

Y = sdpvar(nu,nx,'full');
S12 = sdpvar(nx,nx,'full');

gmm2 = sdpvar(1,1);

% Construct the LMI
Phi11 = A*P2 + P2'*A' + S  + Q - R;
Phi12 = P - P2 + eps*P2'*A';
Phi22 = -eps*P2 - eps*P2' + (h^2)*R;

Xi = [Phi11 , Phi12 , S12 , B2*Y + R - S12;
      Phi12' , Phi22 , zeros(nx) , eps*B2*Y;
      S12' , zeros(nx) , -(S+R) , R-S12';
      (B2*Y + R - S12)' , eps*Y'*B2' , R'-S12 , -(1-d)*Q - 2*R + S12 + S12'];
  
Xi_12 = [B1 , P2'*C1';
         eps*B1 , zeros(nx,nz);
         zeros(nx,nw) , zeros(nx,nz);
         zeros(nx,nw) , Y'*D12'];
     
Xi_22 = [-gmm2*eye(nw,nw) , zeros(nw,nz);
         zeros(nz,nw) , -eye(nz)];

c_mat = [Xi , Xi_12;
         Xi_12' , Xi_22];
    
eta = 10^-5;
cts = [cts, P >= eta*eye(nx), R >= eta*eye(nx), S >= eta*eye(nx)];
cts = [cts, c_mat <= -eta*eye(size(c_mat))];

% Seek a feasible point while minimizing gamma^2
obj = gmm2;
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(cts,obj,ops);

% Display results
if sol.problem == 0
    gamma = sqrt(value(gmm2));
    K = value(Y)/value(P2)
    disp(['The L_2-gain of the closed-loop system is less than ',num2str(gamma)])
else
    disp('No optimal controller attained');
    return
end
