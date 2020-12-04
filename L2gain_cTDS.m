%% Bound on L2 gain for Continuous-Time TDS with Slowly-Varying Delay
% Example implementation of LMI from page 151 of:
%   Fridman E. 2014. Introduction to Time-Delay Systems, Analysis and Control. Springer. ISBN: 978-3-319-09392-5.
% Concerns a simple example of a system:
%   \dot{x}(t) = A x(t) + A1 x(t-tau(t)) +B0 w(t)
%   z(t) = C0 x(t) + C1 x(t-tau(t))
% tau_k can vary in time, assuming integer values in [0,h] for some nonnegative integer h
% A bound d on \dot{\tau}<=1 can be implemented, or set greater than or equal to 1 if unknown

% Define example system
A = [-2, 0 ; 0, -0.9];
A1 = [-1, 0; -1, -1];
B0 = [-0.5; 1];
C0 = [1, 0];
C1 = [0, 0];
h = 0.846;
d = 0.5;

nx = size(A,1);
nw = size(B0,2);
nz = size(C0,1);

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
P3 = sdpvar(nx,nx,'full');
S12 = sdpvar(nx,nx,'full');

gmm2 = sdpvar(1,1);

% Construct the LMI
Phi11 = A'*P2 + P2'*A + S  + Q - R;
Phi12 = P - P2' + A'*P3;
Phi22 = -P3 - P3' + (h^2)*R;

Xi = [Phi11 , Phi12 , S12 , R-S12+P2'*A1;
      Phi12' , Phi22 , zeros(nx) , P3'*A1;
      S12' , zeros(nx) , -(S+R) , R-S12';
      R'-S12'+A1'*P2 , A1'*P3 , R'-S12 , -2*R+S12+S12'-(1-d)*Q];
  
Xi_12 = [P2'*B0 , C0';
         P3'*B0 , zeros(nx,nz);
         zeros(nx,nw) , zeros(nx,nz);
         zeros(nx,nw) , C1'];
     
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
    disp(['The L_2-gain of the system is less than ',num2str(gamma)])
else
    disp('No L_2-bound obtained');
    return
end
