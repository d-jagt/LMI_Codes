%% This script provides a demo of the implementation of D-stability condtions for full-state-feedback
% The D-stability region is based on desired transient response metric values
% The transient response conditions are based on those presented in:
%   G.F. Franklin, J.D. Powell, and A. Emami-Naeini, Feedback Control of Dynamic Systems. Upper Saddle River, NJ: Prentice Hall, 2002. 
% The D-region is approximated as the intersection of several LMI regions for efficient implementation

clear

%% Construct the system:

A = [1,0,0,-1,0;1,-1,-1,0,1;0,0,-1,0,1;1,1,1,2,1;1,2,0,2,1];
B1 = [0;0;0;0;1];
B2 = [0;0;0;0;1];
C1 = [1,0,0,0,0];

n = size(A,1);
nw = size(B1,2);
nu = size(B2,2);
nz = size(C1,1);

sys = ss(A,B1,C1,0,1);
%step(sys)


%% Define the necessary parameters
kr = 1;     % Minimal rise time
ks = 25;    % Maximal settling time
Mp = 0.4;   % Maximal overshoot

p_s = exp(-4.6/ks);

p_r = sinh(1.8/kr);
q_r = cosh(1.8/kr);

wmax = atan(-pi/log(Mp));
xmax = (Mp^(wmax/pi))*cos(wmax);
ymax = (Mp^(wmax/pi))*sin(wmax);
a_M = 0.5*(1+Mp);
b_M = ymax;
q_M = 0.5*(1-Mp);
c_M = ymax/sqrt(1-xmax);


%% Build the LMI

% Construct the variables
P = sdpvar(n,n);
Z = sdpvar(nu,n);

% Define the constraints
eta = 10^-5;

s_mat = [-p_s*P , A*P+B2*Z;
         (A*P+B2*Z)' , -p_s*P];
     
r_mat = [-p_r*P , -q_r*P + A*P+B2*Z;
         -q_r*P + (A*P+B2*Z)' , -p_r*P];
     
M_mat1 = [-2*a_M*b_M*P , -2*b_M*q_M*P + (b_M+a_M)*(A*P+B2*Z) + (b_M-a_M)*(A*P+B2*Z)';
          -2*b_M*q_M*P + (b_M+a_M)*(A*P+B2*Z)' + (b_M-a_M)*(A*P+B2*Z) , -2*a_M*b_M*P];
      
M_mat2 = [-2*P + (A*P+B2*Z) + (A*P+B2*Z)' , (A*P+B2*Z) - (A*P+B2*Z)';
          (A*P+B2*Z)' - (A*P+B2*Z) , -2*c_M*c_M*P];

cts = [P >= eta*eye(n)];
cts = [cts, s_mat <= -eta*eye(size(s_mat))];
cts = [cts, r_mat <= -eta*eye(size(r_mat))];
cts = [cts, M_mat1 <= -eta*eye(size(M_mat1))];
cts = [cts, M_mat2 <= -eta*eye(size(M_mat2))];

% Solve the LMI
obj = [];
ops = sdpsettings('solver','mosek','verbose',0);
sol = optimize(cts,obj,ops);

% Display results
if sol.problem == 0
    K = value(Z)/value(P)
else
    disp('No feasible controller attained');
    return
end


%% Plot the step response

Acl = A+B2*K;
sysK = ss(Acl,B1,C1,0,1);
step(sysK);


%% Plot the locations of the poles if xtra==1

xtra = 1;

if xtra==1
    
z = eigs(Acl);

% Plot the poles along with the actual D-region

wdneg = linspace(-pi,0,100);
wdpos = linspace(0,pi,100);
wd = [wdneg,wdpos];
wd_max = min(pi,(1.8/kr));
wd_alt = linspace(-wd_max,wd_max,200);

bnd_r = exp(sqrt((1.8/kr)^2 - wd_alt.^2));
bnd_rm = exp(-sqrt((1.8/kr)^2 - wd_alt.^2));
bnd_s = exp(-4.6/ks);
bnd_M = Mp.^(abs(wd)/pi);

Re_bnd_r = bnd_r.*cos(wd_alt);
Im_bnd_r = bnd_r.*sin(wd_alt);
Re_bnd_rm = bnd_rm.*cos(wd_alt);
Im_bnd_rm = bnd_rm.*sin(wd_alt);
Re_bnd_s = bnd_s.*cos(wd);
Im_bnd_s = bnd_s.*sin(wd);
Re_bnd_M = bnd_M.*cos(wd);
Im_bnd_M = bnd_M.*sin(wd);

fig2 = figure(2);
plot(cos(wd),sin(wd),'k--');
hold on
xline(0,'k');
yline(0,'k');
plot(Re_bnd_r(1:end),Im_bnd_r(1:end),'-','color',[0 0.2 0.8],'Displayname','k_r');
plot(Re_bnd_rm(1:end),Im_bnd_rm(1:end),'-','color',[0 0.2 0.8]);
plot(Re_bnd_s(1:2:end),Im_bnd_s(1:2:end),'-','color',[0 0.8 0.2],'Displayname','k_s');
plot(Re_bnd_M(1:2:end),Im_bnd_M(1:2:end),'-','color',[0.8 0.1 0.1],'Displayname','M_p');
scatter(real(z),imag(z),'ko','Filled');
hold off
xlabel('Re(z)');
ylabel('Im(z)');
g = get(fig2,'Children');
f = g.Children;
legend([f(5),f(3),f(2)]);
axis([-1 2.5 -1.75 1.75]);
fig2.Position = [450 500 550 500];


% Plot the poles along with the imposed LMI region

Re_bnd_r_LMI = p_r.*cos(wd)+q_r;
Im_bnd_r_LMI = p_r.*sin(wd);

x_M_LMI1 = linspace(q_M-a_M,q_M+a_M,100);
y_M_LMI1 = ymax*sqrt(1-((x_M_LMI1-q_M)/a_M).^2);
y_M_LMI2 = c_M*sqrt(1-x_M_LMI1);

fig3 = figure(3);
plot(cos(wd),sin(wd),'k--');
hold on
xline(0,'k');
yline(0,'k');
plot(Re_bnd_r_LMI,Im_bnd_r_LMI,'-','color',[0 0.2 0.8],'Displayname','k_r LMI');
plot(Re_bnd_s,Im_bnd_s,'-','color',[0 0.8 0.2],'Displayname','k_s LMI');
plot(x_M_LMI1,y_M_LMI1,'-','color',[0.8 0.1 0.1],'Displayname','M_p LMI');
plot(x_M_LMI1,-y_M_LMI1,'-','color',[0.8 0.1 0.1],'Displayname','M_p LMI');
plot(x_M_LMI1,y_M_LMI2,'-','color',[0.8 0.1 0.1],'Displayname','M_p LMI');
plot(x_M_LMI1,-y_M_LMI2,'-','color',[0.8 0.1 0.1],'Displayname','M_p LMI');
scatter(real(z),imag(z),'ko','Filled');
hold off

xlabel('Re(z)');
ylabel('Im(z)');
g = get(fig3,'Children');
f = g.Children;
legend([f(7),f(6),f(2)]);
axis([-1 2.5 -1.75 1.75]);
fig3.Position = [450 500 550 500];

end