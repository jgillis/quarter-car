import casadi.*

tf = 1;                 % a variable to denote the final time in sec
N = 100*tf;             % 100Hz
speedkmh = 60;          % speed in km/h
speed = speedkmh/3.6;   % speed in m/s
horizon = tf*speed;     % preview length

% parameters
ms = 400;   	% sprung mass -> car chassis
mus = 40;       % unsprung mass -> wheel assembly
ks = 15791;     % spring coefficient suspension system
kus = 157910;   % compressibility pneumatic tire
cs = 1508;      % damping coefficient suspension system
cus = 400;      % damping pneumatic tire (41% of cs)
Ns = 1600;      % nonlinear stiffness property of primary suspension

% u(1) = active control force suspension
% u(2) = derivative of road displacement input -> road velocity disturbance
% x(1) = tire deflection
% x(2) = unsprung-mass velocity
% x(3) = suspension stroke
% x(4) = sprung-mass velocity

% ode
ode = @(x,u) [x(2)-u(2);
              -kus/mus*x(1)-(cs+cus)/mus*x(2)+ks/mus*x(3)+cs/mus*x(4)+1/mus*u(1)+cus/mus*u(2)-Ns/mus*(x(3))^3;
              -x(2)+x(4);
              cs/ms*x(2)-ks/ms*x(3)-cs/ms*x(4)-1/ms*u(1)-Ns/mus*(x(3))^3];

% states
X = optivar(N+1,4,'x');

% control: 
U = optivar(N,1,'U');

% X0 = optipar('X0',4,1);
% wegdek = optipar('w',N,1);

% step
% z0      = @(t) 0.*t.*(t<0.1)-10.*(t-0.1).*(t>=0.1 & t<0.11)-0.1.*(t>=0.11);
% z0_punt = @(t) 0.*t.*(t<0.1)-10.*(t>=0.1 & t<=0.11)-0.*(t>0.11);

% 1 bump
% z0      = @(t) 0.*t.*(t<0.25)+0.05*(1-cos(8*pi.*t))/2.*(t>=0.25 & t<0.5)+0.*t.*(t>=0.5);
% z0_punt = @(t) 0.*t.*(t<0.25)+0.05*(8*pi*sin(8*pi.*t))/2.*(t>=0.25 & t<0.5)+0.*t.*(t>=0.5);

% sinus van 1.5Hz
f = 1.5;       %Hz
z0      = @(t) 0.015*(1-cos(2*pi.*t*f))/2;
z0_punt = @(t) 0.015*2*pi*f*sin(2*pi.*t*f)/2;

% 1 bump (verkeersdrempeltje)
% z0      = @(t) 0.*t.*(t<0.25)+0.05*(1-cos(2*pi.*40.*t))/2.*(t>=0.25 & t<0.275)+0.*t.*(t>=0.275);
% z0_punt = @(t) 0.*t.*(t<0.25)+0.05*(2*40*pi*sin(2*pi.*40.*t))/2.*(t>=0.25 & t<0.275)+0.*t.*(t>=0.275);

% 2 bumps
% z0     = @(t) 0.*t.*(t<0.25)+0.05*(1-cos(8*pi.*t))/2.*(t>=0.25 & t<0.5)+0.*t.*(t>=0.5 & t<0.75)-0.05*(1-cos(8*pi.*t))/2.*(t>=0.75 & t<1);
% z0_punt = @(t) 0.*t.*(t<0.25)+0.05*(8*pi*sin(8*pi.*t))/2.*(t>=0.25 & t<0.5)+0.*t.*(t>=0.5 & t<0.75)-0.05*(8*pi*sin(8*pi.*t))/2.*(t>=0.75 & t<1);


% Create CasADi function, for speed of initialization

x = SX.sym('x',4);
u = SX.sym('u',2);
simulate_one_interval = Function('rk4', {x,u},{rk4(ode,tf/N,x,u)});

y = ode(x,u);
cost = Function('cost',{x,u},{(y(4,1)^2+1100*x(1,1)^2+100*x(3,1)^2)*tf/N});

% Construct list of all constraints
g = {};

f = 0;
for k=1:N
   xk      = X(k,:)';
   xk_plus = X(k+1,:)';
   
   % shooting constraint
   uk = [U(k);z0_punt((k-1)*tf/N)];
   
   xf = simulate_one_interval(xk,uk);
   g = {g{:}, xk_plus==xf};

   f = f + cost(xk,uk);
end

% path constraint
% constr = @(x2) 1-sin(2*pi*x2)/2;

g = {g{:}, X(:,3) <= 0.2, X(:,3) >= -0.2 };  %constr(x2)

U.setLb(-3000);
U.setUb(3000);

g = {g{:}, X(1,:)==0};

disp('solving problem')
optisolve(f,g,struct('expand',true));

% X0.setValue([0;0;0;0])
% sol = optisolve(f,g,struct('expand',true));

% for i=1:1000
%  X0.setValue([0;0;0;0])
%  wegdek.setValue();
%  sol.resolve();
% end

d=[1:N];

figure
plot1 = subplot(3,1,1);
hold on
plot(optival(X(:,2)),'m');
plot(optival(X(:,4)),'r');
plot(d,z0_punt((d-1)*tf/N),'k');
legend('unsprung-mass velocity [m/s]','sprung-mass velocity [m/s]','road velocity disturbance [m/s]');
axis([0 N -inf inf]);

plot2 = subplot(3,1,2);
hold on
plot(optival(X(:,1)),'y');
plot(optival(X(:,3)),'c');
plot(d,z0((d-1)*tf/N),'k');
%plot(0.05*ones(size(optival(x(:,1))))+optival(x(:,1))+optival(x(:,3))+z0(d*tf/N)')
legend('tire deflection [m]','suspension stroke [m]','road disturbance [m]');
axis([0 N -inf inf]);

plot3 = subplot(3,1,3);
hold on
stairs(optival(U))
legend('active control force suspension [N]')
axis([0 N -inf inf]);

first_axis = gca;
sqz = 0.05; %// distance to squeeze the first plot
set(first_axis, 'Position', get(first_axis, 'Position') + [0 sqz 0 -sqz ]);
ax2 = axes('Position', get(first_axis, 'Position') .* [1 1 1 0.001] - [0 sqz 0 0],'Color','none');
scale_factor = horizon/N; %// change this to your satisfaction
xlim(get(first_axis, 'XLim') * scale_factor);
set(ax2, 'XScale', get(first_axis, 'XScale')); %// make logarithmic if first axis is too
xlabel(ax2,'m');

linkaxes([plot1,plot2,plot3],'x');
