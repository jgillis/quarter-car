N = 100;

% parameters
ms = 400;   	% sprung mass -> car chassis
mus = 40;       % unsprung mass -> wheel assembly
ks = 15791;     % spring coefficient suspension system
kus = 157910;   % compressibility pneumatic tire
cs = 1508;      % damping coefficient suspension system
cus = 400;      % damping pneumatic tire (41% of cs)

% u(1) = active control force suspension
% u(2) = derivative of road displacement input -> road velocity disturbance
% x(1) = tire deflection
% x(2) = unsprung-mass velocity
% x(3) = suspension stroke
% x(4) = sprung-mass velocity

% ode
ode = @(x,u) [x(2)-u(2);
              -kus/mus*x(1)-(cs+cus)/mus*x(2)+ks/mus*x(3)+cs/mus*x(4)+1/mus*u(1)+cus/mus*u(2);
              -x(2)+x(4);
              cs/ms*x(2)-ks/ms*x(3)-cs/ms*x(4)-1/ms*u(1)];

% states
x = optivar(N+1,4,'x');

% control: 
U = optivar(N,1,'U');

% X0 = optipar('X0',4,1);
% wegdek = optipar('w',N,1);

% a variable to denote the final time
tf = 1;

% step
% z0      = @(t) 0.*t.*(t<0.1)-10.*(t-0.1).*(t>=0.1 & t<0.11)-0.1.*(t>=0.11);
% z0_punt = @(t) 0.*t.*(t<0.1)-10.*(t>=0.1 & t<=0.11)-0.*(t>0.11);

% 1 bump
z0      = @(t) 0.*t.*(t<0.25)+0.05*(1-cos(8*pi.*t))/2.*(t>=0.25 & t<0.5)+0.*t.*(t>=0.5);
z0_punt = @(t) 0.*t.*(t<0.25)+0.05*(8*pi*sin(8*pi.*t))/2.*(t>=0.25 & t<0.5)+0.*t.*(t>=0.5);

% 2 bumps
% z0     = @(t) 0.*t.*(t<0.25)+0.05*(1-cos(8*pi.*t))/2.*(t>=0.25 & t<0.5)+0.*t.*(t>=0.5 & t<0.75)-0.05*(1-cos(8*pi.*t))/2.*(t>=0.75 & t<1);
% z0_punt = @(t) 0.*t.*(t<0.25)+0.05*(8*pi*sin(8*pi.*t))/2.*(t>=0.25 & t<0.5)+0.*t.*(t>=0.5 & t<0.75)-0.05*(8*pi*sin(8*pi.*t))/2.*(t>=0.75 & t<1);

% Construct list of all constraints
g = {};

f = 0;
for k=1:N
   xk      = x(k,:)';
   xk_plus = x(k+1,:)';
   
   % shooting constraint
   uk = [U(k);z0_punt(k*tf/N)];
   
   xf = rk4(ode,tf/N,xk,uk);
   g = {g{:}, xk_plus==xf};

   y=ode(xk,uk);
   f = f + (y(4,1)^2+1100*xk(1,1)^2+100*xk(3,1)^2)*tf/N;
end

% path constraint
% constr = @(x2) 1-sin(2*pi*x2)/2;

g = {g{:}, x(:,3) <= 1, x(:,3) >= -1 };  %constr(x2)

U.setLb(-3000);
U.setUb(3000);

g = {g{:}, x(1,:)==0, x(end,:)==0};

disp('solving problem')
optisolve(f,g,struct('expand',true));

% X0.setValue([0;0;0;0])
% sol = optisolve(f,g,struct('expand',true));

% for i=1:1000
%  X0.setValue([0;0;0;0])
%  wegdek.setValue();
%  sol.resolve();
% end

figure
subplot(3,1,1)
hold on
plot(optival(x(:,2)),'m');
plot(optival(x(:,4)),'r');
plot(d,z0_punt(d*tf/N),'k');
legend('unsprung-mass velocity [m/s]','sprung-mass velocity [m/s]','road velocity disturbance [m/s]');
axis([0 N -inf inf]);

subplot(3,1,2)
hold on
plot(optival(x(:,1)),'y');
plot(optival(x(:,3)),'c');
plot(d,z0(d*tf/N),'k');
%plot(0.05*ones(size(optival(x(:,1))))+optival(x(:,1))+optival(x(:,3))+z0(d*tf/N)')
legend('tire deflection [m]','suspension stroke [m]','road disturbance [m]');
axis([0 N -inf inf]);

subplot(3,1,3)
hold on
stairs(optival(U))
legend('active control force suspension [N]')
axis([0 N -inf inf]);
