% Anesthesia Code
% By: Alvee Hoque and Kah Young

% Constants
STARTTIME = 0;
STOPTIME = 200; %{minutes}
t=STARTTIME:1:STOPTIME;

% Initial conditions
IC=[0,0,0];

[t_out, x] = ode45(@threecomp,t,IC,[]);% Call to ODE45 using fcn hdl below
[~, x_orig] = ode45(@threecomporig,t,IC,[]);% Call to ODE45 using fcn hdl below

figure('DefaultAxesFontSize',18)
subplot(1,2,1)
% multiply by 100 to convert ug to ng
plot(t_out,x_orig(:,1:3)*100,'LineWidth',2)
ylim([0,8])
% Plotting
grid on
%labels
ylabel('Plasma concentration (ng/mL)')
title('Original PBPK Model')
xlabel('Time (min)')
legend('Central','Highly Perfused','Scacrcely Perfused')


% After modifications
subplot(1,2,2)
plot(t_out,x(:,1:3)*100,'LineWidth',2)
ylim([0,8])
grid on
%labels
ylabel('Plasma concentration (ng/mL)')
title('Modified PBPK Model')
xlabel('Time (min)')
legend('Central','Highly Perfused','Scacrcely Perfused')

% x(201,1) %0.035124360427344
% max(x(:,1)) %0.052405605289555

out_change=(x(201,1)-0.035124360427344)/0.035124360427344;
SOF=abs(out_change)/0.05

out_change=(max(x(:,1))-0.052405605289555)/0.052405605289555;
SOF=abs(out_change)/0.05;


 
%%
function dx=threecomp(t,x)

dx = zeros(3,1);

V1=7.88; %mL
V1=7.88*.75; % elderly V1
V2=23.9;
V3=13.8; 
Cl1=2.08; % mL/min
Cl1=2.08*.67; % elderly Cl1
Cl2=0.828;
Cl3=0.0784;
k10=0.172; %min-1
k12=0.373;
k21=0.103;
k13=0.0367;
k31=0.0124;

% Define infusion rate
if t<=20
    I=.25; %ug/kg*min
else
    I=.15;
end

% Central compartment
dx(1)=(-Cl1*x(1)+k21*V2*x(2)+k31*x(3)*V3-((k12+k13+k10)*x(1))*V1+I)/V1;
% Highly perfused compartment
dx(2)=(k12*x(1)*V1-k21*x(2)*V2-Cl2*x(2))/V2;
% Scarcely perfused compartment
dx(3)=(k13*x(1)*V1-k31*x(3)*V3-Cl3*x(3))/V3;

end

%%
function dx=threecomporig(t,x)

dx = zeros(3,1);

V1=7.88; %mL
V2=23.9;
V3=13.8; 
Cl1=2.08; % mL/min
Cl2=0.828;
Cl3=0.0784;
k10=0.172; %min-1
k12=0.373;
k21=0.103;
k13=0.0367;
k31=0.0124;

% Define infusion rate
if t<=20
    I=.25; %ug/kg*min
else
    I=.15;
end

% Central compartment
dx(1)=(-Cl1*x(1)+k21*V2*x(2)+k31*x(3)*V3-((k12+k13+k10)*x(1))*V1+I)/V1;
% Highly perfused compartment
dx(2)=(k12*x(1)*V1-k21*x(2)*V2-Cl2*x(2))/V2;
% Scarcely perfused compartment
dx(3)=(k13*x(1)*V1-k31*x(3)*V3-Cl3*x(3))/V3;

end