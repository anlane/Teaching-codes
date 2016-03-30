function Lets_Solve_ODE_Systems()


% Simulation Parameters
tFinal = 1.0;   % final time
h = 0.2;        % Step Size
t = 0:h:tFinal; % time vector 
y0(1,1) = 3.0;  % initial value for y1
y0(1,2) = 0.0;  % initial value for y2


% Perform ODE Solves
yEulers = give_Me_Euler_Solution(y0,h,t);
yExact = give_Me_Exact_Solution(t);


% Make Plots of What you Want
please_plot_it_all(t,yExact,yEulers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns Euler Solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yEulers = give_Me_Euler_Solution(y0,h,t)

% yEulers: row- corresponds to time-step
%          col- col1: y1, col2: y2

% Initialize Solution Storage 
yEulers = zeros(length(t),2);
yEulers(1,1) = y0(1,1);
yEulers(1,2) = y0(1,2);

for i=2:length(t)
   RHS = give_Me_RHS_of_System(t(i-1),yEulers(i-1,1),yEulers(i-1,2));
   yEulers(i,1) = yEulers(i-1,1) + h * RHS(1); 
   yEulers(i,2) = yEulers(i-1,2) + h * RHS(2); 
end

yEulers



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns Modified Euler Solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yExact = give_Me_Exact_Solution(t)

% yExact: row- corresponds to time-step
%          col- col1: y1, col2: y2

% #1
yExact(:,1) = 4*exp(t) - exp(-2*t);
yExact(:,2) = exp(t)  - exp(-2*t);

% #2
%yExact(:,1) = 4*exp(-t).*sin(t);
%yExact(:,2) = 4*exp(-t).*cos(t);

% #3
%yExact(:,1) = cos(t/2);
%yExact(:,2) = -0.5*sin(t/2);

% #5
%yExact(:,1) = exp(-t) - t;
%yExact(:,2) = -exp(-t) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns RHS of ODE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RHS = give_Me_RHS_of_System(t,y1,y2)

% #1
RHS(1) = 2*y1 - 4*y2;
RHS(2) = y1 - 3*y2;

% #2
%RHS(1) = -y1+y2;
%RHS(2) = -y1-y2;

% #3
%RHS(1) = 0*y1 + 1*y2;
%RHS(2) = -1/4*y1 + 0*y2;

% #5
%RHS(1) = 0*y1 + 1*y2;
%RHS(2) = 1*y1 + 0*y2 + t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the Exact Solution, Euler's Method Solution, and
% Modified Euler's Solution. Also produces a plot of the Error.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_plot_it_all(t,Exact,yEulers)

lw = 3;  % LineWidth
ms = 10; % MarkerSize
fs = 22;  % FontSize

% Find more highly resolved exact solution
tNew = t(1):( t(2)-t(1) )/100:t(end);
fineExact = give_Me_Exact_Solution(tNew);

%
% Plotting Solutions Against Each Other
%
figure(1)
subplot(1,2,1);
plot(tNew,fineExact(:,1),'k-','LineWidth',lw); hold on;
plot(tNew,fineExact(:,2),'g-','LineWidth',lw); hold on;
plot(t,yEulers(:,1),'b-','LineWidth',lw); hold on;
plot(t,yEulers(:,2),'r-','LineWidth',lw); hold on;
plot(t,Exact(:,1),'ko','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,Exact(:,2),'go','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,yEulers(:,1),'b.','MarkerSize',ms+24); hold on;
plot(t,yEulers(:,2),'r.','MarkerSize',ms+24); hold on;
title('Numerical and Exact Solutions','FontSize',fs);hold on;
xlabel('x','FontSize',fs);
ylabel('Exact Soln. and Numerical Solns.','FontSize',fs);
leg=legend('y1 Exact','y2 Exact','y1 Eulers', 'y2 Eulers');
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);
%
% Plotting Errors Against Each Other
%
subplot(1,2,2);
plot(t,abs(Exact(:,1)-yEulers(:,1)),'b.-','LineWidth',lw,'MarkerSize',ms+50); hold on;
plot(t,abs(Exact(:,2)-yEulers(:,2)),'r.-','LineWidth',lw,'MarkerSize',ms+25); hold on;
title('Error Comparison','FontSize',fs);
xlabel('x','FontSize',fs);
ylabel('| Exact - Numerical |','FontSize',fs);
leg=legend('y1 Error','y2 Error');
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);


% PHASE PLANE
figure(2)
plot(fineExact(:,1),fineExact(:,2),'k-','LineWidth',lw); hold on;
plot(yEulers(:,1),yEulers(:,2),'b-','LineWidth',lw); hold on;
plot(Exact(:,1),Exact(:,2),'ko','LineWidth',lw,'MarkerSize',ms); hold on;
plot(yEulers(:,1),yEulers(:,2),'b.','MarkerSize',ms+24); hold on;
title('Phase Plane: Numerical & Exact Solns','FontSize',fs);hold on;
xlabel('y1','FontSize',fs);
ylabel('y2','FontSize',fs);
leg=legend('Exact','Eulers Soln');
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);
