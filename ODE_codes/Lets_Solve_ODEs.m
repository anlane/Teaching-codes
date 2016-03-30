function Lets_Solve_ODEs()


% Simulation Parameters
tFinal = 1.0;   % final time
h = 0.1;        % Step Size
t = 0:h:tFinal; % time vector 
y0 = 1.0;       % initial value


% Perform ODE Solves
yEulers = give_Me_Euler_Solution(y0,h,t);
yModEulers = give_Me_Modified_Euler_Solution(y0,h,t);
yExact = give_Me_Exact_Solution(t);


% Make Plots of What you Want
please_plot_it_all(t,yExact,yEulers,yModEulers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns Euler Solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yEulers = give_Me_Euler_Solution(y0,h,t)

% Initialize Solution Storage 
yEulers = zeros(1,length(t));
yEulers(1) = y0;

for i=2:length(t)
   yEulers(i) = yEulers(i-1) + h * RHS(t(i-1),yEulers(i-1)); 
end

yEulers'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns Modified Euler Solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yModEulers = give_Me_Modified_Euler_Solution(y0,h,t)

% Initialize Solution Storage 
yModEulers = zeros(1,length(t));
yModEulers(1) = y0;

for i=2:length(t)
   K1 = 0.5*h*RHS(t(i-1),yModEulers(i-1));                             %Usual Euler Step
   K2 = 0.5*h*( yModEulers(i-1) + 0.5*h*RHS(t(i-1),yModEulers(i-1)) ); %Approximate Next Step
   yModEulers(i) = yModEulers(i-1) + K1 + K2; 
end

yModEulers'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns Modified Euler Solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yExact = give_Me_Exact_Solution(t)

% TEST
%yExact = exp(t);

% #2
%yExact = sin(pi*t/2);

% #3
%yExact = ( exp(2*t).*(t-1)+t+1 ) ./ ( exp(2*t) + 1 );

% %5
yExact = exp(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns RHS of ODE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = RHS(t,y)

% TEST
%val = y;

% #2
%val = 0.5*pi*sqrt(1-y^2);

% #3
%val = (y-t)^2;

% #5
val = y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the Exact Solution, Euler's Method Solution, and
% Modified Euler's Solution. Also produces a plot of the Error.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_plot_it_all(t,Exact,yEulers,yModEulers)

lw = 3;  % LineWidth
ms = 10; % MarkerSize
fs = 22;  % FontSize

% Find more highly resolved exact solution
tNew = t(1):( t(2)-t(1) )/100:t(end);
fineExact = give_Me_Exact_Solution(tNew);

%
% Plotting Solutions Against Each Other
%
subplot(1,2,1);
plot(tNew,fineExact,'k-','LineWidth',lw); hold on;
plot(t,yEulers,'b-','LineWidth',lw); hold on;
plot(t,yModEulers,'r-','LineWidth',lw);hold on;
plot(t,Exact,'ko','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,yEulers,'b.','MarkerSize',ms+24); hold on;
plot(t,yModEulers,'r.','MarkerSize',ms+24);hold on;
title('Numerical and Exact Solutions','FontSize',fs);hold on;
xlabel('x','FontSize',fs);
ylabel('Exact Soln. and Numerical Solns.','FontSize',fs);
leg=legend('Exact','Eulers','Mod. Eulers');
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);
%
% Plotting Errors Against Each Other
%
subplot(1,2,2);
plot(t,abs(Exact-yEulers),'b.-','LineWidth',lw,'MarkerSize',ms+25); hold on;
plot(t,abs(Exact-yModEulers),'r.-','LineWidth',lw,'MarkerSize',ms+25); hold on;
title('Error Comparison','FontSize',fs);
xlabel('x','FontSize',fs);
ylabel('| Exact - Numerical |','FontSize',fs);
leg=legend('Eulers Method Error','Mod. Eulers Method Error');
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);
