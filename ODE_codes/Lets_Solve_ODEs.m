function Lets_Solve_ODEs()

% Create vector of Step Sizes
hVec = [1e-7:1e-7:9e-7 1e-6:1e-6:9e-6 1e-5:1e-5:9e-5 1e-4:1e-4:9e-4 1e-3:1e-3:9e-3 1e-2:1e-2:9e-2 0.1 0.2 0.25 0.5];

% Allocate space for Error Storage
EulerErrMaxVec = zeros(1,length(hVec)); ModEulerErrMaxVec=EulerErrMaxVec; EulerErrL2Vec=EulerErrMaxVec; ModEulerErrL2Vec=EulerErrMaxVec;


for i=1:length(hVec)

    % Simulation Parameters
    tFinal = 1.0;   % final time
    h = hVec(i);    % Step Size
    t = 0:h:tFinal; % time vector 
    y0 = 1.0;       % initial value

    % Perform ODE Solves
    yEulers = give_Me_Euler_Solution(y0,h,t);
    yModEulers = give_Me_Modified_Euler_Solution(y0,h,t);
    yExact = give_Me_Exact_Solution(t);
    [EulerErrMax,ModEulerErrMax,EulerErrL2,ModEulerErrL2] = give_Me_Errors_Please(yEulers,yModEulers,yExact);

    % Save inf-norm and L2-norm errors
    EulerErrMaxVec(i) =    EulerErrMax;
    ModEulerErrMaxVec(i) = ModEulerErrMax;
    EulerErrL2Vec(i)  =    EulerErrL2;
    ModEulerErrL2Vec(i) =  ModEulerErrL2;
    
    % Make Plots of What you Want
    % please_plot_it_all(t,yExact,yEulers,yModEulers);

end

% Plots the error! ("Convergence Plot")
please_plot_the_error(hVec,EulerErrMaxVec,EulerErrL2Vec,ModEulerErrMaxVec,ModEulerErrL2Vec);



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

yEulers';

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
   K1 = 0.5*h*RHS(t(i-1),yModEulers(i-1));                         % Usual Euler Step
   K2 = 0.5*h*( yModEulers(i-1) + h*RHS(t(i-1),yModEulers(i-1)) ); % Approximate Next Step
   yModEulers(i) = yModEulers(i-1) + K1 + K2; 
end

yModEulers';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns Exact Solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yExact = give_Me_Exact_Solution(t)

% TEST
%yExact = exp(t);

% #2
%yExact = cos(pi*t/2);

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
%val = -0.5*pi*sqrt(1-y^2);
%val = y^2 - (pi/2)*sin(pi/2*t) - cos( (pi/2)*t )^2;

% #3
%val = (y-t)^2;

% #5
val = y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes the Error between all cases, returns infinite-norm
% error and L2-norm error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EulerErrorMax,ModEulerErrorMax,EulerErrorL2,ModEulerErrorL2] = give_Me_Errors_Please(yEulers,yModEulers,yExact)

EulerError = abs( yEulers - yExact );
ModEulerError = abs( yModEulers - yExact );


EulerErrorMax = max(EulerError);
EulerErrorL2 = sqrt( EulerError*EulerError' );

ModEulerErrorMax = max(ModEulerError);
ModEulerErrorL2 = sqrt( ModEulerError*ModEulerError' );


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the Convergence Rate of Euler's and Modified Euler's for
% the inf-norm and L2-norm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_plot_the_error(hVec,EulerErrMaxVec,EulerErrL2Vec,ModEulerErrMaxVec,ModEulerErrL2Vec)

lw = 3;  % LineWidth
ms = 10; % MarkerSize
fs = 22;  % FontSize

%
% Plots Inf-Norm Error against each other
%
subplot(1,2,1);
loglog(hVec,EulerErrMaxVec,'r-','LineWidth',lw); hold on;
loglog(hVec,ModEulerErrMaxVec,'b-','LineWidth',lw); hold on;
loglog(hVec,EulerErrMaxVec,'ro','LineWidth',lw,'MarkerSize',ms); hold on;
loglog(hVec,ModEulerErrMaxVec,'b.','MarkerSize',ms+24); hold on;
title('Convergence Plot: inf-norm','FontSize',fs);hold on;
xlabel('h','FontSize',fs);
ylabel('Inf-Norm Error','FontSize',fs);
leg=legend('Eulers','Mod. Eulers');
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);

%
% Plots L2-Norm Error against each other
%
subplot(1,2,2);
loglog(hVec,EulerErrL2Vec,'r-','LineWidth',lw); hold on;
loglog(hVec,ModEulerErrL2Vec,'b-','LineWidth',lw); hold on;
loglog(hVec,EulerErrL2Vec,'ro','LineWidth',lw,'MarkerSize',ms); hold on;
loglog(hVec,ModEulerErrL2Vec,'b.','MarkerSize',ms+24); hold on;
title('Convergence Plot: L2-norm','FontSize',fs);hold on;
xlabel('h','FontSize',fs);
ylabel('L2-Norm Error','FontSize',fs);
leg=legend('Eulers','Mod. Eulers');
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);

