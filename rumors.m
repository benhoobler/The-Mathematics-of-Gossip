function [x,t] = SIR_Adapt(Kappa_input, Lambda_input, Alpha_input, Susceptibles, Spreaders, Stiflers, Plot_Option)
close all;
% Ben Hoobler
% The Mathematics of Gossip


%% The Mathematics of Gossip

% In this project, the goal is to create a numerical model to investigate
% the spread of rumors, lies, and misinformation that take place throughout
% the world today. The model will be an adapted SIR model used similar to
% the model of the transmission and vaccination of the Ebola virus.
% Different parameters will be tested, such as susceptibility of
% "infection" depending on the source of the information and the
% probability that an infected host can overcome the rumor and believe it
% to be false.


%% Recommended Initials for t = 0

% Susceptibles = .99
% Spreaders = .01
% Hibernators = 0
% Stiflers = 0


%%%%%%%%%%%%
%Quick Copy
%%%%%%%%%%%%
%SIR_Adapt(.3, .8, .3, .99, .01, 0,1);

%% Variables
global Alpha;
Alpha = Alpha_input;
global Kappa;
Kappa = Kappa_input;
global Lambda;
Lambda = Lambda_input;

%% Solver Input
Initials = [Susceptibles, Spreaders, Stiflers]; %for t = 0
TimeScale = [0:.1:100];

%% Solver
[t,x] = ode45( @system, TimeScale, Initials);

%% Plot
if(Plot_Option)
figure
plot(t,x(:,1),'g',t,x(:,2),'b', t,x(:,3),'k');
title('SIR Model of Rumor Spreading')
legend('Ignorant', 'Spreader', 'Stifler');
xlabel('time [hours]')
ylabel('Density of Classes in SIR Model')
end

end

function output = system(t,x)
% This function is to store all the equations
    
 %% Factors of Spreading
% Kappa is the degree of the network. If the origin came from a Twitter
% account with a very large following, Kappa would be much larger than an
% interation between friends on Facebook. 

% Lamda is the probability that the person disseminates the rumor and
% changes into a Spreader. Lamda is the spreading rate of a rumor

% Alpha is the stifling rate. If a Spreader encounters another Spreader,
% Hibernator, or Stifler, the initiating Spreader becomes a Stifler at
% probability Alpha.

    
% Standard ISR Model
%     Kappa = .6;           % Degree of Network Reach                      
% 
%     Lamda = .6;           % Rate of Ign-to-Spr          
% 
%     Alpha = .3;      	    % Stifling rate

%% Variables
global Alpha
global Kappa
global Lambda
                                                 
    %% x Dump
    I = x(1);
    S = x(2);
    R = x(3); 
    
    %% ODE System
    dI = -Lambda*Kappa*I*S;
    dS =  Lambda*Kappa*I*S - Alpha*Kappa*S*(S+R);
    dR =  Alpha*Kappa*S*(S+R);
    

    %% Format output
    output = [dI; dS; dR];
    
end



%%Rumor Script
clear all; close all;
co = [0         0.4470    0.7410;
      0.6350    0.0780    0.1840;
      0.9290    0.6940    0.1250;
      0.4660    0.6740    0.1880;
      0.4940    0.1840    0.5560;
      0.8500    0.3250    0.0980;
      0.3010    0.7450    0.9330];
set(groot,'defaultAxesColorOrder',co)
%% For Loops
Storage = zeros(1001,3,5);
Storage2 = zeros(1001,3,5);
Lambda = 0:.2:.8;

for i = 1:length(Lambda)
   [Storage(:,:,i), Time] = SIR_Adapt(.6, Lambda(i), .3, .99, .01, 0, 0);
end

for i = 1:length(Lambda)
   [Storage2(:,:,i), Time] = SIR_Adapt(Lambda(i), .6, .3, .99, .01, 0, 0);
end

%% Variable Lambda Spreaders
Fig1 = figure;
for i = 2:length(Lambda)
plot(Time,Storage(:,2,i));
hold on;
end
%savefig(Fig1,'Variable_Lambda_Spreaders')


%% Variable Lambda Ignorants 
Fig2 = figure;
for i = 2:length(Lambda)
plot(Time,Storage(:,1,i));
hold on;
end
savefig(Fig2,'Variable_Lambda_Ignorants')

%% Variable Lambda Stiflers
Fig3 = figure;
for i = 2:length(Lambda)
plot(Time,Storage(:,3,i));
hold on;
end
savefig(Fig3,'Variable_Lambda_Stiflers')

%% Variable Kappa on Stiflers
Fig4 = figure;
for i = 2:length(Lambda)
plot(Time,Storage(:,3,i));
hold on;
end
savefig(Fig4,'Variable_Kappa_Stiflers')


% %SIR_Adapt(Kappa_input, Lambda_input, Alpha_input, Susceptibles, Spreaders, Stiflers,Plot_Option)
x = 0:.1:100;

%% Normal ISR
normal = SIR_Adapt(.6, .6, .3, .99, .01, 0, 0); % L = .6,   A =.3

figure(1),
plot(x, normal)
title('ISR Model of Rumor Spreading')
xlabel('Time')
ylabel('Density of Classes')
legend('Ignorant', 'Spreader', 'Stifler');

%% Alpha = Lambda
about_same = SIR_Adapt(.6, .5, .49, .99, .01, 0, 0); % L = .5,   A =.49

figure(2),
plot(x, about_same);
title('\lambda = .5, \alpha = .49, \kappa =.6')
xlabel('Time')
ylabel('Density of Classes')
legend('Ignorant', 'Spreader', 'Stifler');

%% Lambda >>> Alpha
lambda_over_alpha = SIR_Adapt(.6, .95, .05, .99, .01, 0, 0); % L = .95,   A =.05

figure(3),
plot(x, lambda_over_alpha);
title('\lambda = .95, \alpha = .05, \kappa =.6')
xlabel('Time')
ylabel('Density of Classes')
legend('Ignorant', 'Spreader', 'Stifler');

%% Alpha >>> Lambda
alpha_over_lambda = SIR_Adapt(.6, .05, .95, .99, .01, 0, 0); % L = .05,   A =.95

figure(4),
plot(x, alpha_over_lambda);
title('\lambda = .05, \alpha = .95, \kappa =.6')
xlabel('Time')
ylabel('Density of Classes')
legend('Ignorant', 'Spreader', 'Stifler');

%% Normal ISR with almost all stiflers at end
lambda8 = SIR_Adapt(.6, .8, .3, .99, .01, 0, 0); % L = .05,   A =.95

figure(5),
plot(x, lambda8);
title('\lambda = .8, \alpha = .3, \kappa =.6')
xlabel('Time')
ylabel('Density of Classes')
legend('Ignorant', 'Spreader', 'Stifler');

