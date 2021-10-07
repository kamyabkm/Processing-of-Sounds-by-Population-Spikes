

function single_column

N_E = 5;
N_I = 5;
J_ee = 6; %synaptic efficacy for exc to exc 
J_ei = -4; %synaptic efficacy for inh to exc 
J_ie = 0.5; %synaptic efficacy for exc to inh 
J_ii = -0.5; %synaptic efficacy for inh to inh 
s=0; %sensory 
U = 0.5;
tau = 0.001; % tau constant
tau_ref = 0.003;
tau_rec = 0.8; 


%% background synaptic inputs
% excitatory
e_temp = rand(N_E,1); % initial random sampling of N_E
e_temp_min = min(e_temp);
e_temp_max = max(e_temp);
e_min = -10; % parameters
e_max = 10; % parameters
e = ((e_max-e_min)*(e_temp-e_temp_min)/(e_temp_max-e_temp_min))+e_min; % scaling to e_min + e_max 
ee_final = sort(e); % uniform distribution

% inhibitory
e_temp = rand(N_I,1); % initial random sampling of N_I
e_temp_min = min(e_temp);
e_temp_max = max(e_temp);
e_min = -10;
e_max = 10;
e = ((e_max-e_min)*(e_temp-e_temp_min)/(e_temp_max-e_temp_min))+e_min; % scaling
ei_final = sort(e);
    
    
vs0 =[5*ones(N_E,1);5*ones(N_I,1);0.5*ones(N_E,1);0.5*ones(N_I,1)]; % order is [E,I,x,y]

final= rate_auditory(0,vs0);
size(final)

tspan = [-10 2];
[tt,xx] = ode45(@rate_auditory,tspan,vs0);

OE = xx(1:length(tt),1:N_E);
OI = xx(1:length(tt),N_E+1:N_E+N_I);
Ox = xx(1:length(tt),N_E+N_I+1: 2*N_E+N_I);
Oy = xx(1:length(tt),2*N_E+N_I+1:end);

mOE = mean(OE,2);
mOI = mean(OI,2);
close all;
figure;
it = floor(length(tt)/2);
plot(mOE,mOI);
hold on
title('phase-plane diagram')

figure;
subplot(2,2,1); plot(tt,OE);
hold on; plot(tt,mOE,'linewidth',2);
subplot(2,2,2); plot(tt,OI);
subplot(2,2,3); plot(tt,Ox);
subplot(2,2,4); plot(tt,Oy);
% plot(tt,xx);


    
%% nested function
function out = rate_auditory(t,vs)
    E = vs(1:N_E);
    I = vs(N_E+1:N_E+N_I);
    x = vs(N_E+N_I+1: 2*N_E+N_I);
    y = vs(2*N_E+N_I+1:end);

    if t > 0
        s = 4;
    end

    sum_E = J_ee/N_E * sum(U.*x.*E) + (J_ei/N_I * (sum(U.*y.*I)));
    out_E = max(0,sum_E + ee_final + s); %relu
    sum_I= J_ie/N_E * sum(E) + J_ii/N_I * sum(I);
    out_I = max(0,sum_I+ei_final); %relu

    dEdt = (-E + (1-tau_ref*E).*out_E)/tau;
    dIdt = (-I + (1-tau_ref*I).*out_I)/tau;
    dxdt = (1-x)/tau_rec - U*x.*E;
    dydt = (1-y)/tau_rec - U*y.*I;

    out = [dEdt; dIdt; dxdt; dydt]; % single column 
    end

end