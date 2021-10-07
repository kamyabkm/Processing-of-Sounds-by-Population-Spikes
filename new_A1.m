
function A1
N_E = 3;
N_I = 2;
P = 5;

global J_ee0 J_ei0 J_ii0 J_ie0 J_ee1 J_ee2 J_ie1 J_ie2;

J_ii0 = -0.5;
J_ie0 = 0.5;
s=0; %sensory 
U = 0.5;
tau = 0.001; % tau constant
tau_ref = 0.003;
tau_rec = 0.8; 


%% background synaptic inputs
% excitatory
ee_test = [];
ei_test = [];
for i = 1:P
    e_temp = rand(N_E,1); % initial random sampling of N_E
    e_temp_min = min(e_temp);
    e_temp_max = max(e_temp);
    e_min = -10; % parameters
    e_max = 10; % parameters
    e = ((e_max-e_min)*(e_temp-e_temp_min)/(e_temp_max-e_temp_min))+e_min; % scaling to e_min + e_max 
    ee_final = sort(e); % uniform distribution
    ee_test = [ee_test,ee_final];

% inhibitory
    e_temp = rand(N_I,1); % initial random sampling of N_I
    e_temp_min = min(e_temp);
    e_temp_max = max(e_temp);
    e_min = -10;
    e_max = 10;
    e = ((e_max-e_min)*(e_temp-e_temp_min)/(e_temp_max-e_temp_min))+e_min; % scaling
    ei_final = sort(e);
    ei_test = [ei_test,ei_final];
end
    

lambdaX = 0.25;
alpha = 2;
delta_left = 5;
delta_right = 5;

vs= [];
for i=1:P
    vs0 =[rand(N_E,1);rand(N_I,1);rand(N_E,1);rand(N_I,1)]; % order is [E,I,x,y]
    vs = [vs; vs0];
end

% rate_auditory(0,vs)
% tspan = [-10 2];
% [tt,xx] = ode45(@rate_auditory,tspan,vs);

E_range = N_E;
I_range = N_E+N_I;
x_range = 2*N_E+N_I;
% OE = xx(1:length(tt),1:E_range*P);
% OI = xx(1:length(tt),E_range*P+1:I_range*P);
% Ox = xx(1:length(tt),I_range*P+1: x_range*P);
% Oy = xx(1:length(tt),x_range*P+1:end);
% 
% mOE = mean(OE,2);
% mOI = mean(OI,2);
% close all;
% figure;
% it = floor(length(tt)/2);
% plot(mOE,mOI);
% hold on
% title('phase-plane diagram')
% 
% figure;
% subplot(2,2,1); plot(tt,OE);
% hold on; plot(tt,mOE,'linewidth',2);
% subplot(2,2,2); plot(tt,OI);
% subplot(2,2,3); plot(tt,Ox);
% subplot(2,2,4); plot(tt,Oy);
% %% nested function
function out = rate_auditory(t,vs)
    
    % state variables in matrix and vector form

    E = vs(1:E_range*P);
    E_mat = reshape(E,N_E,P);
    I = vs(E_range*P+1:I_range*P);
    I_mat = reshape(I, N_I, P);
    x = vs(I_range*P+1:x_range*P);
    x_mat = reshape(x, N_E,P);
    y = vs(x_range*P+1:end);
    y_mat = reshape(y, N_I,P);
%     E_mat = vs(1:N_E,1:P);
%     I_mat = vs(N_E+1:N_E+N_I,1:P);
%     x_mat = vs(N_E+N_I+1: 2*N_E+N_I,1:P);
%     y_mat = vs(2*N_E+N_I+1:end,1:P);

    sum1_E=[];
    sum1_I= [];
    for q=1:P
        switch q
            case 1
                R_range = 0:2;
            case 2
                R_range = -1:2;
            case P-1
                R_range = -2:1;
            case P
                R_range = -2:0;
            otherwise
                R_range = -2:2;
        end
        q_sumE=[];
        q_sumI = [];
        for R = R_range
            var1 = j_ee(abs(R))/N_E;
            var2 = sum(U*x_mat(:,q+R).*E_mat(:,q+R));
            final = var1.*var2;
            q_sumE = [q_sumE,final];
            sum_e = sum(E_mat(:,q+R)); %sum of E for the I rate
            sum_I = (j_ie(abs(R))/N_E) * sum_e;
            q_sumI = [q_sumI,sum_I];
        end
        sum1_E = [sum1_E,sum(q_sumE)];
        sum1_I = [sum1_I, sum(q_sumI)];
        
    end
     
    sum2_E = (J_ei0/N_I) * sum(U.*y_mat.*I_mat);
    %s=0 BUT double check with markus
    
    out_E = max(0,sum1_E + sum2_E + ee_test); %relu
    

    sum_I = sum1_I +J_ii0/N_I * sum(I_mat);    
    out_I = max(0,sum_I+ei_final); %relu

    dEdt = (-E_mat + (1-tau_ref*E_mat).*out_E)/tau;
    dIdt = (-I_mat + (1-tau_ref*I_mat).*out_I)/tau;
    dxdt = (1-x_mat)/tau_rec - U*x_mat.*E_mat;
    dydt = (1-y_mat)/tau_rec - U*y_mat.*I_mat;
    

%     out = [dEdt; dIdt; dxdt; dydt]; % single column
    out= [dEdt(:);dIdt(:);dxdt(:);dydt(:)];
end



function out = j_ee(R)
        if R == 0
            out = 6;
        elseif R == 1
            out = 0.045;
        else
            out = 0.015;
        end
end
function out = j_ie(R)
    if R == 0
        out = 6;
    elseif R == 1
        out = 0.0035;
    else
        out = 0.0015;
end
function out_h = heavf(x) %heaviside function    
    for i = 1:length(x)
        if x(i)>0
            out_h(i)=1;
        else
            out_h(i)=0;
        end
    end      
end
end
end



