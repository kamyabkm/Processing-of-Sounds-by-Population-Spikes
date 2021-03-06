function A1
N_E = 100;
N_I = 100;
P = 15;

E_range = N_E;
I_range = N_E+N_I;
x_range = 2*N_E+N_I;

J_ei0 = -4;
J_ie0 = 0.5;
J_ii0 = -0.5;
J_ee0 = 6;
J_ee2 = 0.015;
J_ee1 = 0.045;
J_ie1 = 0.0035;
J_ie2 = 0.0015;

lambdaC = 0.25;
D_left = 5;
D_right = 5;

alpha = 2;
A = zeros(P,1);
A(8) = 12;
% A(3) = 4;


z = 0; %sensory
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

ee_test = repmat(ee_final,1,P);
ei_test = repmat(ei_final,1,P);

vs= [];
for i=1:P
    vs0 =[rand(N_E,1);rand(N_I,1);rand(N_E,1);rand(N_I,1)]; % order is [E,I,x,y]
    vs = [vs; vs0];
end

            
rate_auditory(0,vs)
tspan = [-5 0];
% run simulation to time zero
E_non_zero = true(N_E,P);
[tt,xx] = ode45(@rate_auditory,tspan,vs);

% run simulation from time zero
vs = xx(end,:)';
tspan = [0 0.8];
% E_non_zero = reshape(xx(end,1:N_E*P)>0,N_E,P); % <<<<<<<<<
[tt,xx] = ode45(@rate_auditory,tspan,vs);


OE = xx(1:length(tt),1:E_range*P);
OI = xx(1:length(tt),E_range*P+1:I_range*P);
Ox = xx(1:length(tt),I_range*P+1: x_range*P);
Oy = xx(1:length(tt),x_range*P+1:end);

E_range = N_E;
I_range = N_E + N_I;
x_range = 2*N_E + N_I;

I = vs(E_range*P+1:I_range*P);
        I_mat = reshape(I, N_I, P);
        x = vs(I_range*P+1:x_range*P);
        x_mat = reshape(x, N_E,P);
        y = vs(x_range*P+1:end);
        y_mat = reshape(y, N_I,P);
        
% time, columns, neurons
mOE = zeros(length(tt),P);
mOI = zeros (length(tt),P);
mOx = zeros(length(tt),P);
mOy = zeros(length(tt),P);

for i=1:P
    mOE(:,i) = mean(OE(:,(1:N_E)+(i-1)*N_E),2);
    mOI(:,i) = mean(OI(:,(1:N_I)+(i-1)*N_I),2);
    mOx(:,i) = mean(Ox(:,(1:N_E)+(i-1)*N_E),2);
    mOy(:,i) = mean(Oy(:,(1:N_I)+(i-1)*N_I),2);
end

% close all;
figure(1);
it = floor(length(tt)/2);
plot(mOE,mOI);
hold on
title('phase-plane diagram')
% 
% figure; 
% plot(tt,mOE(:,1),'linewidth',1); hold on; 
% plot(tt,mOI(:,1),'--','linewidth',1); xlim([-0.05 0.3])
% 
% figure; 
% plot(tt,mOE(:,8),'linewidth',1); hold on; 
% plot(tt,mOI(:,8),'--','linewidth',1); xlim([-0.05 0.3])


% 
% figure(2);
% subplot(2,2,1); plot(tt,mOE(:,1),'linewidth',1); title('Excitatory'); hold on; 
% plot(tt,mOI(:,1),'--','linewidth',1);
% subplot(2,2,2); plot(tt,OI); title('Inhibitory')
% subplot(2,2,3); plot(tt,Ox); title('x')
% subplot(2,2,4); plot(tt,Oy); title('y')

figure(3);
for i = 1:P
%     subplot(3,5,i); 
    plot(tt,mOI(:,i)); hold on;
    xlim([0 0.1])
end

figure(4);
uimagesc(1:P, tt, mOE)
ylim([0 0.07])
set(gca,'ydir','normal')
% 
% OE_1 = OE(:,1:4);
% OI_1 = OI(:,1:4);
% MOE_1 = mean(OE_1,2);
% MOI_1 = mean(OI_1,2);
% 
% OE_8 = OE(:,36:40);
% OI_8 = OE(:,36:40);
% MOE_8 = mean(OE_8,2);
% MOI_8 = mean(OI_8,2);
% figure; plot(tt,MOE_1,'--'); hold on; plot(tt,MOI_1);xlim([0 1])
% figure; plot(tt,MOE_8, '--'); hold on; plot(tt,MOI_8); xlim([0 1])


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
                sum_I = (j_ie(abs(R))/N_I) * sum_e;
                q_sumI = [q_sumI,sum_I];
            end
            
            sum1_E = [sum1_E,sum(q_sumE)];
            sum1_I = [sum1_I, sum(q_sumI)];
            mid = round(P/2);
            
            if t > 0
                for mid:
                    z = 1;
                end
                for mid-1 | mid + 1:
                    z = 0.5;
                end
                for mid-2 | mid+2:
                    z = 0.2;
                end
            end
            
            
            h = spatial(q);
            s = z*h;
            sum3_E(q) = sum(s);
        end
        
         sum2_E = (J_ei0/N_I) * sum(U.*y_mat.*I_mat);
        %s=0 BUT double check with markus
        
        
        out_E = max(0,sum1_E + sum2_E + ee_test + sum3_E.*sum3_E); %relu
        
        
        sum_I = sum1_I +J_ii0/N_I * sum(I_mat);
        out_I = max(0,sum_I+ei_test); %relu
        
        dEdt = (-E_mat + (1-tau_ref*E_mat).*out_E)/tau;
        dIdt = (-I_mat + (1-tau_ref*I_mat).*out_I)/tau;
        dxdt = (1-x_mat)/tau_rec - U*x_mat.*E_mat;
        dydt = (1-y_mat)/tau_rec - U*y_mat.*I_mat;
        
        out= [dEdt(:);dIdt(:);dxdt(:);dydt(:)];
    end


    function out = j_ee(R)
        if R == 0
            out = J_ee0;
        elseif R == 1
            out = J_ee1;
        else
            out = J_ee2;
        end
    end

    function out = j_ie(R)
        if R == 0
            out = J_ie0; 
        elseif R == 1
            out = J_ie1;
        else
            out = J_ie2;
        end
    end

    function out = spatial(q)
        lambda_S_left = max(lambdaC,lambdaC+(A-alpha)/D_left); % lambda_S_left for each Amplitude
        lambda_S_right = max(lambdaC,lambdaC+(A-alpha)/D_right); % lambda_S_right for each Amplitude
        M = (1:P)'; % frequencies
        lambda_S = (q < M).*lambda_S_left + (q >= M).*lambda_S_right;
        out = A.*(exp(-abs(q-M)./lambda_S));
    end
end