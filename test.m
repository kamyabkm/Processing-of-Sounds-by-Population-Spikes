
function test
N_E = 2;
N_I = 2;
U=0.5;
P = 5;

vs= [];
for i=1:P
    vs0 =[rand(N_E,1);rand(N_I,1);rand(N_E,1);rand(N_I,1)]; % order is [E,I,x,y]
    vs = [vs, vs0];
end

E = vs(1:N_E,1:P);
I = vs(N_E+1:N_E+N_I,1:P);
x = vs(N_E+N_I+1: 2*N_E+N_I,1:P);
y = vs(2*N_E+N_I+1:end,1:P);

for q = 3:P-2
        for R = -2:2
            var1 = j_ee(abs(R))/N_E;
            var2 = U.*x(1:N_E,1:q+R).*E(1:N_E,1:q+R);
          
%             sum1_temp = sum(get_j(abs(R))/N_E)
            sum1 = [sum1,var1];
        end
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
end
