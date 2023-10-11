clear all
LineWidth = 3; MarkerSize = 8; FontSize = 14;

T = 0.005;         
f0 = 1/T;       
NoOfPointsPerPeriod = 50;
A = 0.2;          
n_max = 5;     
ns = 1:1:n_max; 

t = -T:T/NoOfPointsPerPeriod:T;     
s = (A*square(2*pi/T*(t+T/4))); 

A0 = 0;
An = 4*A./(ns*pi) .*sin(ns*pi/2); 
Bn = zeros(1,size(ns,2));


s_e = A0/2*ones(1,size(t,2)); 

for n=1:1:n_max
   s_e = s_e + An(n)* cos(2*pi*n/T*t) + Bn(n)* sin(2*pi*n/T*t);
end

figure(1)
plot(t, s,'b-', t, s_e,'r--', 'LineWidth', LineWidth);
ylabel('s(t) - volts'); xlabel('time - seconds');

legend('original s(t)', 'up to n = 5');
grid on; set(gca, 'FontSize', FontSize);

P_AC = (An.^2+Bn.^2)/2; 
P_all = [(A0/2)^2 P_AC]; 
n_all = [0 ns];        
freq  = n_all*f0;     

figure(2); 
stem(freq, P_all, 'LineWidth', 2, 'MarkerSize', 8); 
ylabel('PSD (Watts/Hz)'); xlabel('frequency - Hz');
grid on; set(gca, 'FontSize', FontSize);
