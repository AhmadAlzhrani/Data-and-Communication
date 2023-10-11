clear all
LineWidth = 3; MarkerSize = 8; FontSize = 14;
% Section0 {Clear all; initialization of A, T, n_max and other constants}
% define constants - These are PROBLEM SPECIFIC
T = 0.1;          % in sec 
f0 = 1/T;       % in Hz
NoOfPointsPerPeriod = 50; % needed to define t-axis
A = 2;          % amplitude in Volts
fer = 100;
n_max = fer*T;     % max n index in the FSE summation
ns = 1:1:n_max; % define n vector from 1 to n_max

% Section1 {define t-axis and computation of original s(t)}
% define original s(t) using the periodic Matlab square() function
t = -T:T/NoOfPointsPerPeriod:T;     % define the t-axis
s = (A*square(2*pi/T*t)); % This is PROBLEM SPECIFIC


% Section2 {computation of A0, An, Bn}
% Define FSE constants A0, An, and Bn - These are PROBLEM SPECIFIC
A0 = 0;
An = zeros(1,size(ns,2)); % this is valid for all ns
Bn = 4*A./(pi*ns);
Bn(2:2:n_max) = 0;

% Section3 {computation of FSE using for loop}
% Build s(t) using FSE summation - This is NOT PROBLEM SPECIFIC
s_e = [0]; % Initially the DC component
% Now add the AC components till n_max
for n=1:1:n_max
   s_e = s_e + An(n)* cos(2*pi*n*f0*t) + Bn(n)* sin(2*pi*n*t*f0);
end

% Section4 {plot of original s(t) and FSE expansion}
% plot original s(t) with FSE - This is NOT PROBLEM SPECIFIC
figure(1)
plot(t, s_e,'r-', 'LineWidth', LineWidth);
ylabel('s(t) - volts'); xlabel('time - seconds');
%axis([0 1 -0.2 1.2]); % this line is problem specific
legend('original s(t)');
grid on; set(gca, 'FontSize', FontSize);


% Section5 {PSD function computation}
% To plot the PSD function - This is NOT PROBLEM SPECIFIC
P_AC = (Bn.^2)/2; % power of nth harmonic
P_all = [(A0/2)^2 P_AC]; % DC power added as first element of vector
n_all = [0 ns];         % index for DC power is zero - added to the vector
freq  = n_all*f0;       % define frequency axis

% Section6 {plotting PSD function using stem()}
figure(2); % This is NOT PROBLEM SPECIFIC
stem(freq, P_all, 'LineWidth', 2, 'MarkerSize', 8); 
%axis([-1 12 0 0.3]); % this line is problem specific
ylabel('PSD (Watts/Hz)'); xlabel('frequency - Hz');
grid on; set(gca, 'FontSize', FontSize);

totalPower = A^2; % same as a regular cosine
powerFSE = sum(P_AC); % from from the FSE
dynamicP = P_AC./max(P_AC); % normalized to the maximum
PwrRatio = 10*log10(dynamicP);

fprintf('Signal amplitude, A = %7.3f volts\n', A);
fprintf(' Signal period, T = %7.3f sec - fundamental freq, f0 = %7.3f Hz\n',T, f0);
fprintf('Filter cut-off, X = %7.3f Hz - Highest component, n = %3d (f = %7.3f Hz)\n',fer,  n_max-1, fer-10);
fprintf('power s(t)= %7.3f Watts, power s0(t)= %7.3f Watts (or %0.3f %%)\n',totalPower,powerFSE,powerFSE/totalPower*100);
fprintf('Max power at f = %4.1f Hz - Power = %7.3f (%7.3f dBW)\n',freq(2), P_AC(1), 10*log10(P_AC(1)));
fprintf('Dynamic range = %7.3f dB\n', -1*PwrRatio(length(PwrRatio)-1));

fprintf('\n\nPower spectral density function and power ratio function:\n');
nn = [1:2:n_max];
fprintf('Index :'); fprintf(' %6d ',nn );fprintf('\n');
fprintf('freq (Hz) :'); fprintf(' %6.1f ', nn*f0); fprintf('\n');
fprintf('PSD (W) :'); fprintf(' %6.3f ', P_AC(1:2:n_max)); fprintf('\n');
fprintf('PSD-normalized :'); fprintf(' %6.3f ', dynamicP(1:2:n_max)); fprintf('\n');
fprintf('PSD-normalized (dB):'); fprintf(' %6.2f ', PwrRatio(1:2:n_max)); fprintf('\n');

figure(3); % This is NOT PROBLEM SPECIFIC
stem(freq(2:1:n_max+1),PwrRatio, 'LineWidth', 2, 'MarkerSize', 8); 
%axis([-1 12 0 0.3]); % this line is problem specific
ylabel('PSD (Watts/Hz)'); xlabel('frequency - Hz');
grid on; set(gca, 'FontSize', FontSize);
