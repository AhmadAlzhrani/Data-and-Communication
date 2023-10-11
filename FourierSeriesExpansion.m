clear all
FontSize = 14; LineWidth = 2; MarkerSize = 8;
FontSizeLegend = 12;

%
% full-wave rectified cosine
A = 2; % volts
T = 0.1; % second
f0 = 1/T; % fundamental freq in Hz
fer = 40;
num = T*fer;

n = 1:num; % define my n vector
An = zeros(size(n)); % the general An term - needs fixing
Bn = 4*A./(pi*n); % Bn is all zeros for this example
Bn(2:2:num) = 0; % 0 for odd n's

%
% The DC component - from the formula is C0/2
A0 = 0;

DC_Power = (A0/2)^2;
AC_Power = (An.^2+Bn.^2)/2;

Total_Power_Signal = A^2; % same as a regular cosine
Total_Power_FSE = DC_Power + sum(AC_Power); % from from the FSE

%
% How to plot signal components
% want to evaluate s_e(k=0), s_e(k=2), s_e(k=4), and s_e(k=6)
t = -T:T/100:T; % define your time axis
s_e = [0]; % Initially the DC component
% Now add the AC components till n_max
for n=1:1:num
   s_e = s_e + An(n)* cos(2*pi*n*f0*t) + Bn(n)* sin(2*pi*n*t*f0);
end

%
%
% original s_t
s_t = A*sin(2*pi*t/T);

%
% for part (e) building the PSD function and the power ratio function
X = 40; ns_max_9 = X/10; % for the filtered signal
ns_9 = [1:2:ns_max_9]; % multiples of f0 with non-zero powers
freq_9 = ns_9*f0;
temp = [DC_Power AC_Power(1:9)]; % contains 0 powers at 1f0, 3f0, & 5f0 - remove
PSD_9 = temp(1:2:10); % keep only 0, 2f0, 4f0, and 6f0 - no 0 power components
[MaxP, MaxI] = max(PSD_9); % identify the maxmuim
PSDn_9 = PSD_9./MaxP; % normalized to the maximum
PwrRatio = 10*log10(PSDn_9); % in dBs
TotalPower_9 = sum(PSD_9);

fprintf('Signal amplitude, A = %7.3f volts\n', A);
fprintf(' Signal period, T = %7.3f sec - fundamental freq, f0 = %7.3f Hz\n', ...
T, f0);
fprintf('s(t) - Total power = %7.3f Watts or %7.3f dBW\n', Total_Power_Signal, 10*log10(Total_Power_Signal));
fprintf('Filter cut-off, X = %7.3f Hz - Highest component, n = %3d (f = %7.3f Hz)\n', ...
X, ns_max_9, max(freq_9));
fprintf('s_0(t) - Total power = %7.3f Watts - or %5.2f%% of total power\n', ...
TotalPower_9, TotalPower_9/Total_Power_Signal*100);
fprintf('s_0(t) - Max power at f = %4.1f Hz - Power = %7.3f (%7.3f dBW)\n', ...
freq_9(MaxI), MaxP, 10*log10(MaxP));
fprintf('s_0(t) - Dynamic range = %7.3f dB\n', -1*min(PwrRatio));
fprintf('\n\nPower spectral density function and power ratio function:\n');
fprintf('Index :'); fprintf(' %6d ', ns_9); fprintf('\n');
fprintf('freq (Hz) :'); fprintf(' %6.1f ', freq_9); fprintf('\n');
fprintf('PSD (W) :'); fprintf(' %6.3f ', PSD_9); fprintf('\n');
fprintf('PSD-normalized :'); fprintf(' %6.3f ', PSDn_9); fprintf('\n');
fprintf('PSD-normalized (dB):'); fprintf(' %6.2f ', 10*log10(PSDn_9)); fprintf('\n');

tt = t*1000; % in msec
figure(1); % (a) plot for s(t)
h1 = plot(tt, s_t); grid on;
ylabel('Amplitude (volts)'); xlabel('time - msec');
title('The original s(t)');
set(h1, 'LineWidth', LineWidth); set(gca, 'FontSize', FontSize);

figure(2); % (b) some FSE expansions
h2 = plot(tt,s_e);
grid on; ylabel('Amplitude - volts'); xlabel('time - msec');
title('The original s(t) and its FSE');
set(h2, 'LineWidth', LineWidth);
set(gca, 'FontSize', FontSize);
set(hl, 'FontSize', FontSizeLegend, 'Interpreter', 'none');

figure(3); % (c) plot of PSD function
h1 = stem([0:10]*f0, [DC_Power AC_Power(1:10)]); grid on;
ylabel('PSD function (Watts/Hz)'); xlabel('freq (Hz)');
set(h1, 'LineWidth', LineWidth, 'MarkerSize', MarkerSize);
set(gca, 'FontSize', FontSize);

figure(4); % (d) for LPF signal - same as s_e_6
h2 = plot(tt, s_e_5, '-k'); grid on;
title('LPF output s_O(t)'); xlabel('time - msec'); ylabel('Amplitude - volts');
set(h2, 'LineWidth', LineWidth);
set(gca, 'FontSize', FontSize);
figure(5); % (e) plot of dynamic
h1 = stem(freq_9, PwrRatio,'x'); grid on;
ylabel('dynamic power ratio (dB)'); xlabel('freq (Hz)');
set(h1, 'LineWidth', LineWidth, 'MarkerSize', MarkerSize);
set(gca, 'FontSize', FontSize);
