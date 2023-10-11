% Ahmad Fahad Alzhrani 201917030
% Programming Assignment #2 - Semester 211
k = 1.38e-23; % Boltzman's constant
Thetas = [1/4 3/4 -3/4 -1/4]*pi; % M-phases
T = 40e3; % Effective temperature in Kelvin
A = 0.4e-6; % symbol amplitude in volts
fc = 20e3; % Carrier frequency in Hz
Ts = 0.2e-3; %Symbol duration in sec
r = 1.0; % filter efficiency number
% sequence of data bits
Seq = [1 0 1 0 0 0 1 0 0 1 1 1];

Es = A^2 /2 * Ts;    % the needed values in the equations 
N0 = k * T;
EbN0 = (Es/2)/N0;
Rs = 1/Ts;
Rb = 2*Rs;
BER = qfunc(sqrt(2*EbN0));  
BW = (1+r)*Rs;
BWe = Rb/BW;
signal_pwr = A^2/2;
noise_pwr = N0 * BW;
SNR = (Es/Ts)/(noise_pwr);
bitNum = size(Seq);


symb = [0];        % the vector for the symbols
data = ["0"];       % the list for the data
iq = ["0"];          % the list for the quarters
phase = ["0"];        % the list for the phase    

for i = 1:2:bitNum(2)             % loop to iterate in the sequence and iclude the values ..
    one = Seq(i);                   % in its specific list or vector 
    two = Seq(i+1);

    if one==0 && two==0                   % if statement to divide the value into four types
        symb(end+1)= 0;                     
        data(end+1)="0,0";
        iq(end+1)= "-1,-1";
        phase(end+1)= "-3/4";
    elseif one==0 && two==1
        symb(end+1)=1;
        data(end+1)="0,1";
        iq(end+1)= "-1,1";
        phase(end+1)= "3/4";
    elseif one==1 && two==0
        symb(end+1)= 2;
        data(end+1)="1,0";
        iq(end+1)= "1,-1";
        phase(end+1)= "-1/4";
    else
        symb(end+1)= 3;
        data(end+1)="1,1";
        iq(end+1)= "1,1";
        phase(end+1)= "1/4";
    end
end
change = [0];           % list for the 1 pi changes in the phases
phasec = [0];               % list for the change in the phases
for i =(bitNum(2)/2):-1:2    
    i0 = str2num(phase(i+1));
    i1 = str2num(phase(i));   % loop to include the phase changes in the two lists above
    phasec(end+1)= i0 - i1;
    if i0 - i1 == 1
        change(end+1)=i;
    end
end

fprintf('    Es = %.2d joules, N0 = %.2d Watts/Hz \n', Es, N0);
fprintf('    Eb/N0 = %.3f dB, BER = %.2d \n', 10*log10(EbN0), BER);         %printing the results
fprintf('    Rs = %.2f sym/s, Rb = %.2f b/s \n', Rs, Rb);
fprintf('    BW = %.2f Hz, BWe = %.3f b/s/Hz \n', BW, BWe);
fprintf(' Signal pwr = %.4d Watts or (%.3f dBW) \n', signal_pwr, 10*log10(signal_pwr));
fprintf('  Noise pwr = %.4d Watts or (%.3f dBW) \n', noise_pwr, 10*log10(noise_pwr));
fprintf('   SNR = %.3f or (%.3f dB) \n', SNR, 10*log10(SNR));
s = size(change);
fprintf(' There are %.f pi changes ',s(2)-1);
if s(2) == 1
    
else
    fprintf('- occuring at');
end
    
for i = s(2):-1:2
    fprintf(' %.f ',change(i));
end
fprintf('\n');
fprintf(' Number of QPSK symbols Ns = %.f, Number of bits = %.f \n',bitNum(2)/2,bitNum(2));
fprintf(' Detailed Table - I and Q signals \n');
fprintf(' Index: ');
for s = 1:bitNum(2)/2
    fprintf(' |     %.f   ',s);
end
fprintf('\n');
fprintf(' Symbol:');
for s = 2:(bitNum(2)/2+1)
    fprintf(' |     %.f   ',symb(s));
end
fprintf('\n');
fprintf(' Data:  ');                        % loops for the values in the lists
for s = 2:(bitNum(2)/2+1)
    fprintf(' |   %s   ',data(s));
end
fprintf('\n');
fprintf(' (I,Q):');
for s = 2:(bitNum(2)/2+1)
    fprintf('  | (%5s)',iq(s));
end
fprintf('\n');
fprintf(' Phase:  ');
for s = 2:(bitNum(2)/2+1)
    fprintf('|%6s pi ',phase(s));
end
fprintf('\n');
fprintf(' PhaseC: |         ');
for s = (bitNum(2)/2):-1:2
    fprintf(' |  %4.1f pi',phasec(s));
end
fprintf('\n');
