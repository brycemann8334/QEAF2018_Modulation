%% OOK Implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Binary data Modulation
load msg.mat
len_vec = length(bin_vec)*10000;
bin_wave = ones(1,len_vec);
for i = 1:length(bin_vec) %create stacks of 10,000
    if i == 1
        bin_wave(i:i+9999) = bin_vec(i)*ones(1,10000);
    else
        bin_wave((10000*(i-1))+1:10000*i) = bin_vec(i)*ones(1,10000);
    end
end

t = linspace(0,length(bin_wave)-1,length(bin_wave))/48000; %t vector for
%modulated wave
Ac = 1;
fc =  1000; %carrier params
c = Ac * cos(2*pi*fc*t);
s = c.*bin_wave; %OOK Modulation

figure(2)
clf(2)
hold on
%visualize modulated wave to be sent
plot(t,bin_wave)
plot(t,s)
hold off
%% TRANSMIT THE SIGNAL ACOUSTICALLY
%% Load recieved signal
[y, fs] = audioread('text_1k_10000bit.wav');
figure(1)
clf(1)
%visualize raw data
plot(y);
title('Raw Recieved Signal')
xlabel('Signal index'); ylabel('Amplitude')
%% Bandpass Filter, OOK recieved signal

bpFilt = designfilt('bandpassfir','FilterOrder',1024, ...
         'CutoffFrequency1',990,'CutoffFrequency2',1010, ...
         'SampleRate',fs);
fvtool(bpFilt) %visualize the mag and phase response of our bandpass

out = filter(bpFilt,y); %filter data with bandpass filter

figure(1)
clf(1)
%visualize raw vs filtered data
subplot(2,1,1)
plot(data);
title('Raw Recieved Signal')
xlabel('Signal Index'); ylabel('Amplitude')
ylim([-.5 .5])
subplot(2,1,2)
plot(out);
title('Bandpass (990 Hz to 1010 Hz) Filtered Signal')
xlabel('Signal Index'); ylabel('Amplitude')
ylim([-.5 .5])

% Frequency domain anlysis
N = length(out);
signal1 = fft(out); % filtered
shifted_signal1 = fftshift(signal1);
signal2 = fft(y); % unfiltered
shifted_signal2 = fftshift(signal2);
freqs = linspace(-pi, pi-2/N*pi, N) + pi/N*mod(N,2);
freqs_hz = (freqs*fs)/(2*pi);

figure(2)
clf(2)
%visualize raw and filtered in the frequency domain (not interesting)
subplot(2,1,1)
plot(freqs_hz,abs(shifted_signal1))
subplot(2,1,2)
plot(freqs_hz,abs(shifted_signal2))
%% Demodulate, OOK Vis
fc = 1000;
fs = 48000;
t = linspace(0,length(out)-1,length(out))'/fs;
demod = (out.*cos(2*pi*fc*t-2)); %include phase shift term from bandpass

N = length(out);

%low pass filter
lpFilt = designfilt('lowpassfir','PassbandFrequency',1/48, ...
         'StopbandFrequency',1/24, 'PassbandRipple',0.5, ...
         'StopbandAttenuation',100,'DesignMethod','kaiserwin');
fvtool(lpFilt)

dataOut = filter(lpFilt,demod); %filter

signal1 = fft(demod); % unfiltered
shifted_signal1 = fftshift(signal1);  
signal2 = fft(dataOut); % filtered
shifted_signal2 = fftshift(signal2);
freqs = linspace(-pi, pi-2/N*pi, N) + pi/N*mod(N,2);
freqs_hz = (freqs*fs)/(2*pi);

figure(111)
clf(111)
%visualize frequency domains before and after filtering
hold on;
plot(freqs_hz,abs(shifted_signal1))
plot(freqs_hz,abs(shifted_signal2))
legend('unfiltered', 'filtered')
xlabel('Frequency (Hz)'); ylabel('Amplitude')
%% Demodulate, OOK, Integration
t2 = linspace(0,10000,10000)/fs; %time vector for each stack
out2 = out(9.563e4:9.563e4+293e4-1); %truncate the data to a reasonable range
demod = zeros(1,293); %length of original message (doesn't have to be preset)
step = 10000; %stack size
final_step = step*293;
iter = step:step:final_step;
for i = iter
    y = cos(2*pi*fc*t2-2); %multiply bandpass filtered stacks by the carrier wave
    %including the phase shift correction
    z = y'.*out2((i-(step-1)):i);
    int = sum(abs(trapz(t2,z))); %compute the integral of each stack and 
    %compare to the threshold 
    %this operation is mimicing the low pass filter we used to visualize
    %the demodulated signal above
    if int>.1 %say whether the stack is a 1 or 0
        a = 1;
    else
        a = 0;
    end
    demod(i/10000) = a; %rebuild the demodulated signal
end
%% Compare Demodulated data to the original message

disp('A 1 means they are equivalent, A 0 means at least one digit is wrong')
all(demod==bin_vec)


























