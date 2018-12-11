%% Binary ASK, FSK, PSK simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upsampling a binary signal
bin_vec = [1 0 1 1 0 0 1]; % short test signal

len_vec = length(bin_vec)*100;
bin_wave = ones(1,len_vec);
%create stacks of 1000
for i = 1:length(bin_vec)
    if i == 1
        bin_wave(i:i+99) = bin_vec(i)*ones(1,100);
    else
        bin_wave((100*(i-1))+1:100*i) = bin_vec(i)*ones(1,100);
    end
end

bp = .001; %bit period, how many seconds it takes to transmit a single stack
t1 = bp/100:bp/100:len_vec*(bp/100);

figure(2)
clf(2)
plot(t1,bin_wave)

%% ASK Modulation
A1 = 5; A2 = 0; br = 1/bp; fc = br*10; %carrier wave params
t2 = bp/100:bp/100:bp; %stacks of 100 instead of 1000 for simulation
t3 = bp/100:bp/100:bp*length(bin_vec); %time vector of full upsampled signal
s = zeros(1,length(bin_vec)*100);
s1 = zeros(1,length(bin_vec)*100);
for i = 1:1:length(bin_vec)
    if bin_vec(i) == 1 
        y = A1*cos(2*pi*fc*t2); %if a 1, multiply by carrier wave, A1 = 5
    else
        y = A2*cos(2*pi*fc*t2);
    end
    if i == 1
        s(i:i+99) = y; %if a 0, multiply by 0, A2 = 0
    else
        s(((i-1)*100)+1:(i*100)) = y; %concatenate the waves
    end
end

%visualizing the base carrier wave and the modulated wave
figure(10)
clf(10)
hold on;
plot(t3,A1*cos(2*pi*fc*t3),'r--')
plot(t3,s, 'b')
legend('Carrier Wave', 'OOK Modulated Wave', 'Location', 'southoutside')
hold off;
%% ASK Demodulation
demod = zeros(1,length(bin_vec));
step = length(t2);
final_step = length(s);%defining step sizes for looping through stacks
iter = step:step:final_step;

for i = iter
    y = cos(2*pi*fc*t2); 
    z = y.*s((i-(step-1)):i);%multiply each stack by the carrier wave
    int = trapz(t2,z);%integrate
    int_round = round(2*int/bp); %eliminate bp and 2 terms to isolate amplitude
    if int_round>2.5 %threshold of (A1+A2)/2
        a = 1;
    else
        a = 0;
    end
    demod(i/100) = a;
end
disp('OOK Demod Test. 1=pass, 0=fail')
all(demod==bin_vec)
%% FSK Modulation
A = 5; br = 1/bp; f1 = br*8; f2 = br*2; % carrier wave params
t2 = bp/100:bp/100:bp; %time vector for stacks
t3 = bp/100:bp/100:bp*length(bin_vec); %time vector for all stacks
s = zeros(1,length(bin_vec)*100);
for i = 1:1:length(bin_vec)
    if bin_vec(i) == 1
        y = A*cos(2*pi*f1*t2); %different frequencies for 1s and 0s
    else
        y = A*cos(2*pi*f2*t2);
    end
    if i == 1
        s(i:i+99) = y; %put the modulated stacks together
    else
        s(((i-1)*100)+1:(i*100)) = y;
    end
end

figure(10)
clf(10)
hold on;
plot(t3,s)
plot(t3,bin_wave)
legend('FSK Modulated Wave', 'Binary Data')
hold off
%% FSK Demodulation
demod = zeros(1,length(bin_vec));
step = length(t2);
final_step = length(s);%step sizes to loop through stacks
iter = step:step:final_step;
for i = iter
    y1 = cos(2*pi*f1*t2); %2 possible carrier waves
    y2 = cos(2*pi*f2*t2);
    z1 = y1.*s((i-(step-1)):i);
    z2 = y2.*s((i-(step-1)):i); %multiply s by each
    int1 = trapz(t2,z1);
    int2 = trapz(t2,z2); %integrate
    int1_round = round(2*int1/bp); %isolate amplitude
    int2_round = round(2*int2/bp);
    if(int1_round>A/2); %threshold of A/2
        a = 1;
    elseif(int2_round>A/2);
        a = 0;
    end
    demod(i/100) = a;
end
disp('FSK Demod Test. 1=pass, 0=fail')
all(demod==bin_vec)
%% PSK Modulation
A = 5; br = 1/bp; fc = br*10; %carrier params
t2 = bp/100:bp/100:bp; %stack time vector
t3 = bp/100:bp/100:bp*length(bin_vec);%full wave time vector
s = zeros(1,length(bin_vec)*100);

for i = 1:1:length(bin_vec)
    if bin_vec(i) == 1
        y = A*sin(2*pi*fc*t2); %180 degree shift between a 1 and 0
    else
        y = A*cos((2*pi*fc*t2)+pi); %-A*sin(2*pi*fc*t)
    end
    if i == 1
        s(i:i+99) = y;
    else
        s(((i-1)*100)+1:(i*100)) = y;
    end
end

figure(10)
clf(10)
% visualize PSK and binary data
plot(t3,s)
hold on;
plot(t3,bin_wave)
title('Binary Data and PSK Modulated Wave')
xlabel('Time (s)'); ylabel('Magnitude')
legend('PSK Wave','Binary Data')

%% PSK Demodulation
demod = zeros(1,length(bin_vec));

step = length(t2);
final_step = length(s);
iter = step:step:final_step;
for i = iter
    y = sin(2*pi*fc*t2);
    z = y.*s((i-(step-1)):i); %multiply the stack by the carrier for a 1
    int = trapz(t2,z); %integrate and isolate amplitude
    int_round = round(2*int/bp);
    if(int_round>0); %if a 1 the integral will be positive
        a = 1;
    else
        a = 0; %if a 0 it will be negative
    end
    demod(i/100) = a;
end
disp('PSK Demod Test. 1=pass, 0=fail')
all(demod==bin_vec)
