%% Soft Acoustic Modem Receiver 
% Adam Hess


%% Transmitter Side
% Define variables
clear all;
close all
clc;
A = 1; %absolute amplitude of NRZ value 
Rb = 100; % Bit rate
Tb = 1/Rb; %period for 100bits/sec
fs = 44000; %sampling frequency, also used for number of points transmitted per second;
PPB = fs/Rb;%number of Points per bit;
msgSize = 25; %characters
startBits = 200; %must be multiple of 2
zeroBits=startBits; %must be multiple of two
% prepare the 2 streams of binary sequences for differential encoding
%% Collect Users Message
message = input('Please input your message you wish to transmit. (Max 15 chars)\n You must transmit at least one alphabetical character\n', 's');
%FIXED SIZED MESSAGE
%IF MESSAGE IS TOO SMALL THEN CONCATINATE SPACES TO IT
if (length(message) < (msgSize+1))
    tmplength = length(message);
    for z = 1:(msgSize-tmplength)
        message(tmplength + z) = ' ';
    end
end

message( msgSize + 1) = A; %elminiates 6 bit entries 
message = abs(message); %converts to ascii (remember to add a zero at the beginig of every bit
message = dec2bin(message);

%% Linerizer
%NxM matrix and turns it into a Z long array that is a multiple of eight by
%appending a zero to the start of each seven bit message
linearize = 1;
message_lin = [];
for i = 1:msgSize
    message_lin =cat(2, message_lin, '0', message(i, 1:7));
end

%message_lin
%% Convolutional Encouder
for i = 1:length(message_lin)
    message_num(i) = str2num(message_lin(i));
end

figure(1)
subplot(321)
plot(message_num);
title('Original message')

clear message_lin;
%zero initial conditions
prev1 = 0;
prev2 = 0;
count_conv = 1;
for i = 1:length(message_num)
    convE(count_conv) = prev1 + prev2 + message_num(i);
    count_conv = count_conv + 1;
    convE(count_conv) = prev2 + message_num(i);
    count_conv = count_conv + 1;
    prev2 = prev1;
    prev1 = message_num(i);
end

convE = mod(convE,2);

subplot(322)
plot(convE)
title('Codewords from convolutional encoder')

%% Frame Formating
% concatinates startBit, barkersequence and endingsequence

for i = 1:startBits
    formatOnes(i) = 1;
end
formatZeros = 0;

for i = 1:zeroBits
    formatZeros(i) = 0;
end
frameFormat = 0;

barkerSeq = [ 1 1 1 1 1 0 0 1 1 0 1 0 1];
frameFormat = cat(2, formatOnes, formatZeros, barkerSeq);
endFrameFormat = [ 1 1 1 1 1 1 1 1 1 ]; %ending sequence 1111_1111_0101_0101 the 
message_format = cat(2, frameFormat, convE, endFrameFormat);

subplot(323)
plot(message_format);
title('formatted message')

%% Parser I
%Seperates the two bits into repsective unencoded I and Q channels
% The I channel contanins all Odd numbered bits (1, 3, 5)
% the Q channel contains all even Numbered bits (2, 4, 6)
I_parse = 1;
Q_parse = 1;
count_real = 1;
count_imag =1;
for i = 1:length(message_format)
    if (mod(i,2) == 1)
        I_parse(count_real) = message_format(i);
        count_real = count_real + 1;
    else
        Q_parse(count_imag) = message_format(i);
        count_imag = count_imag + 1;
    end
end



%% NRZ Level Generator
% converts message from a string to  +/- 1 value (non-return to zero)

I_NRZ = 0;
Q_NRZ = 0;
for i = 1:length(I_parse)
            if (I_parse(i) == 1)
                I_NRZ(i)=1;
                        
        elseif (I_parse(i) == 0)
            I_NRZ(i)=-1;             
        else
            error('INVALID MESSAGE ENTRY');
        end
end
for i = 1:length(Q_parse)
            if (Q_parse(i) == 1)
                Q_NRZ(i)=1;
                        
        elseif (Q_parse(i) == 0)
            Q_NRZ(i)=-1;             
        else
            error('INVALID MESSAGE ENTRY');
        end
end

subplot(324)
plot(I_NRZ);
title('I Channel NRZ');
subplot(325);
plot(Q_NRZ);
title('Q Channel NRZ')

%% DIFFERENTIAL ENCODER BLOCK
I=I_NRZ; %binary sequence u1
Q=Q_NRZ;  %binary sequence u2
s_k=I+j*Q;
s_k=s_k/sqrt(2); %normalize S_k by square root of 2
s_k; %input to the diff. encoder
x_k=s_k(1);
for t=2:length(s_k)
    x_k(t)=x_k(t-1)*s_k(t);
end
x_k; %output of diff. encoder

%separate the real and imaginary terms of X_k for IQ decompositon
I_d=real(x_k);
Q_d=imag(x_k);




%% START_FRAME FORMATTING BLOCK
% taken care of already. Need to differentially encode the ENTIRE signal
% not just the message bits. Need to search for barker sequence

%% Waveform Generator
I_wave = [];
Q_wave = [];
for i=1:length(I_d)
    I_wave=cat(2,I_wave,I_d(i));
    I_wave=cat(2,I_wave,zeros(1,PPB-1));   
end
for k=1:length(Q_d)
    Q_wave=cat(2,Q_wave,Q_d(k));
    Q_wave=cat(2,Q_wave,zeros(1,PPB-1));
end

figure(2)
subplot(321)
plot(I_wave);
title('I channel oversampling');
subplot(322);
plot(Q_wave);
title('Q channel oversampling');


%% RRF ROOT-RAISED-COSINE FILTER BLOCK
t_T=(-3:Rb/fs:3)+1e-8;
r = 0.35; %roll off factor;
RRCF =(cos((1+r)*pi*t_T)+sin((1-r)*pi*t_T)./(4*r*t_T))./(1-(4*r*t_T).^2);
RRCF=RRCF/max(RRCF);

%Convolve with NRZ_I to get message signal for modulator
I_message =  conv(RRCF, I_wave); %I channel

subplot(323)
plot(I_message);
title('I Channel Filter Output');

%Convolve with NRZ_Q to get message signal for modulator
Q_message =  conv(RRCF, Q_wave); %Q channel

subplot(324)
plot(Q_message);
title('Q Channel Filter Output');

%clear I_wave Q_wave;


%% BANDPASS MODULATOR BLOCK
totalmessage= Q_message + I_message; %need total length of the block, i.e., the longer one is the value
t_c = 0:1/fs:(length(totalmessage)/fs - 1/fs);
fc = 1000; %carrier frequency
w = 2*pi*fc;

S_transmit = Q_message.*sin(w*t_c) + I_message.*cos(w*t_c);
%transmitted message for QPSK is the sum of the I and Q channels multiplied
%by cosine and sine respectivly

subplot(325)
t_transmit = 0:1/fs:(length(S_transmit)-1)/fs;
plot(t_transmit,S_transmit);
title('transmitted signal');

%% MATLAB AUDIO OUTPUT


S_transmit = S_transmit';
disp('Now, get the Receiver ready.');
disp('Since the Receiver has delay, you need to trigger the Receiver first.');
disp('Press any key RIGHT AFTER trigger the Receiver To Transmit');
pause;
wavplay(S_transmit, fs, 'async');
