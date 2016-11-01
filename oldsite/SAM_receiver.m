%% Soft Acoustic Modem Receiver 
% Adam Hess
%  RECEIVER
%% Frequency Recovery
clear all
close all
clc
A = 10; %absolute amplitude of NRZ value 
Rb = 100; % Bit rate
Tb = 1/Rb; %period for 100bits/sec
fs = 44000; %sampling frequency, also used for number of points transmitted per second;
PPB = fs/Rb;%number of Points per bit;
startBits = 60;
numChar = 30; %need to have larger transmission than reception;
numBits = numChar*8 +startBits;
numSamples = fs*3;
msgSize = 25;

t_receive = 5;
disp('Press any key to Receive Message');
pause; 
Received_m = wavrecord(fs*t_receive, fs);


figure(1)
t_received = 0:1/fs:(length(Received_m)-1)/fs;
plot(t_received,Received_m);
title('Received Signal'),
axis([0 (length(Received_m)-1)/fs -1.2 1.2]);

%% Non-Coherent Demodulation block
% assume Fc remains stable in the system
fc = 1000;
w = 2*pi*fc;
t_rec = 0:1/fs:length(Received_m)/fs - 1/fs;
Received_m = Received_m';


I_r = Received_m.*cos(w*t_rec);
Q_r = Received_m.*sin(w*t_rec);

figure(2)
subplot(321)
plot(I_r)
title('I channel Passband Demodulation')
subplot(322)
plot(Q_r)
title('Q channel Passband Demodulation')

%clear Received_m;

I_r = I_r';
Q_r = Q_r';

%% MATCHED FILTER BLOCK
t_T=(-5:Rb/fs:5)+1e-8;
r = 0.35; %roll off factor;
RRCF =(cos((1+r)*pi*t_T)+sin((1-r)*pi*t_T)./(4*r*t_T))./(1-(4*r*t_T).^2);
RRCF=RRCF/max(RRCF);

R_Ifilt = conv(I_r,RRCF);     % Signal after Matched Filter
R_Qfilt = conv(Q_r,RRCF);     % Signal after Matched Filter

clear I_r Q_r;

subplot(323)
plot(R_Ifilt)
title('I channel MF output')
subplot(324)
plot(R_Qfilt)
title('Q channel MF output')


%% SAMPLER BLOCK
%I Channel 
z_I = R_Ifilt;
K=length(Received_m);                % use 50 bit intervals to calculate average energy
mag=zeros(1,fs/Rb);  % find fs/Rb=440 average magnitudes for each bit
for i=1:fs/Rb        % interval has 440 samples
    mag(i)=sum(abs(z_I(i:fs/Rb:length(z_I)))); %sum up the magnitudes
    mag(i)=mag(i)/K;     % find the average magnitude
end
Pos=find(mag==max(mag));  % maximum position as optimal sampling time
Z_IT=z_I(Pos:fs/Rb:length(z_I)); %the output of sampler

%Q Channel
z_Q = R_Qfilt;
K=length(Received_m);                % use 50 bit intervals to calculate average energy
mag=zeros(1,fs/Rb);  % find fs/Rb=440 average magnitudes for each bit
for i=1:fs/Rb        % interval has 440 samples
    mag(i)=sum(abs(z_Q(i:fs/Rb:length(z_Q)))); %sum up the magnitudes
    mag(i)=mag(i)/K;     % find the average magnitude
end
Pos=find(mag==max(mag));  % maximum position as optimal sampling time
Z_QT=z_Q(Pos:fs/Rb:length(z_Q)); %the output of sampler

%clear R_Ifilt R_Qfilt;

subplot(325)
plot(Z_IT)
title('I channel Sampler output')
subplot(326)
plot(Z_QT)
title('Q channel Sampler output')

%% Length Checker (APPENDS ZEROS TO END OF MESSAGE IF TOO SMALL)
LenQ = length(Z_QT);
LenI = length(Z_IT);
if (LenQ > LenI)
    for i = 1:(LenQ-LenI)
        Z_IT(LenI + i) =0;
    end
end
if (LenI > LenQ)
    for i = 1:(LenI-LenQ)
        Z_QT(LenQ + i) =0;
    end
end

%% DIFFERENTIAL DECODER
y_k=Z_IT+j*Z_QT;  % combine the I Channel and Q Channel
z_k=y_k(1);
for t=2:length(y_k)
    z_k(t)=y_k(t)/y_k(t-1);
end
z_k; %output of the diff. decoder
%z_k=z_k*sqrt(2); %normal z_k by square root of 2
%separate the real term and imaginary terms
I_diffD=real(z_k);
Q_diffD=imag(z_k);
clear z_k Z_QT Z_IT;

%% Threshold DETECTOR BLOCK
%  converted back to string
I_bin = '1';
Q_bin = '1';
for i = 1:length(I_diffD)
    if(I_diffD(i) > 0)
        I_bin(i) = '1';
    else
        I_bin(i) = '0';
    end
end
for k = 1:length(Q_diffD)
    if(Q_diffD(k) > 0)
       Q_bin(k) = '1';
    else
        Q_bin(k) = '0';

    end
end


%% Channel Recombination Block
recoverd_message = '0'; 

for i= 1:length(I_diffD)
    recovered_message(2*i -1) = I_bin(i);
end
for i = 1:length(Q_diffD)
    recovered_message(2*i) = Q_bin(i);
end
%% remove startbits
clear I_bin Q_bin;

 %% Numberizer
for i = 1:length(recovered_message)
    message_num(i) = str2num(recovered_message(i));
end


barkerSeq = [ 1 1 1 1 1 0 0 1 1 0 1 0 1];
endFrameFormat = [ '1' '1' '1' '1' '1' '1' '1' '1']; %ending sequence 111



 weight = 0;
 for i = 1:(length(recovered_message)-length(barkerSeq))-2      
     weight(i) = sum(xor(barkerSeq,message_num(i:i+length(barkerSeq)-1)));
 end
 [mag startPos] = min(weight);


Z = message_num(startPos+length(barkerSeq): startPos + length(barkerSeq)+ 2*8*msgSize - 1);


   

%%  VITERBI DECODER
% Initialize all the parameters
Sa.state = [0 0];
Sa.state_metric = 0;
Sa.prev_state = [0 0];
Sa.path_metric = [0 0];

Sb.state = [1 0];
Sb.state_metric = 20;
Sb.prev_state = [1 0];
Sb.path_metric = [0 0];

Sc.state = [0 1];
Sc.state_metric = 20;
Sc.prev_state = [0 1];
Sc.path_metric = [0 0];

Sd.state = [1 1];
Sd.state_metric = 20;
Sd.prev_state = [1 1];
Sd.path_metric = [0 0];

% define a metric matrix for the entire process
metric = zeros(4,length(Z)/2);
metric(1,1) = 0;
metric(2,1) = 20;
metric(3,1) = 20;
metric(4,1) = 20;

% define a prev_state matrix for the entire process
prev_state = zeros(4,length(Z)/2);
prev_state(1,1) = 1;
prev_state(2,1) = 2;
prev_state(3,1) = 3;
prev_state(4,1) = 4;

for t=1:length(Z)/2
    
    % State a
    Sa(t+1).path_metric(1) = pdist([Z(2*t-1:2*t); 0 0],'hamming')*2+Sa(t).state_metric;
    Sa(t+1).path_metric(2) = pdist([Z(2*t-1:2*t); 1 1],'hamming')*2+Sc(t).state_metric;
    j = find(Sa(t+1).path_metric==min(Sa(t+1).path_metric));
    j = j(1); % if two path_metrics are the same, pick the first one
    if j==1
        Sa(t+1).state_metric = Sa(t+1).path_metric(1);
        metric(1,t+1) = Sa(t+1).state_metric;
        Sa(t+1).prev_state = [0 0];
        prev_state(1,t+1) = 1;
    else
        Sa(t+1).state_metric = Sa(t+1).path_metric(2);
        metric(1,t+1) = Sa(t+1).state_metric;
        Sa(t+1).prev_state = [0 1];
        prev_state(1,t+1) = 3;
    end
    
    
    % Sate b
    Sb(t+1).path_metric(1) = pdist([Z(2*t-1:2*t); 1 1],'hamming')*2+Sa(t).state_metric;
    Sb(t+1).path_metric(2) = pdist([Z(2*t-1:2*t); 0 0],'hamming')*2+Sc(t).state_metric;
    j = find(Sb(t+1).path_metric==min(Sb(t+1).path_metric));
    j = j(1); % if two path_metrics are the same, pick the first one
    if j==1
        Sb(t+1).state_metric = Sb(t+1).path_metric(1);
        metric(2,t+1) = Sb(t+1).state_metric;
        Sb(t+1).prev_state = [0 0];
        prev_state(2,t+1) = 1;
    else
        Sb(t+1).state_metric = Sb(t+1).path_metric(2);
        metric(2,t+1) = Sb(t+1).state_metric;
        Sb(t+1).prev_state = [0 1];
        prev_state(2,t+1) = 3;
    end

    
    % State c
    Sc(t+1).path_metric(1) = pdist([Z(2*t-1:2*t); 1 0],'hamming')*2+Sb(t).state_metric;
    Sc(t+1).path_metric(2) = pdist([Z(2*t-1:2*t); 0 1],'hamming')*2+Sd(t).state_metric;
    j = find(Sc(t+1).path_metric==min(Sc(t+1).path_metric));
    j = j(1); % if two path_metrics are the same, pick the first one
    if j==1
        Sc(t+1).state_metric = Sc(t+1).path_metric(1);
        metric(3,t+1) = Sc(t+1).state_metric;
        Sc(t+1).prev_state = [1 0];
        prev_state(3,t+1) = 2;
    else
        Sc(t+1).state_metric = Sc(t+1).path_metric(2);
        metric(3,t+1) = Sc(t+1).state_metric;
        Sc(t+1).prev_state = [1 1];
        prev_state(3,t+1) = 4;
    end

    
    % State d
    Sd(t+1).path_metric(1) = pdist([Z(2*t-1:2*t); 0 1],'hamming')*2+Sb(t).state_metric;
    Sd(t+1).path_metric(2) = pdist([Z(2*t-1:2*t); 1 0],'hamming')*2+Sd(t).state_metric;
    j = find(Sd(t+1).path_metric==min(Sd(t+1).path_metric));
    j = j(1); % if two path_metrics are the same, pick the first one
    if j==1
        Sd(t+1).state_metric = Sd(t+1).path_metric(1);
        metric(4,t+1) = Sd(t+1).state_metric;
        Sd(t+1).prev_state = [1 0];
        prev_state(4,t+1) = 2;
    else
        Sd(t+1).state_metric = Sd(t+1).path_metric(2);
        metric(4,t+1) = Sd(t+1).state_metric;
        Sd(t+1).prev_state = [1 1];
        prev_state(4,t+1) = 4;
    end

end

%% Trace back (Viterbi con't)
N=length(Z)/2+1;
j = find(metric(:,N)==min(metric(:,N))); % j represents the state
j = j(1);
decoded_mesg = zeros(1,length(Z)/2);
path = zeros(1,length(Z)/2);

for i=N:-1:2
    path(i-1) = prev_state(j,i);
    
    if j==1
        decoded_mesg(i-1) = 0;
    elseif j==2
        decoded_mesg(i-1) = 1;
    elseif j==3
        decoded_mesg(i-1) = 0;
    elseif j==4
        decoded_mesg(i-1) = 1;
    end
    
    j = prev_state(j,i);
end

message_str = '0';
for i = 1:length(decoded_mesg)
    if (decoded_mesg(i) == 1)
        message_str(i) = '1';
    else 
        message_str(i) = '0';
    end
end

%% Convert message back into ASCII
count8 = 1;
countAscii = 1;
holdCharacter = '0';
recov_inMesg = 0;
for i = 1:length(message_str)
    holdCharacter(count8) = message_str(i);
    count8 = count8 +1; 
    if (count8 == 9)
    count8 = 1;
    recov_inMesg(countAscii) = bin2dec(holdCharacter);
    countAscii = countAscii +1; 
    end
end

recovered_mesg = char(recov_inMesg)

