close all;
clc; clear;
format longg;

Fs = 160e6; % system sample rate 160MHz
Fc = 2e6;   % central frequency (mixer frequency) 2MHz
R = 1e6;    % data rate 1Mbit/s (actually Msymbol/s)
OSR = 16;    % oversampling rate

phaseOffsetDegr=20;
phaseOffset=phaseOffsetDegr*pi/180;
fft_prec = 2^12;
samplesnum = 2^16; % simulation time in smples
bits_to_send = 2^12;
modulation_order = 2; % 0-BPSK, 1-QPSK, 2-QAM16, 3-QAM64, ..., QAM-2^(modulation_order*2) 
bits_per_symb = 2*modulation_order;
%symbols_to_send = ceil(bits_to_send / bits_per_symb);

% define simulation time vector
t=0:1/Fs:samplesnum/Fs;

sinw = sin(2*pi*Fc*t);
CW   = exp(2*pi*Fc*t*1j); % generating complex sinewave  by exponent 
CW   = cos(2*pi*Fc*t)+1j*sin(2*pi*Fc*t); % generating complex sinewave algebraic form

CW1  = exp(-2*pi*Fc*t*1j);

noise = .0001*rand(1,numel(CW));

%CW = CW - mean(CW);
%CW1 = CW1 - mean(CW1);

% first subplot: sin and cos
figure; 
subplot(1,2,1);
plot(t(1:1e3),[real(CW(1:1e3));imag(CW(1:1e3))]');
legend("sin(2\pi)","cos(2\pi)");
title("Complex sinwave in quadrature representation");

subplot(1,2,2); 
plot3(t(1:1e3),CW(1:1e3),'.-b');
hold on;
plot3(t(1:1e3),CW1(1:1e3),'.-r');
hold off;

% print('ScreenSizeFigure','-dpdf')

CW_10bit = floor(CW*2^10)/2^10;
CW1_14bit = floor(CW*2^14)/2^14;

figure; hold on;
pwelch(CW, abs(blackman(samplesnum)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); xlabel('Mhz');
pwelch(CW1 , abs(blackman(samplesnum)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); 
pwelch(CW_10bit, abs(blackman(samplesnum)),.2,2^10,Fs/1e6,'db','centerdc','none'); 
pwelch(CW1_14bit , abs(blackman(samplesnum)),.2,2^10,Fs/1e6,'db','centerdc','none'); 
legend("floatpoint CW","floatpoint -CW","fixedpoint 10bit CW","fixedpoint 14bit -CW");
hold off;



% subplot(3,1,2);pwelch(abs(fft(CW1+CW)));


% Generate uniformly distributed random binary data to transmitt
d_tx = randi([0 1], 1, bits_to_send);
%figure; plot (d_tx(1:100));
figure; stem(2*d_tx(1:100)-1);


% map bits to symbol
if modulation_order==0
	symb1=d_tx;
  symb=d_tx;
else
	s=1;
	d_tx = [d_tx,zeros(1,2*modulation_order-mod(numel(d_tx),2*modulation_order))]; % tail zeros
	
	d_tx1 = reshape(d_tx,bits_per_symb/2,ceil(numel(d_tx)/bits_per_symb*2));
	symb1 = bin2dec(num2str(d_tx1(:,1:end).'));
for k=1:bits_per_symb:bits_to_send
%		 binary value of symbol is 
		num2str(d_tx(k:k+modulation_order-1));
		symb(s) = bin2dec (num2str(d_tx(k:k+modulation_order-1)));
		s=s+1;
	endfor
endif

% grey coding
symb_gray = bitxor(symb1,floor(symb1/2)); % https://www.dsprelated.com/showthread/comp.dsp/96917-1.php

%tx_inphase   = zeros(1,ceil(numel(symb))/2);
%tx_quadrture = zeros(1,ceil(numel(symb))/2);

if modulation_order==0 % BPSK
    tx_inphase   = symb_gray*2-1;
	tx_quadrture = symb_gray*2-1;
else % QAM
    tx_inphase   = symb_gray(1:2:end)*2-2^(modulation_order)+1;
    tx_quadrture = symb_gray(2:2:end)*2-2^(modulation_order)+1;	% map even symbols t
endif

signal = tx_inphase + 1j*tx_quadrture;
figure; plot (signal ,'o');
axis('equal')

signalPhaseOffset = signal*(cos(phaseOffset)+1j*sin(phaseOffset));
figure; plot (signalPhaseOffset ,'o');
axis('equal')


%% interpolating filter design
ORDER = 62;
norm_match_fir = fir1 (ORDER, 1/OSR);
norm_match_fir = norm_match_fir/max(norm_match_fir);

upsampled_signal = upsample(signal,OSR);


% simulate FIR filter using input vector
complx_signal = filter(norm_match_fir,1,upsampled_signal);

%% simulate FIR filter using sample-by-sample loop
states=zeros(1,numel(norm_match_fir)-1);
for k_iter=1:numel(upsampled_signal)
  [complx_signal2(k_iter) , states] = filter (norm_match_fir, 1, upsampled_signal(k_iter), states);
end

zero_hold_tx_inphase = repelems(signal,[1:numel(signal); OSR*ones(1,numel(signal))]);


%figure; pwelch(complx_signal , abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); xlabel('Mhz');
figure; hold on; 
pwelch(resample(signal,OSR,1) , abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); 
pwelch(upsample(resample(signal,4,1),OSR/4) , abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); xlabel('Mhz'); 
xlabel('Mhz'); title("Compare resample vs upsample"); legend("resample", "upsample");
hold off;

figure; plot(complx_signal); hold on;
plot(tx_inphase + 1j*tx_quadrture,'o'); hold off;

figure; plot(imag(complx_signal)); hold on; 
plot([zeros(ORDER /2,1);real(upsampled_signal)],'-o');
plot([zeros(1,ORDER /2) real(zero_hold_tx_inphase)],'.-');
legend ("resampled","upsampled","zero order hold upsampled");
title ("Resampling in time domain");
hold off

figure; hold on;
pwelch(resample(signal,OSR,1), abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); xlabel('Mhz');
pwelch(upsample(signal,OSR) , abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none');
pwelch(zero_hold_tx_inphase , abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none');
legend ("resampled","upsampled","zero order hold upsampled");
title ("Resampling in frequency domain");
hold off

eyediagram(complx_signal(floor(ORDER/2):end),OSR,[],OSR-mod(ORDER,OSR/2)-1);


%% Quadrature modulation (frequency shift)
Fc = 50e6; 
x = resample(signal,OSR,1).';
samplenum = numel(x);
t=(0:samplenum-1)/Fs;
figure; hold on;
pwelch(x, abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); 
pwelch(exp(2j*pi*Fc*t).*x, abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); 
pwelch(real(exp(2j*pi*Fc*t).*x), abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); 
pwelch(imag(exp(2j*pi*Fc*t).*x), abs(blackman(fft_prec)),.2,fft_prec,Fs/1e6,'db','centerdc','none'); 
legend ("baseband signal","quadrature modulated signal", "real part of quadrature modulated signal","imag part of quadrature modulated signal");
hold off;
figure; plot( [real(x); real(exp(2j*pi*Fc*t).*x)].' );
legend ("baseband signal","quadrature modulated signal");


% animation of phase shift 
figure; 
qam=qammod(0:2^(2*modulation_order)-1,2^(2*modulation_order));
for phsh = 0:1:45;
  clf; hold on; title(sprintf ("animation of phase shift QAM%d on %d deg",2^(2*modulation_order),phsh));
  plot(qam,'ob');
  plot(exp(2j*pi/360*phsh)*qam,'or'); 
  pause(.01);
end  
hold off;