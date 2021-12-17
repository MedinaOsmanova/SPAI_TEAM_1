%Wout Peeters & Jeroen Coppens SPAI R&D Experience
%Measuring the frequency response of the microphone 
%% SECTION1
clear, clc, close all 

Fs = 44100;
t=0:(1/Fs):11;                    % 2 secs @ 8kHz sample rate
%Fs = 1/mean(diff(t));
%y=chirp(t,20,50,20000,'linear')';  % Start @ 200Hz, cross 100Hz at t=1sec
y = sweeptone(2,1,Fs,'SweepFrequencyRange',[20 20e3]);
figure(2)
plot(y)
%soundsc(y,Fs);
%Connecting to the audio devices
recObj = audiorecorder(Fs,16,1);
disp('Start speaking/playing sound.')
recordblocking(recObj, 12);
disp('End of Recording.');
%Speel af
%play(recObj);
%Store in array 
y_hat = getaudiodata(recObj);
figure(1);
plot(y);
%% SECTION2 
y_hat_down = y_hat(t<(12-(1/Fs)));
ir = impzest(y,y_hat_down);
figure(3)
plot((1:length(ir))/Fs, ir);

IR = fft(ir);
figure(4)
semilogx(linspace(0,Fs/2,length(IR)/2+1), db(abs(IR(1:end/2+1))))
xlim([30 Fs/2]), grid on
title('Impulse response (frequency)')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
%irEstimate = irEstimate(1:101);
% fft1 = fft(y_hat_down');
% fft2 = fft(y);
% h = (fft1/fft2);
% sys = tf(h);
% figure(3)
% step(sys)