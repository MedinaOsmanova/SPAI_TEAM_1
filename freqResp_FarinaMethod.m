%Wout Peeters & Jeroen Coppens SPAI R&D Experience
%Measuring the frequency response of the microphone 
%% SECTION 1 - use the measured data to get x and y
clear, clc, close all;
[x,Fsx] = audioread('input short sine sweep.wav');    % is the input sine sweep
[y,Fsy] = audioread('sinereq_10s(cut)_ana.wav');    % measure the output through the mic

%% SECTION2 - listen to the data
%soundsc(x,Fsx);

%% SECTION3 - listen to the data
%soundsc(y,Fsy);

%% SECTION4 - read in the data
close all;
Nx=length(x);
Ny=length(y);
tmax = (Nx-1)/Fsx;
tx=0:1/Fsx:(Nx-1)/Fsx;
ty=(0:1/Fsy:(Ny-1)/Fsy)';
Ny=length(y);
durationx = Nx/Fsx;         % takes 17.96 s
durationy = Ny/Fsy;% takes 20 s
figure(1);
plot(tx,x);
title("Clean input sweep")
xlabel("time(s)")
figure(2);
plot(ty,y);
title("Recorded sine sweep")
xlabel("time(s)")
%% SECTION5 - cut the data into right pieces
x_mono =x(:,1);
y_mono =y(:,1);
x_rev= x_mono(end:-1:1);
% zero padd x to y 
x_monoPad=[ x_mono' zeros(1,90000)]';
x_revPad=[zeros(1,90000) x_rev']';
% https://www.thinksrs.com/downloads/pdfs/applicationnotes/SR1_SweptSine.pdf
% (pagina4)
convTest = conv(x_rev,y_mono); 
%plotting a spectrogram if desired
%figure(3) 
%pspectrum(convTest,Fsx,'spectrogram')
%% SECTION6 - display the spectrogram of the clean and recorded sweep 
%figure(4)
%spectrogram(x_mono);
%title("Clean input sweep");
%figure(5)
%spectrogram(y_mono);
%title("Recorded input sweep");
%% SECTION7 - execute formula to get mic response h
% see: https://www.researchgate.net/publication/2456363_Simultaneous_Measurement_of_Impulse_Response_and_Distortion_With_a_Swept-Sine_Technique
fft1 = fft(y_mono);
fft2 = fft(x_monoPad);
H = fft1./fft2;
figure(6);
Hdb = 20*log10(abs(H));
N = length(H);
ftest = ((0:N-1)*Fsx/N)';
figure(5)
semilogx(ftest,Hdb);
title("Frequency response")
xlabel("Frequency(Hz)")
ylabel("Magnitude(dB)")
%showing the response between 25 and 19500 Hz
xlim([25 19500])
%% SECTION8 - smoothing the response
L = 1001;
ham = hamming(L)/((L-1)/2);
Hdb = 20*log10(abs(H));
H_smooth = conv(Hdb,ham,'same');
figure(6)
semilogx(ftest, H_smooth);
xlim([25 19500]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('smoothing with hamming window');
%% SECTION9 - extracting impulse response
impulse_t = ifft(H_smooth);
t_im=0:1/Fsx:(length(impulse_t)-1)/Fsx;
figure(7)
plot(t_im,impulse_t)
title("Impulse Response") 
ylabel("Magnitude")
xlabel("Time(s)")
xlim([0 0.035])
