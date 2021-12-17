%Wout Peeters & Jeroen Coppens SPAI R&D Experience
%Measuring the frequency response of the microphone 
%% SECTION 1 - use the measured data to get x and y
clear, clc, close all;
[x,Fsx] = audioread('input short sine sweep.wav');    % is the input sine sweep
[y,Fsy] = audioread('sinereq_10s(cut)_ana.wav');    % measure the output through the mic

%% listen to the data input
%soundsc(x,Fsx);

%% listen to the data recorded
%soundsc(y,Fsy);
%% SECTION2
% plot the data
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
figure(2);
plot(ty,y);
%% SECTION3
%x_mono is a logarithmic sine sweep --> needs equalisation
x_mono =x(:,1);
y_mono =y(:,1);
x_rev= x_mono(end:-1:1);
% zero padd x to y 
x_revPad=[zeros(1,90000) x_rev']';
[Hx,win] = freqz(x_mono,1,2000);
Hxdb = 20*log10(abs(Hx));
figure(3)
semilogx(win/pi*Fsx/2, Hxdb);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
% plotting y_mono
[Hy,win] = freqz(y_mono,1,2000);
Hydb = 20*log10(abs(Hy));
figure(4)
semilogx(win/pi*Fsy/2, Hydb);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%% SECTION4
%proberen te compenseren voor logar sine
HRES = (Hy./Hx);
HRESdb = 20*log10(abs(HRES));
figure(5)
plot(win/pi*Fsy/2, HRESdb);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%% SECTION5 - smoothing met hamming window size = 1011
% smoothing after dB scaling!
L = 101;
ham = hamming(L)/((L-1)/2);
Hdb = 20*log10(abs(HRES));
H_smooth = conv(Hdb,ham,'same');
figure(6)
freqs=win/pi*Fsy/2;
plot(freqs, H_smooth);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('smoothing with hamming window');

%% SECTION 6 - knippen van data van frequentierespons
% 20000 samples met frequentieresolutie van win*Fsy/(2*pi) = win*7.0187e+03
% Eindpunt voor X = 1.9415e+04, heeft index = 17611 (export cursorinfo van grafiek)
% Beginpunt voor X = 25.3575 (dichtbij 20 Hz + marge) 
% => check freqs variabele, komt overeen met index 24
freqs_cut=freqs(24:1761);
H_cut=H_smooth(24:1761);
figure(7)
plot(freqs_cut, H_cut);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('plot after cutting off edges');
Ncut=length(H_cut);
bool=(Ncut==length(freqs_cut));
%% SECTION 7 - testen van gekozen transferfunctie microfoon
% Testmethode: convolutie van gekozen transferfunctie met outputsignaal zou
% normaal gezien een dirac functie (irl een sinc functie) moeten geven.
gimmesinc=conv(H_cut,y_mono,'same');
figure(8)
plot(freqs_cut, gimmesinc);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%% SECTION8 Applying the response
hfinal_t = abs(ifft(H_cut)); 
%% save the response
save('FinalRespons.mat','hfinal_t')
%% 
[testSig,Fst] = audioread('sev5.wav');
testSigMono =testSig(:,1)';
result = conv(testSigMono,hfinal_t);
figure(9)
plot(testSigMono)
figure(10)
plot(result)
%% play input
soundsc(testSig,Fst);
figure(11)
%pspectrum(testSigMono,Fst,'spectrogram')
%% play result
soundsc(result,Fst);
figure(12)
%pspectrum(result,Fst,'spectrogram')