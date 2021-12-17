%Wout Peeters & Jeroen Coppens SPAI R&D Experience
%Measuring the frequency response of the microphone 
%% SECTION 1 - use the measured data to get x and y
clear, clc, close all;
[x,Fsx] = audioread('input short sine sweep.wav');    % is the input sine sweep
[y,Fsy] = audioread('sinereq_10s(cut)_ana.wav');    % measure the output through the mic

%% listen to the data
soundsc(x,Fsx);

%%
soundsc(y,Fsy);

%% plot the data
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
%% cut the data into right pieces
% delay vinden tussen twee signalen ( gewenste en opgenomen )
x_mono =x(:,1);
y_mono =y(:,1);
x_rev= x_mono(end:-1:1);
% soundsc(x_rev,Fsx);
% zero padd x to y 
x_monoPad=[ x_mono' zeros(1,90000)]';
x_revPad=[zeros(1,90000) x_rev']';
% Bij het uitvoeren van onderstaande convolutie krijg je inderdaad een plot
% met verschillende pieken, het is met niet duidelijk hoe dat je dan verder
% moet mss hier:
% https://www.thinksrs.com/downloads/pdfs/applicationnotes/SR1_SweptSine.pdf
% (pagina4)
convTest = conv(x_rev,y_mono); 
%figure(3) 
%pspectrum(convTest,Fsx,'spectrogram')
%% test0
Delay = finddelay(x_mono,y_mono);
timeDelay = Delay/Fsx;
%% test1
figure(4)
spectrogram(x_mono);
title("Clean input sweep");
figure(5)
spectrogram(y_mono);
title("Recorded input sweep");
%% test2
x_mono = [x_mono(256:end)];
x_mono = [x_mono' zeros(1,90255)]';
%% execute formula to get mic response h
fft1 = fft(y_mono);
fft2 = fft(x_monoPad);
h = ifft(fft1./fft2);
figure(6);
H = 20*log10(abs(h));
N = length(h);
ftest = ((0:N-1)*Fsx/N)';
semilogx(ftest,H);
%% using freqz function of matlab to obtain repsonse
% Choosing 4400 points ( can be increased or decreased) 
[H,w] = freqz(y_mono,1,4400);
Hdb = 20*log10(abs(H));
figure(8)
semilogx(w/pi*Fsx/2, Hdb);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
% inverting the response 
Hdbin = -20*log10(abs(H));
figure(9)
semilogx(w/pi*Fsx/2, Hdbin);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
% compensating for response
Hdbfin = -20*log10(abs(H))+20*log10(abs(H));
figure(10)
semilogx(w/pi*Fsx/2, Hdbfin);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%% Testing it by performing a linear convolution (? correct) 
% inladen audio file, ik had eentje genomen van severity 5
[testSig,Fst] = audioread('sev5.wav');
testSigMono =testSig(:,1)';
h_t =abs(ifft(-abs(H)));
result = conv(h_t,testSigMono);
figure(11)
plot(testSigMono)
figure(12)
plot(result)
%% Deconv
H_norm= H/norm(H); 
[Q,R] = deconv(convTest,Hdb);
plot(abs(Q))
%% play input
soundsc(testSig,Fst);
figure(13)
%pspectrum(testSigMono,Fst,'spectrogram')
%% play result
soundsc(result,Fst);
figure(14)
%pspectrum(result,Fst,'spectrogram')


