%% testmethode 2: berekende transferfunctie gebruiken om van y_mono terug te
% construeren naar perfecte sine sweep
% y=F(x)=conv(H,x) <=> x = F^-1(y) = deconv(y,H) 
% dus x^ = F_cut^-1(y) = deconv(y_mono,H_cut)
y_monods = downsample(y_mono,10);
H_cutds = downsample(H_cut,10);
x_est = deconv(y_monods,H_cutds);
%% 
Nx2=length(x_est);
tx2=0:1/Fsx:(Nx2-1)/Fsx;
figure(10)
plot(tx2, x_est);     
xlabel('time'); 
ylabel('Magnitude');
h_ifft=ifft(H_cutds);
%%
%proberen te compenseren voor logar sine
HRES = (fft(y_monods)/H_cut);
HRESdb = 20*log10(abs(HRES));
figure(5)
plot(win/pi*Fsy/2, HRESdb);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%%
x_est3 = deconv(y_monods,h_ifft);
%% 
Nx3=length(x_est3);
tx3=0:1/Fsx:(Nx3-1)/Fsx;
figure(11)
plot(tx3, x_est3);     
xlabel('time'); 
ylabel('Magnitude');

%%
x_est4 = deconv(fft(y_monods),H_cutds);
%% 
Nx4=length(x_est4);
tx4=0:1/Fsx:(Nx4-1)/Fsx;
figure(12)
plot(tx4, x_est4);     
xlabel('time'); 
ylabel('Magnitude');