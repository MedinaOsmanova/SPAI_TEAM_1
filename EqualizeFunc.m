function [y,Fs] = EqualizeFunc(x,plt)
%Help: 
%Input of this function is a filename: e.g.: x = 'file.wav'
%----------------
%On succes the function will return the equalized audio and the sampling
%frequency: you can listen to it using soundsc(y,Fs)
%----------------
%Plotting of the signals can be done by using the boolean plt
%plt = 1 : showing signals
%plt = 0 : hide the plots
[testSig,Fss] = audioread(x);
% Returning the sampling frequency
Fs = Fss;
% Make the signal mono
testSigMono =testSig(:,1)';
m = matfile('FinalRespons.mat');
hfinal_t = m.hfinal_t;
% Apply a convoltion with the response
result = conv(testSigMono,hfinal_t);
if plt == 1
    figure(1)
    plot(testSigMono)
    title("Non-equalized signal")
    figure(2)
    plot(result)
    title("Equalized sinal")
end
y = result; 
end

