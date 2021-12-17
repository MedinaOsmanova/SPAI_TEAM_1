% Randall Ali = randall.ali@esat.kuleuven.be
% Created:           Oct 29, 2018
% Last update:       Nov 2, 2019
%
%
%
% For SPAI course 2019 - A framework for carrying out noise reduction in
% the STFT domain
%% Setup Configuration
% Here we set up the initial parameters for the simulations

addpath(genpath(fullfile('..')));   % add the previous folders 

[y_TD, fs] = audioread('noisy_speech_s0_n30_multi.wav');    % Read in the noisy signal
[x_TD, ~] = audioread('clean_speech_s0_single.wav');       % This is the clean speech only - use for comparison
M = size(y_TD,2);               % number of microphones
micsp = 0.05;                   % microphone spacing (known a priori)
source_angle_assumed = 0;       % Assumed source angle location (may or may not be used). 0 deg. corresponds to an endfire direction
cair = 340;                     % speed of sound in air (m/s)


%% Transform to STFT domain
% See functions - This implements WOLA method

%%%%%%% WOLA parameters - For STFT processing %%%%%%%

N_fft = 512;                        % number of FFT points
R_fft = N_fft/2;                    % shifting (50% overlap)
win = sqrt(hann(N_fft,'periodic')); % analysis window
N_half = floor(N_fft/2)+1;          % number of bins in onsided FFT 
freqs = 0:fs/N_fft:fs/2;            % one-sided set of actual frequencies

%%

% compute the STFT (frequency x frame x channel) and plot
y_STFT = calc_STFT(y_TD, fs, win, N_fft, N_fft/R_fft, 'onesided'); % noisy microphone signals
x_STFT = calc_STFT(x_TD, fs, win, N_fft, N_fft/R_fft, 'onesided'); % clean microphone signal
[N_freqs, N_frames] = size(y_STFT(:,:,1)); % total number of freq bins + time frames  
figure; imagesc(1:N_frames, freqs, mag2db(abs(y_STFT(:,:,1))), [-65, 10]); colorbar; axis xy; set(gcf,'color','w'), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Spectrogram of noisy speech'); 
figure; imagesc(1:N_frames, freqs, mag2db(abs(x_STFT(:,:,1))), [-65, 10]); colorbar; axis xy; set(gcf,'color','w'), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Spectrogram of clean speech'); 

% Compute the Speech Presence Probability on the first mic (reference)
[noisePowMat, SPP] = spp_calc(y_TD(:,1),N_fft,R_fft);
figure; imagesc(1:N_frames, freqs,SPP); colorbar; axis xy; set(gcf,'color','w');  set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Speech Presence Probability for ref mic');


%% Define a priori RTFs using the direct path model (if desired, comment out if not used).

a_pr = zeros(M,N_half);   % a-priori ATF
h_pr = zeros(M,N_half);   % a-priori RTF

for f = 1:(N_half)
    for m = 1:M
        % compute the corresponding time delay and ATF
        
        
    end
    h_pr(:,f) = a_pr(:,f)/a_pr(1,f);                % make into an RTF
end   
         
%% Define necessary constants

Rnn = cell(N_freqs,N_frames);  Rnn(:) = {zeros(M,M)};      % Noise Only (NO) corr. matrix. Initialize to zeros
Ryy = cell(N_freqs,N_frames);  Ryy(:) = {zeros(M,M)};      % Speech + Noise (SPN) corr. matrix. Initialize to zeros
lambda = 0.995;                                     % Forgetting factor for computing correlation matrices - change values to observe effects on results
SPP_thr = 0.8;                                      % Threshold for the SPP - also change values to observe effects on results

% Single Channel
sig_s = zeros(N_freqs,N_frames);
sig_n = zeros(N_freqs,N_frames);
G_sc_stft = zeros(N_freqs,N_frames);          % Single Channel gain
S_sc_stft = zeros(N_freqs,N_frames); 
Xi_min = 1e-6;  
alpha_n = 0.9;  alpha_s = 0.92;               % speech and noise forgetting factors in single channel reduction - change values to observe effects on results


% With A Priori RTF
S_mc_stft = zeros(N_freqs,N_frames);          % output speech estimate STFT domain 
W_mc = (1/M)*ones(M,N_freqs);                 % multi-channel filter


% STFT Processing framework
% compute filter, filter signal, save output in array
tic
for l=2:N_frames % Time index - start from 2 since things would involve "l-1" - this is just a convenience for the first iteration 
    
    for k = 1:N_freqs % Freq index
                
         
        %%%%%%%%%%%%%----- SINGLE CHANNEL PROCESSING ----------%%%%%%%%%%%%
        
        % Use SPP thresholding to compute the noise power
        % Compute speech power using the decision directed (or any other approach)
        % Compute the gain for each time-freq bin
        % Apply the gain to the noisy signal
        
        
        
        %%%%%%%%%%%%%----- MULTI-CHANNEL PROCESSING ----------%%%%%%%%%%%%
        
        % Use SPP thresholding to compute the Speech+Noise and Noise-Only Correlation matrices, i.e. Ryy and Rnn
        % Use either the a priori RTF vector or estimate the RTF vector
        % Compute the MVDR/MWF beamformer - you may also compute the MVDR beamformer first 
        % and then perform a single channel enhancement on the output from this.
                       
    end % end freqs
end % end time frames
toc



%% Transform back to time domain


% PLOT YOUR ENHANCED STFTs
times_stft = ((1:N_frames));
figure; imagesc(times_stft,freqs,mag2db(abs(x_STFT(:,:,1))), [-65, 15]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('speech component (target), 1st mic');
figure; imagesc(times_stft,freqs,mag2db(abs(y_STFT(:,:,1))), [-65, 10]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('microphne signal, 1st mic');
figure; imagesc(times_stft,freqs,mag2db(abs(S_sc_stft(:,:))), [-65, 10]); colorbar; axis xy; set(gcf,'color','w'); set(gca,'Fontsize',14);xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Single Channel Enhancement');
figure; imagesc(times_stft,freqs,mag2db(abs(S_mc_stft(:,:))), [-65, 10]); colorbar; axis xy; set(gcf,'color','w'); set(gca,'Fontsize',14);xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Multi Channel Enhancement');

% compute the iSTFT - convert back to time domain (plot to also have a look)
time = 0:1/fs:(length(y_TD(:,1))-1)/fs;
s_sc_TD = calc_ISTFT(S_sc_stft, win, N_fft, N_fft/R_fft, 'onesided');
s_mc_TD = calc_ISTFT(S_mc_stft, win, N_fft, N_fft/R_fft, 'onesided');


% Write audio files 
audiowrite('../audio_processed/enhanced_SC.wav',s_sc_TD(:,1),fs);
audiowrite('../audio_processed/enhanced_MC.wav',s_mc_TD(:,1),fs);
















