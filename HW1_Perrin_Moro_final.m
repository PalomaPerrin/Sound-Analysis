function  exportPGFPlots( filename, varargin )
clc; clear; close all;

%% Your name(s), student ID number(s)
%-------------------------------------------------------------------------%
% Paloma Perrin, 10873558
% 
%-------------------------------------------------------------------------%

%% Here are some MATLAB function that you might find useful:
% audioread, soundsc, flipud, fliplr, xcorr, eig, eigs, filter, toeplitz,
% fft, ifft, pause, disp, ...

%% Read the wav files 
% load the source signal x
[s_x, Fs_x] = audioread('x.wav');

% load the microphone signal y
[s_y, Fs_y] = audioread('y.wav');


%% Parameters
% number of filter taps

M = 4000; 

% length of the source signal x

N = length(s_x);

%% Wiener-Hopf solution
% compute the autocorrelation matrix
r = xcorr(s_x);
r_p = r((length(r)/2):(length(r)/2)+(M-1))/length(s_x);%we take all the positive lags to have a positive diagonal in the matrix
R = toeplitz(r_p);%we take toeplitz matrix of the autocorrelation
%the number of thap is the number of equations that we want to solve


% compute the cross-correlation vector
p = xcorr(s_x,s_y);
p = p(length(s_y):-1:length(s_y)-M+1)/length(s_y); %we take all the negative lags to respect the formula


% compute the Wiener-Hopf solution
w_o = R^(-1)*p;


%% Steepest Descent
% determine the step-size parameter 

% The necessary & sufficient condition for convergence or stability is
% govern by the largest eigenvalue : 
% 0<mu<[2/(largest eigenvalue)]

e=eig(R);
e_max= max(e);
e_lim=2/e_max; %=13.0456
mu = 13;  

% determine the global time constant of the Steepest Descent algorithm

%The convergence time is governed by the smallest eigenvalue : 
%tau=1/(2*mu*smallest eigenvalue)

e_min=min(e);
tau=1/(2*mu*e_min);

% initialize the vector of filter taps
w = zeros(M,1);

% perform the iterative update

% we have to compute the iteration for at least 2000 gradients updates
for ii=1:2000
    w=w+mu*(p-R*w);
end

% compute the theoretical minimum MSE, we applied the formula
J_min = var(s_y)-p'*R^-1*p;

% compute the MSE associated with the current taps estimate, we applied the
% formula
J_w = J_min+(w-w_o)'*R*(w-w_o);

% compute and discuss the cost functions ratio 
% the ratio is around 0,99; we expected 1 to have the best filter response
% it validates our reasoning
ratio = J_w / J_min;


%% Apply time-domain filter

% we compute the convolution of the input signal x with the Steepest
% Descent filter w, the filter taps vector

y_hat_td = filter(w,1,s_x);


%% playback: source signal x
disp('Press any key to start playback the input signal x') 
pause()
soundsc(s_x,Fs_x);


%% playback: microphone signal y
disp('Press any key to start playback the output signal y') 
pause()
soundsc(s_y,Fs_y);

%% playback: time-domain filtered signal y_hat
disp('Press any key to start playback the time-domain filtered signal y_hat') 
pause()
soundsc(y_hat_td,16000)  

%% Filtering in the frequency domain
% determine the length of the signal after zero-padding

%The anti-aliasing condition says that the frequency of the signal must be
%at least twice the sampling frequency
%lenght(fft(s_x))=14437, 2^15=32768, fs=16000

L_w=length(w); %length of the Wiener taps vector
N_y=N+L_w-1; %minimal length of the fft to respect anti-aliasing condition
s_zeropad=cat(1,s_x,zeros(18331,1)); %zero-padding the input signal x to the length 2^15 
wpad=cat(1,w,zeros(28768,1)); %zero-padding the filter taps vector x to the length 2^15
L = length(s_zeropad); 

% compute the spectra
W = fft(wpad);
X = fft(s_zeropad);

% perform frequency-domain filtering
Y = X.*W; %convolution in the time domain = multiplication in the frequency domain

% transform back to time domain
y_hat_fd = ifft(Y);


%% playback: freqeuncy-domain filtered signal
disp('Press any key to start playback the frequency-domain filtered signal y_hat_fd') 
pause()
soundsc(y_hat_fd,16000);

%% OLA Filtering
% window length
wlen = 256; %=M

% hop-size (must respect the COLA condition)
hop = wlen/2; %hop size of a triangular window

% define a tapered window
win = triang(wlen);

% determine the length of the windowed signal after zero-padding
%To avoid time-domain aliasing, we need to zero-pad each block up to a
%length superior to L_ola = wlen+L_w-1
%L_ola=4255 and the closest power of two is 2^13=8192
L_ola = 2^(ceil(log2(wlen+L_w-1))); 

% compute the fft of the Wiener filter
W = fft(w,L_ola);

% compute the total number of frames to be processed via OLA
n_frames = floor((N-wlen)/hop)+1;

% initialize the output filtered signal
y_hat_ola = zeros(N+L_ola,1);

% implement the OLA algorithm 
for n=0:n_frames-1
    index = n*hop+1:n*hop + wlen; %indices for the nth frame of the input signal
    x_bck = s_x(index).*win; %triangle windowed nth frame of input signal x
    X_ft = fft(x_bck,L_ola); % computes FFT of the input signal including zero padding to L_ola
    o_index = n*hop+1:n*hop + L_ola; %indices for the nth frame of the output signal
    y_bck = ifft(X_ft.*W); %computes ifft of each block filtered by the wiener filter
    y_hat_ola(o_index) = y_hat_ola(o_index) + y_bck;  %reconstruction of the signal in time domain
    %by adding every filtered samples
end
 

%% playback: OLA filtered signal y_hat_ola
disp('Press any key to start playback the OLA filtered signal y_hat_ola') 
pause()
soundsc(y_hat_ola,16000);

%% Plot the estimated filter taps against the true RIR
% load the RIR (g.wav)
[s_g, Fs_g] = audioread('g.wav');

%time defining
%it permits to plot the graph having the amplitude being function of time
dt=1/Fs_x; %time step
t=0:dt:(length(s_g)-1)/Fs_g; %time vector


% produce the plot
figure(1)
subplot(3,1,1)
plot (t,s_g,'r')
title('True RIR')
ylabel('Amplitude')
xlabel('Time [s]')
hold on
subplot(3,1,2)
plot(t,w_o,'b')
title('Wiener-Hopf filter taps')
ylabel('Amplitude')
xlabel('Time [s]')
hold on
subplot(3,1,3)
plot(t,w,'m')
title('Steepest Descent filter taps')
ylabel('Amplitude')
xlabel('Time [s]')
hold on
grid on

% save current figure as png
saveas(gcf, 'Figure_Perrin_Moro.jpg')

%Others plot that helped us understanding the homework 

%Plot of the time-domain filtered signal
figure(2)
plot(s_y)
hold on
plot(y_hat_td)
grid on
set(gca,'Xlim',[0,400])
legend('x','x filtered')
ylabel('Amplitude')
xlabel('Samples')
title('Input signal X  filtered in the Time Domain')

%Plot of the FFT filtered signal in the frequency domain
figure(3)
plot(abs(Y),'b')
grid on
ylabel('Amplitude')
xlabel('Samples')
title('Absolute value of the FFT of the imput signal x filtered in the frequency domain')

%Plot the frequency-domain filtered signal 
figure(4)
plot(y_hat_fd,'b')
hold on
plot(s_y,'r')
grid on
set(gca,'Xlim',[0,400])
ylabel('Amplitude')
xlabel('Samples')
title('Y hat fd : Input signal x filtered in the frequency domain')

%Plot the OLA filtered signal 
figure(5)
plot(y_hat_ola,'b')
hold on
plot(s_y,'r')
grid on
set(gca,'Xlim',[0,400])
ylabel('Amplitude')
xlabel('Samples')
title('Y hat ola : Input signal x filtered with OLA')


% ----EOF----