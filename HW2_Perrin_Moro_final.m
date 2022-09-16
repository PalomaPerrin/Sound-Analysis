clc; clear; close all;

%% Your name(s), student ID number(s)
%-------------------------------------------------------------------------%
% Paloma PERRIN, 10873558
% 
%-------------------------------------------------------------------------%

% import the array data
load("array_data_64_mics.mat")

%% Parameters

% sampling frequency
Fs = 8000;
% speed of sound [m/s]
c = 340;
% number of sources
N_src = 2;

% number of microphones
M = 64;

% microphone signal length
N = size(y,2); %There are 64 microphones having 10000 components each

% determine the distance between two mics (see anti-aliasing condition)
%Anti aliasing condition : d < pi*c/(2*pi*Fs)=c/(fc*2)
%Microphone spacing d was selected to avoid spatial aliasing for every angle 
% of arrival and for every frequency up to the Nyquist frequency meaning
% that fc=Fs/2 to respect the sampling theorem

d=c/Fs;

%% Frequency estimation

% plot the spectrum of a microphone signal (use fftshift)

%We created a frequency axis which step is discretized with the microphone signal length
freq_axis = 0:Fs/N:Fs - Fs/N;
freq_axis = freq_axis - Fs/2; % We centered the data arounf the frequency 0 in order to have
%two symmetrical pics at f1= 500Hz and f2=700Hz

magnitude_spectrum = abs(fftshift(y(51,:))); %we chose randomly one microphone signal to plot  
%because their spectrum are all the same 

figure(1)
plot(freq_axis, magnitude_spectrum)
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Magnitude spectrum of a microphone signal')

%The fonction findpeaks helps us to know the frequency of the pics f1=500Hz
%& f2=700Hz
[PKS,LOCS]= findpeaks(magnitude_spectrum,Fs,'MinPeakHeight',500);

% estimate of the angular frequency of source signal 1
w_c1 = 2*pi*500; %w_c1=2*pi*f1

% estimate of the angular frequency of source signal 2
w_c2 = 2*pi*700; %w_c2=2*pi*f2

% row vector containing the estimated angular frequencies
w_c = [w_c1 w_c2];

%% Delay-and-sum beamformer

% angluar step between candidates for our DOA estimate
candidate_angle_step = 1*pi/180; %conversion of 1 degree step in radian

% the angles to be examined
candidate_angles = -pi/2:candidate_angle_step:pi/2; %angles from -90 to 90 degrees

% the number of the angles to be examined
%n_angles = length(candidate_angles);
n_angles = 181;

% pre-allocate the array for the DAS pseudo-spectrum
DAS_pseudo = zeros(n_angles,1);

% pre-allocate the array for the DOA estimates of the two sources 
DAS_DOA = zeros(N_src,1);

% pre-allocate the steering vector for the two sources
a = zeros(M, n_angles, N_src);

% pre-allocate the spatial filter coefficients for the two sources
h_DOA = zeros(M, N_src);

% sample estimate of the covariance matrix of the array data
R=(y*y')/N;
%Ms is the vector containing the steps of the propagation vector
Ms = 0:M-1;
Ms = Ms';

for ii=1:N_src % each i-th source...
    
    for j= 1:n_angles % every j-th angle to be examined...
        
% compute the spatial frequency
%we have a spatial frequency for the two sources for every angle defined
        w_s(j,ii) = (w_c(ii)*d*sin(candidate_angles(j)))/c;
       
% compute the steering vector for the given source and angle
%the propagation vector is then given by the formula :
         a(:,j,ii) = exp(-1i*Ms*w_s(j,ii));
         Mt= a(:,j,ii)'* a(:,j,ii); 
% save the value of the pseudo-spectrum for the given frequency and angle
        DAS_pseudo(j,ii) = (a(:,j,ii)'*R*a(:,j,ii))/Mt^2; 
        
    end    

    
% remove the spurious imaginary part from the DAS pseudo-spectrum
    DAS_pseudo = real(DAS_pseudo); 
    
% estimate the DOA as the angle for which the DAS pseudo-spectrum has
% the most prominent peak; DS=Delay and Sum 
%The most prominent peaks of the pseudo spectrum give us the angles towards
%which the microphones need to be steered to hear the better version of the
%source noise
     [maxDS, iDS] = max(DAS_pseudo);
     DAS_DOA(ii) = candidate_angles(iDS(ii));
    
% compute the spatial frequency associated to the DOA
%We now take only the spatial frequency of the two angles defining the
%sources
      w_s_DOA(ii) = w_c(ii)*d*sin(DAS_DOA(ii))/c;

% compute the steering vector associated to the DOA
%We now compute the propagation vector of the 64 microphones only for the
%two DoAs
      a_DOA = exp(-1i*Ms*w_s_DOA);

% compute and save the spatial filter associated to the angle of the DOA
% for the given source
      h_DOA = a_DOA/M;

% compute the spatial response using the pre-computed steering
% vectors and the spatial filter associated to the estimated DOA
    DAS_spatial_response = h_DOA(:,ii)'*a(:,:,ii);
    
        
% %plot the DAS pseudo-spectrum for the given source frequency
    figure(ii+1)
    subplot(211)
    plot(rad2deg(candidate_angles), DAS_pseudo(:,ii))
    xlabel('Angle [deg]')
    ylabel('Pseudo-spectrum')
    title("Delay-and-sum beamformer: pseudo-spectrum at " + num2str(w_c(ii)/(2*pi)) + " Hz")
     
% %plot the DAS spatial response (beam pattern) for the given source frequency
    figure(ii+1)
    subplot(212)
    polarplot(candidate_angles, real(DAS_spatial_response));
    title("Delay-and-sum beamformer: beam pattern at " + num2str(w_c(ii)/(2*pi)) + " Hz")
    thetalim([-90, 90])
    
end

%Display the DOA estimate (in degrees)
DAS_DOA_1_degrees = rad2deg(DAS_DOA(1));
DAS_DOA_2_degrees = rad2deg(DAS_DOA(2));

disp('DELAY-AND-SUM:')
disp(DAS_DOA_1_degrees)
disp(DAS_DOA_2_degrees)

%% Apply spatial filtering to the array data
%We computed the 2 dimension ifft of the signal and shifted it  
mic_time_domain_signal = real(ifft(fftshift(y),[],2));  
%we filtered the signal using the spatial filter defined for each DOA
% and took the real parts of their ifft
s_hat_1 =  real(ifft(h_DOA(:,1)'*y));
s_hat_2 = real(ifft(h_DOA(:,2)'*y)); 


%% Play a microphone signal and the two beamformer outputs

% % %play the mixture (mic signal)
disp('Press any key to start playback the mic_time_domain_signal') 
pause()
soundsc(mic_time_domain_signal(1,:), Fs)

% % %play the filtered singal (beamformer steered towards source 1)
disp('Press any key to start playback the y spatial filtered signal s_hat_1 of the first source') 
pause()
soundsc(s_hat_1, Fs)

% % %play the filtered singal (beamformer steered towards source 2)
disp('Press any key to start playback the y spatial filtered signal s_hat_2 of the seconf source') 
pause()
soundsc(s_hat_2, Fs) 


%% Parametric methods

%we have N_src(2)<M(64) and ws1 different from ws2
% eigenvalue decomposition of the covariance matrix R
[V,D] = eig(R);
e = eig(R);
%V is the matrix filled with the eigenvectors
%D is the matrix which diagonal is filled with eigenvalues

% sort the eigenvalues in descending order
%e_sort are the eigenvalues sorted in descending order
%Sort is the index associated to each eigenvalue sorted
[e_sort, Sort] =sort(e,'descend');  
e_mat = diag(e_sort); %diagonal matrix with sorted eigenvalues

% permute the columns of the eigenvector matrix accordingly
Q = V(:,Sort);


%% MUSIC

% retain the matrix whose columns span the noise subspace (matrix V)
V = Q(:,3:end); 

% pre-allocate the array for the MUSIC pseudo-spectrum
MUSIC_pseudo = zeros(n_angles,1);
% pre-allocate the array for the DOA estimates of the two sources
MUSIC_DOA = zeros(N_src,1);

for ii=1:N_src % each i-th source...
    
    for j= 1:n_angles % every j-th angle to be examined...
        
        % compute the MUSIC pseudo-spectrum using the precomputed steering vector
        MUSIC_pseudo(j,ii) =  1/((a(:,j,ii)')*(V*V')*a(:,j,ii));
        
    end
    
    % remove the spurious imaginary part
     MUSIC_pseudo = real(MUSIC_pseudo);
        
    % estimate the DOA as the angle for which the MUSIC pseudo-spectrum 
    % has the most prominent peak
    [maxmusic, iMUSIC] = max(MUSIC_pseudo);
     MUSIC_DOA(:,ii) = candidate_angles(iMUSIC);
        
    % plot the MUSIC pseudo-spectrum for the given source
    figure(4)
    subplot(2,1,ii)
    plot(rad2deg(candidate_angles), MUSIC_pseudo(:,ii))
    xlabel('Angle [deg]')
    ylabel('Pseudo-spectrum')
    title("MUSIC pseudo-spectrum at " + num2str(w_c(ii)/(2*pi)) + " Hz")
    
end

% display the MUSIC DOA estimate (in degrees)
MUSIC_DOA_1_degrees = rad2deg(MUSIC_DOA(1,1));
MUSIC_DOA_2_degrees = rad2deg(MUSIC_DOA(2,2));

disp('MUSIC:')
disp(MUSIC_DOA_1_degrees)
disp(MUSIC_DOA_2_degrees)

% EOF

%% QUESTION 3 

%Creation of the two sine signals
fc=20;
amp=1;
dt=1/Fs;
wc = fc*2*pi;
n=dt:dt:1;
s_1=amp*cos(2*pi*fc*n);
s_2=amp*cos(2*pi*fc*n+10);

figure('Name','Question 3.1')
subplot(3,1,1)
plot(n,s_1,'b');
title('Sine signal s_1(n)')
ylabel('amplitude')
xlabel('time[s]')
legend('s_1(n)')
subplot(3,1,2)
plot(n,s_2,'r');
legend('s_2(n)')
title('Sine signal s_2(n)')
ylabel('amplitude')
xlabel('time[s]')
subplot(3,1,3)
plot(n,s_1,'b',n,s_2,'r');
legend('s_1(n)','s_2(n)')
title('Sine signals s_1(n) and s_2(n)')
ylabel('amplitude')
xlabel('time[s]')

%fft of the two sine signals
Sw_1=fft(s_1);
Sw_2=fft(s_2);
S=[Sw_1 ; Sw_2];


%Choosing the DoAs of  the two sources
theta_1=deg2rad(50);
theta_2=deg2rad(-25);

%Creating the noise signal 
e3=rand(3,length(n))*0.6;
ef=fft(e3,[],1)

%we have a spatial frequency for the two sources with the angles that we
%know (chosen angles)
ws_1 = ((2*pi*fc)*d*sin(theta_1))/c;
ws_2 = ((2*pi*fc)*d*sin(theta_2))/c;
ws =[ws_1 ; ws_2]

% compute the steering vector for the given source and angle
%Given that we have 3 microphones, the propagation vector is :
M3=0:2;
M3=M3';
a_1 = exp(-1i*M3*ws(1)*sin(theta_1)/c);
a_2 = exp(-1i*M3*ws(2)*sin(theta_2)/c);
A3 = [a_1 , a_2];

y3 = A3*S + ef;

N3 =  size(y3,2);

freq_axis3 = 0:Fs/N3:Fs - Fs/N3;
freq_axis3 = freq_axis3 - Fs/2;
magnitude_spectrum3 = abs(fftshift(y3(1,:)));

figure('Name','Question 3.2')
subplot(211)
stem(freq_axis3, magnitude_spectrum3)
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Magnitude spectrum of a microphone signal')
subplot(212)
plot(real(ifft(fftshift(y3(1,:)))))
xlabel('samples')
ylabel('Magnitude')
title('Time domain signal y')

