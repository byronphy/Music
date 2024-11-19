clear; close all; clc

%% Input mp3
filename = 'JingleBells.mp3';
[y, fs] = audioread(filename); % y=mxn, for n-channel m timesteps; Fs=sampling rate [Hz]
y = y/max(y(:)); % normalization
[Nt,Nc] = size(y);
dt = 1/fs; % [s]
orig_duration = dt*length(y); % [s]
disp(['Original duration = ', num2str(orig_duration), ' s'])

% Spectrum
df = fs/Nt; %[Hz]
freq_all = (-Nt/2:Nt/2-1)*df;  % [Hz]
spec = fftshift(fft(y,[],1),1); 

% Plot
figure(1);
subplot(221)
plot((1:Nt)*dt, y(:,1),'r','LineWidth',1); 
xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 1')

subplot(222)
plot(freq_all, abs(spec(:,1)),'r','LineWidth',1); 
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 1')

subplot(223)
plot((1:Nt)*dt, y(:,2),'b','LineWidth',1); 
xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 2')

subplot(224)
plot(freq_all, abs(spec(:,2)),'b','LineWidth',1); 
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 2')

%% Low-Pass Filter
fc = 2000; % [Hz]
spec_lp = spec;
spec_lp(abs(freq_all) > fc,:) = 0;
y_lp = real(ifft(ifftshift(spec_lp,1),[],1));
y_lp = y_lp/max(y_lp(:)); % normalization

% Plot
figure(2);
subplot(221)
plot((1:Nt)*dt, y_lp(:,1),'r','LineWidth',1); 
xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 1')

subplot(222)
plot(freq_all, abs(spec_lp(:,1)),'r','LineWidth',1); 
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 1')

subplot(223)
plot((1:Nt)*dt, y_lp(:,2),'b','LineWidth',1); 
xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 2')

subplot(224)
plot(freq_all, abs(spec_lp(:,2)),'b','LineWidth',1); 
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 2')

%% Resample and weak out
fs_final = 100000; % [Hz]
y_lp_re = resample(y_lp, fs_final, fs);
N_duration = size(y_lp_re,1); % [t]
dt_final = 1/fs_final; % [t]

% Cut
duration = 40; % [s]
N_duration = round(duration/dt_final);
y_final = y_lp_re(1:N_duration,:);
t_final_all = (1:N_duration)*dt_final; % [t]

% Weak out
weak_duration = 3; % [s]
null_duration = 2; % [s]
N_weak = round(weak_duration/dt_final);
N_null = round(null_duration/dt_final);

% Scale
scale = ones(size(y_final));
scale(N_duration-N_weak-N_null+1:end,:)=repmat([linspace(1,0,N_weak),zeros(1,N_null)]',[1,2]);
y_final = y_final.*scale;
disp(['Total duration = ', num2str(N_duration*dt_final), ' s'])

% Spectrum
df_final = fs_final/N_duration; %[Hz]
freq_final_all = (-N_duration/2:N_duration/2-1)*df_final;  % [Hz]
spec_final = fftshift(fft(y_final,[],1),1); 
%spec_new = spec_new./max(spec_new); 

% Stop: clear sound
% sound(y_lp_re_wo,fs_new);  % LP resampled weak out
% sound(y_lp_re_wo(:,1),fs_new);  % channel 1
% sound(y_lp_re_wo(:,2),fs_new);  % channel 2

% Plot
figure(3);
subplot(221)
plot(t_final_all, y_final(:,1),'r','LineWidth',1); 
xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 1')

subplot(222)
plot(freq_final_all, abs(spec_final(:,1)),'r','LineWidth',1); 
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 1')

subplot(223)
plot(t_final_all, y_final(:,2),'b','LineWidth',1); 
xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 2')

subplot(224)
plot(freq_final_all, abs(spec_final(:,2)),'b','LineWidth',1); 
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 2')


%% Gradients
G_str_lim = 20; % [mT/m]
Gwave = y_final*G_str_lim; % [mT/m]
max_G_str = max(Gwave,[], 1);   % max strength, [mT/m]
disp(['Maximum Gx = ', num2str(max_G_str(1)), ' mT/m']);
disp(['Maximum Gy = ', num2str(max_G_str(2)), ' mT/m']);

Gslr = cat(1,diff(Gwave,1)/dt_final/1000,[0,0]); % slew rate, [mT/m/ms]
max_G_slr = max(Gslr,[], 1);   % max slew rate, [mT/m/ms]
disp(['Maximum dGx/dt = ', num2str(max_G_slr(1)), ' mT/m/ms']);
disp(['Maximum dGy/dt = ', num2str(max_G_slr(2)), ' mT/m/ms']);

disp(['Total duration = ', num2str(N_duration*dt_final), ' s']);

Gwave_xyz = cat(2, Gwave, zeros(N_duration,1));

% Plot
figure(4);
subplot(221)
plot(t_final_all, Gwave(:,1),'r','LineWidth',1); 
xlabel('Time (s)'); ylabel('Gx (mT/m)'); title('Gradient Strength')

subplot(222)
plot(t_final_all, Gslr(:,1),'r','LineWidth',1); 
xlabel('Time (s)'); ylabel('dGx/dt (mT/m/s)'); title('Gradient Slew Rate')

subplot(223)
plot(t_final_all, Gwave(:,2),'b','LineWidth',1); 
xlabel('Time (s)'); ylabel('Gy (mT/m)'); title('Gradient Strength')

subplot(224)
plot(t_final_all, Gslr(:,2),'b','LineWidth',1); 
xlabel('Time (s)'); ylabel('dGy/dt (mT/m/s)'); title('Gradient Slew Rate')


% %save('Gwave_xyz.dat','Gwave_xyz','-ascii')
% 
% output_dir = 'JingleBellsGwave_xyz.dat';
% % fid = fopen(output_dir, 'wb');
% % fwrite(fid, Gwave_xyz, 'double');
% % fclose(fid);
% [~, fn, ~] = fileparts(output_dir);
% dlmwrite([fn,'_',num2str(duration_sec),'s_',num2str(fc),'Hz.dat'],Gwave_xyz,'delimiter','\t');
% 
% % writematrix(Gwave_xyz, output_dir);

%% Pulseq
% Define system properties
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6) ;
% Create a new sequence object
seq = mr.Sequence(system) ;
% Set sequence parameters
Nrep = 1 ;% number of repetitions
Nx = 6400 ; % number of samples
adcDur = 256e-3 ; % ADC duration
rfDur = 500e-6 ; % RF duration
TR_fid = 500e-3; % unit: s
TR = 5000e-3 ; % unit: s
TE = 30e-3 ; % unit: s
% Create a non-selective pulse 
rf_ex = mr.makeBlockPulse(pi/2, 'Duration', rfDur, 'system', system) ;
% Define ADC event
adc = mr.makeAdc(Nx,'Duration',adcDur, 'system', system) ;
% Define delays
delayTE = TE - rf_ex.shape_dur/2 - rf_ex.ringdownTime - adc.delay ;
delayTR_fid = TR_fid - mr.calcDuration(rf_ex) - delayTE - mr.calcDuration(adc) ;
assert(delayTE >= 0) ;
assert(delayTR_fid >= 0) ;
% Define Gradients
g_duration = 40e-3; % [s]
N_g = g_duration/dt_final;
% gx = mr.makeArbitraryGrad('x', Gwave(1:N_g,1)*system.gamma/1e3);
% gy = mr.makeArbitraryGrad('y', Gwave(1:N_g,2)*system.gamma/1e3);
% % Define delay
% delayTR = TR - TR_fid - mr.calcDuration(gy);
N_g_rep = round(N_duration/N_g);
%N_g_rep = 300;
% Loop over repetitions and define sequence blocks
for i = 1:Nrep
    seq.addBlock(rf_ex) ;
    seq.addBlock(mr.makeDelay(delayTE)) ;
    seq.addBlock(adc) ;
    seq.addBlock(mr.makeDelay(delayTR_fid)) ;
    for j = 1:N_g_rep-1
        gx = mr.makeArbitraryGrad('x', Gwave(round(j*N_g):round((j+1)*N_g),1)*system.gamma/1e3, 'first',0,'last',0);
        gy = mr.makeArbitraryGrad('y', Gwave(round(j*N_g):round((j+1)*N_g),2)*system.gamma/1e3, 'first',0,'last',0);
        seq.addBlock(gx,gy);
    end
    %seq.addBlock(mr.makeDelay(delayTR));
end
% Check whether the timing of the sequence is compatible with the scanner
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
% Write to .seq
seq.write('fid_40000ms.seq')       % Write to pulseq file
% Plot sequence diagram
seq.plot() ;
% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points