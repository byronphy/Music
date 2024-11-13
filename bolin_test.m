clear; close all; clc

filename = 'JingleBells.mp3';
[y, fs] = audioread(filename); % y=mxn, for n-channel m timesteps; Fs=sampling rate [Hz]
dt = 1/fs; % [s]
y = y(1:2e6,:);
y = y/max(y(:)); % normalization
[Nt,Nc] = size(y);
t_all = [1:Nt]*dt; % [s]

df = fs/Nt; %[Hz]
freq_all = (-Nt/2:Nt/2-1)*df;  % [Hz]
spec = fftshift(fft(y,[],1),1); 
spec = spec./max(spec); 

% Stop: clear sound
% sound(y,fs);  % original
% sound(y(:,1),fs);  % channel 1
% sound(y(:,2),fs);  % channel 2

% % Original waveforms
% figure(1);
% subplot(221)
% plot(t_all, y(:,1),'r','LineWidth',1); 
% xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 1')
% 
% subplot(222)
% plot(freq_all, abs(spec(:,1)),'r','LineWidth',1); 
% xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 1')
% 
% subplot(223)
% plot(t_all, y(:,2),'b','LineWidth',1); 
% xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 2')
% 
% subplot(224)
% plot(freq_all, abs(spec(:,2)),'b','LineWidth',1); 
% xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 2')


%% Low-Pass Filter
fc = 2000; % [Hz]
spec_lp = spec;
spec_lp(abs(freq_all) > fc,:) = 0;
y_lp = real(ifft(ifftshift(spec_lp,1),[],1));
y_lp = y_lp/max(y_lp(:)); % normalization

% Stop: clear sound
% sound(y_lp,fs);  % LP
% sound(y_lp(:,1),fs);  % channel 1
% sound(y_lp(:,2),fs);  % channel 2

% % LP waveforms
% figure(2);
% subplot(221)
% plot(t_all, y_lp(:,1),'r','LineWidth',1); 
% xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('LP Channel 1')
% 
% subplot(222)
% plot(freq_all, abs(spec_lp(:,1)),'r','LineWidth',1); 
% xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('LP Channel 1')
% 
% subplot(223)
% plot(t_all, y_lp(:,2),'b','LineWidth',1); 
% xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('LP Channel 2')
% 
% subplot(224)
% plot(freq_all, abs(spec_lp(:,2)),'b','LineWidth',1); 
% xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('LP Channel 2')


%% Resample and weak out
fs_new = 100000; % [Hz]
y_lp_re = resample(y_lp, fs_new, fs);
y_lp_re = y_lp_re(1:4e6,:);
Nt_new = size(y_lp_re,1); % [t]
dt_new = 1/fs_new; % [t]
t_new_all = [1:Nt_new]*dt_new; % [t]

% Weak out
N_weak = 3e5;
N_null = 2e5;
scale = ones(size(y_lp_re));
scale(Nt_new-N_weak-N_null+1:end,:)=repmat([linspace(1,0,N_weak),zeros(1,N_null)]',[1,2]);
y_lp_re_wo = y_lp_re.*scale;

df_new = fs_new/Nt_new; %[Hz]
freq_new_all = (-Nt_new/2:Nt_new/2-1)*df_new;  % [Hz]
spec_new = fftshift(fft(y_lp_re_wo,[],1),1); 
spec_new = spec_new./max(spec_new); 

% Stop: clear sound
% sound(y_lp_re_wo,fs_new);  % LP resampled weak out
% sound(y_lp_re_wo(:,1),fs_new);  % channel 1
% sound(y_lp_re_wo(:,2),fs_new);  % channel 2

% New waveforms
figure;
subplot(221)
plot(t_new_all, y_lp_re_wo(:,1),'r','LineWidth',1); 
xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 1')

subplot(222)
plot(freq_new_all, abs(spec_new(:,1)),'r','LineWidth',1); 
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 1')

subplot(223)
plot(t_new_all, y_lp_re_wo(:,2),'b','LineWidth',1); 
xlabel('Time (s)'); ylabel('Normalized Amplitude'); title('Channel 2')

subplot(224)
plot(freq_new_all, abs(spec_new(:,2)),'b','LineWidth',1); 
xlabel('Frequency (Hz)'); ylabel('Amplitude'); title('Channel 2')


%% Gradients
max_G_str = 20; % [mT/m]
Gwave = y_lp_re_wo*max_G_str; % [mT/m]
max_G_slr = max(diff(Gwave,1)/dt,[], 1);   % max slew rate, [mT/m/s]
duration_sec = Nt_new*dt_new; % [s]

Gwave_xyz = cat(2, Gwave, zeros(Nt_new,1));
%save('Gwave_xyz.dat','Gwave_xyz','-ascii')

output_dir = 'JingleBellsGwave_xyz.dat';
% fid = fopen(output_dir, 'wb');
% fwrite(fid, Gwave_xyz, 'double');
% fclose(fid);
[~, fn, ~] = fileparts(output_dir);
dlmwrite([fn,'_',num2str(duration_sec),'s_',num2str(fc),'Hz.dat'],Gwave_xyz,'delimiter','\t');

% writematrix(Gwave_xyz, output_dir);
