order    = 10; % order of the butterworth
fcutlow  = 1/140.; %in Hz
fcuthigh = 1/40.; %in Hz

trange_start = '2013-09-30/01:00:01';
trange_end = '2013-09-30/01:39:30';
trange_show = [datenum(trange_start), datenum(trange_end)];

%% %%%%%% load data
data_original = datain_simple('ths_fgs_fac_4matlab.dat', [4, 21600], 'double');
data_filtered = datain_simple('ths_fgs_fac_bpfilt_4matlab.dat', [4, 21600], 'double');

t = data_filtered(1,:); %% serial date number
ts = t*24*3600.; %% change t to seconds
s = data_filtered(end,:); %% Bz
fs = 1./(ts(2)-ts(1));

x_smoothfilter = data_filtered(end,:); %% filtered Bz for comparison

%% %%%%%% Use butterworth filter
% [b,a]    = butter(order,[fcutlow,fcuthigh]/(fs/2), 'bandpass'); % fs is the sampling frequency in Hz
% x        = filter(b,a,s); %x is the filtered time-series data, s is the input time-series data

%% %%%%%% Use Bessel
[x, b, a] = besselfilter(order,fcutlow,fcuthigh,fs,s);

%% Examine the results
i_plot = find((t > trange_show(1)) & (t < trange_show(2)));

figure(1);
subplot(2,1,1);
plot(t(i_plot), x_smoothfilter(i_plot));
xlim(trange_show);
datetick('x','HH:MM:SS')
subplot(2,1,2);
plot(t(i_plot), x(i_plot));
xlim(trange_show);
datetick('x','HH:MM:SS')