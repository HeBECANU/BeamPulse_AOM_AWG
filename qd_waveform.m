%% User settings

% Use this to pick waveform type
wf_mode = 'burst';
% 'sweep' or 'burst'
%%
% General settings
sample_rate=1e9;
max_points=double(4e6);
points_min=double(32);
repeats_max=1e6;

% % Sweep settings
centre_freq = 2.1e6;%Hz
span_freq = -0.3e6;% Hz
amp_sweep = 0.15;% Vrms?
phase_sweep = 0;
% % PAL settings
% freq = 1.8e6; %Hz
amp_PAL = sqrt(2)*0.450; %Vrms
phase_PAL = 0;
cycles = 6;
dur = cycles/freq;
% tight PAL settings
% freq = 1.3e6; %Hz
% freq = e6; %Hz
% amp_PAL = sqrt(2)*0.650; %Vrms
% phase_PAL = 0;
% cycles = 6;
% dur = cycles/freq;
%%
% {'sine',freq(hz),phase(rad),amplitude,Gauss mod,sample rate,duration}
burst={
        {'sine',freq,phase_PAL,amp_PAL*sqrt(2),nan,sample_rate,dur},...
        {'const',0,sample_rate,1e-6}};
% 'sweep_sine', freq_start, freq_end, phase, amp, sample_rate, duration
sweep = {{'sweep_sine',centre_freq-0.5*span_freq,centre_freq+0.5*span_freq,phase_sweep,amp_sweep,sample_rate,1e-3}};
null_wf = {{'const',0,sample_rate,1e-6}};

if strcmp(wf_mode,'burst')
    wave_out = burst;
elseif strcmp(wf_mode,'sweep')
    wave_out = sweep;
end

cli_header('Transmitting waveform : %s',wf_mode);

chanels_dev1={ch_to_waveforms(wave_out),ch_to_waveforms(null_wf)};
send_segments(chanels_dev1,1)

