%BRYCE Wafeworm generatormag
%a simple and easy to use waveform programer for the 33600Amag
%able to produce long constant sections without using up valuable waveform space
%resonable documentation
%ability to repeat sub wavefroms
%built uig sections of code from  Roman Khakimovs 2015 OOP program

%to improve
%padd the waveform with a zero
%determine if it is possible to have sub waveforms with differing sample rates
%fixed

%% Sequence selection
% NOTE: choose one sequence from below
% NOTE: K is equal to the momentum recoil of the laser detuning
%--------------------------------------------------------------------------
% Segments:
%'mag_transfer':  Magnetic transfer
%'k=0,-1,-2':     delay + Bragg (|k=0> |--> |k=0> + |k=-1K> + |k=-2K>)
%'k=0,-1':        delay + Bragg (|k=0> |--> |k=0> + |k=-1K>)
%'k=-1,-2':       delay + Bragg (|k=0> |--> |k=-1K> + |k=-2K>)
%'mirror':        delay + Mirror pulse
%'splitter':      delay + 50:50 Beam Splitter
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Common Sequences:
% Mach-Zender:    {'mag_transfer','k=0,-1','splitter','mirror','splitter_MZ'}
% Rarity-Tapster: {'mag_transfer','k=0,-1,-2','mirror','splitter'}
%--------------------------------------------------------------------------
% sequence = {'k=0,-1,-2','splitter','mirror','splitter_MZ','mag_transfer'}; %Construct desired experimental sequence from segments above
% sequence = {'k=-1,-2','mirror','const','mag_transfer'};
%equence = {'const','mag_transfer',};%{,'const'};%{,'const','mag_transfer'};%{'k=0,-1'};
% sequence = {'k=0,-1'}%,'const','splitter','splitter_MZ','mag_transfer'};
 sequence = {'mag_k=-1'}%,'k=0,-1'}%mag_k=-1','const','k=0,-1','mirror,'splitter'}%,'splitter'}%,'mag_transfer'};%,'splitter','const','mag_transfer'}%,'mirror','splitter','const','mag_transfer'}%,'const'}%'const','mag_transfer'}%,'mag_transfer'};%,'mirror'};%,'splitter','mag_transfer', 
%% START user settings
doNormalization     =   false;
useVoltageScaling   =   true;
if_plot             =   true;
Trigg_Delay         =   0.0e-3;

%% END user settings
mat_fname='mat_wf/Raman_Bragg.mat';
global max_points srate_max points_min repeats_max
max_points=double(4e6);
srate_max=1e9;              % [S/s] max is 1e9
points_min=double(32);
repeats_max=1e6;

%% General configs
f0_AOM=80e6;%80.94e6;            % [Hz]     Central AOMs frequency
Ek=84.9e3;              % beam geometry (90 deg)
ampfun = @(b,x) b(1).*x(:,1).^b(2);
ampfun_sin = @(b,x) asin(sqrt(x(:,1)./b(1)))./b(2);
nullfun = @(b,x) b(1).*b(2).*x(:,1).*0;

%% PULSE SETTINGS

%%% MAGNETIC TRANSFER pulse
%--------------------------------------------------------------------------
% optimised for 1.35 ms
% B_trap_bottom=1.283e6;%298e6;%1.096e6*sqrt(2)-00e3;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;% %detuning (MHz)
% del = 1.8e3;%3e3;%-4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
% del2 = 0e3;%1.57e3;%0e3;
% T_pulse_del = 0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=92.5e-6; %95e-6;%80e-6;%80e-6;%25e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
% Gs_mod_R_mix=3.35;%2.7;%2.5;%3.2;%2.95;%0.95;%0.95;%0.015;%0.8;%1.8;%2.69400894;
% phi1_mix=pi;
% phi2=0;
% K_R_mix=0.1625;%0.23;%0.22;%0.22%;0.45;%0.21;%0.2;%0.4;%0.338;%0.338;%0.34063277;%0.31 %power

%optimised for 23.462 shunt=0.65 quad=2.0
B_trap_bottom=1.291e6;%0.99e6;%0.992e6;%298e6;%1.096e6*sqrt(2)-00e3;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;% %detuning (MHz)
del = 2.0e3;%2.5e3;%3e3;%-4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
del2 = 0e3;%1.57e3;%0e3;
T_pulse_del = 0.0e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
T_Raman_mix=160e-6; %95e-6;%80e-6;%80e-6;%25e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
Gs_mod_R_mix=0.5;%3.035;%5.5;%3.735;%2.7;%2.5;%3.2;%2.95;%0.95;%0.95;%0.015;%0.8;%1.8;%2.69400894;
phi1_mix=pi;
phi2=0;
K_R_mix=0.3485;%0.105;%0.09;%0.09;%0.23;%0.22;%0.22%;0.45;%0.21;%0.2;%0.4;%0.338;%0.338;%0.34063277;%0.31 %power


%optimised for 23.46135 shunt=0.65 quad=2.0, 1/2 wave plate at 240 degrees
% B_trap_bottom=0.49e6;%0.992e6;%298e6;%1.096e6*sqrt(2)-00e3;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;% %detuning (MHz)
% del = 0.0e3;%2.5e3;%3e3;%-4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
% del2 = 0e3;%1.57e3;%0e3;
% T_pulse_del = 0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=124e-6; %95e-6;%80e-6;%80e-6;%25e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
% Gs_mod_R_mix=0.6;%0.6;%3.035;%5.5;%3.735;%2.7;%2.5;%3.2;%2.95;%0.95;%0.95;%0.015;%0.8;%1.8;%2.69400894;
% phi1_mix=pi;
% phi2=0;
% K_R_mix=0.188;%0.188%0.0485;%0.105;%0.09;%0.09;%0.23;%0.22;%0.22%;0.45;%0.21;%0.2;%0.4;%0.338;%0.338;%0.34063277;%0.31 %power

%optimised for 23.46135 shunt=0.75 quad=1.0, 1/2 wave plate at 245 degrees
% p_bottom=0.48e6;%0.992e6;%298e6;%1.096e6*sqrt(2)-00e3;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;% %detuning (MHz)
% del = 0.0e3;%2.5e3;%3e3;%-4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
% del2 = 0e3;%1.57e3;%0e3;
% T_pulse_del = 0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=100e-6; %95e-6;%80e-6;%80e-6;%25e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
% Gs_mod_R_mix=0.2;%0.6;%3.035;%5.5;%3.735;%2.7;%2.5;%3.2;%2.95;%0.95;%0.95;%0.015;%0.8;%1.8;%2.69400894;
% phi1_mix=pi;
% phi2=0;
% K_R_mix=0.2;%0.188%0.0485;%0.105;%0.09;%0.09;%0.23;%0.22;%0.22%;0.45;%0.21;%0.2;%0.4;%0.338;%0.338;%0.34063277;%0.31 %power


%optimised for 23.46135 shunt=0.65 quad=2.0, 1/2 wave plate at 210 degrees
% B_trap_bottom=0.49e6;%0.992e6;%298e6;%1.096e6*sqrt(2)-00e3;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;% %detuning (MHz)
% del = 0.0e3;%2.5e3;%3e3;%-4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
% del2 = 0e3;%1.57e3;%0e3;
% T_pulse_del = 0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=124e-6; %95e-6;%80e-6;%80e-6;%25e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
% Gs_mod_R_mix=0.6;%0.6;%3.035;%5.5;%3.735;%2.7;%2.5;%3.2;%2.95;%0.95;%0.95;%0.015;%0.8;%1.8;%2.69400894;
% phi1_mix=pi;
% phi2=0;
% K_R_mix=0.2;%0.188%0.0485;%0.105;%0.09;%0.09;%0.23;%0.22;%0.22%;0.45;%0.21;%0.2;%0.4;%0.338;%0.338;%0.34063277;%0.31 %power


%optimised for 23.4614 shunt=0.65 quad=2.0
% B_trap_bottom=0.49e6;%0.32e6;%0.992e6;%298e6;%1.096e6*sqrt(2)-00e3;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;% %detuning (MHz)
% del = -1.0e3;%3.0e3;%2.5e3;%3e3;%-4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
% del2 = -0.0e3;%1.57e3;%0e3;
% T_pulse_del = 0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=100e-6; %95e-6;%80e-6;%80e-6;%25e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
% Gs_mod_R_mix=0.01;%0.5;%3.035;%5.5;%3.735;%2.7;%2.5;%3.2;%2.95;%0.95;%0.95;%0.015;%0.8;%1.8;%2.69400894;
% phi1_mix=pi;
% phi2=0;
% K_R_mix=0.35;%0.0485;%0.105;%0.09;%0.09;%0.23;%0.22;%0.22%;0.45;%0.21;%0.2;%0.4;%0.338;%0.338;%0.34063277;%0.31 %power

B_trap_bottom=3.375e6%419e3%3.375e6%480e3%419e3%3.37e6%419e3%394e3;%109e3%0.929e6;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;%
del = 84.96e3%-11e3;%-10e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
del2 = 0%34e3/%0e3;

T_pulse_del = 0.0e-6;%delay between pulses
dF_Raman=(B_trap_bottom) - del;     %[Hz]    Raman detuning
T_Raman_mix= 200e-6;%11.5e-6;%45e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
Gs_mod_R_mix= 10e-6;%0.6;%0.015;%0.8;%1.8;%2.69400894;
phi1_mix=pi;
phi2=0;
K_R_mix=0.62;%0.095;%%0.095;%0.4;%0.338;%0.338;%0.34063277;%0.31

f1_Raman_mix=f0_AOM+dF_Raman/2;%-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM-dF_Raman/2;%-del2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"


B_trap_bottom_intrap=4.2e6;%1.025e6;%1.0603782e6;
dF_Raman_in_trap=-(B_trap_bottom_intrap);     %[Hz]    Raman detuning
f1_Raman_mix_in_trap=f0_AOM-dF_Raman_in_trap/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix_in_trap=f0_AOM+dF_Raman_in_trap/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"


%%% DELAYs between pulses
%--------------------------------------------------------------------------
% T_delay_mix=3000e-6;      % Delay between the SRC and MIX pulse
T_delay_mix=0e-6;%700e-6;      % Delay between the MAG and Bragg pulse
T_delay_mix_1=0e-6
T_delay_mirror=220e-6;%210e-6%460e-6; %-40e-6;%835e-6;%750e-6;%75e-6;%500e-6; % Delay between the Bragg pulse and Mirror pulse
T_delay_mirror_2=5e-6; %75e-6;%500e-6; % Delay between the Bragg pulse and Mirror pulse
T_delay_splitter=150e-6;%460e-6;  %400e-6+40e-6;%440e-6-71.2e-6;%3000e-6;%500e-6; % Delay between the Bragg pulse and Mirror pulse
T_delay_splitter_MZ=210e-6;  %0e-6;%
T_delay=2.50e-3;%-470e-6;%1.55e-3;  %50e-6;%200e-6+40e-6;%+45e-6;%240e-6;

%%% phases
%--------------------------------------------------------------------------
% NOTE: separated out since we find zero gives fine results
phi1_Bragg=0;
phi2_Bragg=0;

phi1_mirror=0;
phi2_mirror=0;

phi1_splitter=0;
phi2_splitter=0;

phi1_splitter_MZ=0;%4.5;%-pi/4;
phi2_splitter_MZ=0;

%%
%%% MIRROR pulse
%--------------------------------------------------------------------------
dF_Bragg_1=245e3/2; %29e3/2%17.5e3;%18e3%0.064e6;%0.098e6-Ek;%0.103e6;%Ek;%0.112e6-Ek;%0.112e6/2;%0.112e6;%0.0858e6;%0.112e6+Ek/2;%0.115e6;%
dF_Bragg_2=245e3/2; %29e3/2%17.5e3;%18e3%0.064e6;%0.098e6-Ek;%0.064e6;%0.103e6;%Ek;%0.112e6-Ek;%0.112e6/2;%0.112e6;%0.0858e6;%0.112e6+Ek/2;%0.115e6;%0.0858e6;%
f1_Bragg_mirror=f0_AOM-dF_Bragg_1%1-3.85e6/2;
f2_Bragg_mirror=f0_AOM+dF_Bragg_2;

T_Bragg_mirror=200e-6%40e-6;%10E-6;%35E-6;%15E-6;%32E-6;
t0_Bragg_mirror=nan;%3.9895e-6;

% sinc_scale_Bragg_mirror=7.9e-6;%6.5e-6;%3.9e-6;%8e-6;%1.5e-6;%5.5e-6;%5.3e-6;%5.3e-6;
sinc_scale_Bragg_mirror_1=20e-6%7.68e-6;%7.286e-6;%12.000e-6;%8.129e-6;%7.843e-6;%3.9e-6;%6.1e-6;
sinc_scale_Bragg_mirror_2=20e-6%9.4e-6%7.68e-6;%5.714e-6;%6.65e-6;%7.15e-6;%3.9e-6;%6.1e-6;

% Amp_sinc_Bragg_mirror_1=0.329;%0.214;%sqrt(29.0);%sqrt(50.0);%sqrt(26.0);%sqrt(30.0);%
% Amp_sinc_Bragg_mirror_2=0.307;
Amp_sinc_Bragg_mirror_1=0.4%0.26;%0.18+0.02;%0.154;%0.144;%0.329;%sqrt(25.0);
Amp_sinc_Bragg_mirror_2=0.4%0.26;%0.18+0.02;%0.199;0.189;%0.307; %sqrt(25.0);5.714, 12.000


Chirp_grad_mirror = 0;%7.6e6;%

% wf_mirror_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);%.*cos(pi*(t-b(1)/2)/(b(1))).^2;
wf_mirror_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%sinc((t-b(1)/2)./b(4)).*b(3);%.*cos(pi*(t-b(1)/2)/(b(1))).^2;

% wf_mirror_pulse_1 = @(b,t) ampfun([60.117, 0.5638],abs(wf_mirror_pulse(b,t)).^2)/2e3.*sin(2*pi*b(2)*t).*sign(wf_mirror_pulse(b,t));%sinc((t-b(1)/2)./b(4)).*b(3).*sin(2*pi*b(2)*t);%
% wf_mirror_pulse_2 = @(b,t) ampfun([132.62, 0.5283],abs(wf_mirror_pulse(b,t)).^2)/2e3.*sin(2*pi*b(2)*t).*sign(wf_mirror_pulse(b,t));%sinc((t-b(1)/2)./b(4)).*b(3).*sin(2*pi*b(2)*t);%

%normal sinc pulse
% wf_mirror_pulse_1 = @(b,t) ampfun([136.5, 0.5234],abs(wf_mirror_pulse(b,t)).^2)/2e3.*sin(2*pi*b(2)*t+b(5)).*sign(wf_mirror_pulse(b,t));%sinc((t-b(1)/2)./b(4)).*b(3).*sin(2*pi*b(2)*t);%
% wf_mirror_pulse_2 = @(b,t) ampfun([62.91, 0.5534],abs(wf_mirror_pulse(b,t)).^2)/2e3.*sin(2*pi*b(2)*t+b(5)).*sign(wf_mirror_pulse(b,t));%sinc((t-b(1)/2)./b(4)).*b(3).*sin(2*pi*b(2)*t);%

%chirped sinc pulse
% wf_mirror_pulse_1 = @(b,t) ampfun([136.5, 0.5234],abs(wf_mirror_pulse(b,t)).^2)/2e3.*sin(2*pi*(b(2)-b(6).*t).*t+b(5)).*sign(wf_mirror_pulse(b,t));%sinc((t-b(1)/2)./b(4)).*b(3).*sin(2*pi*b(2)*t);%
% wf_mirror_pulse_2 = @(b,t) ampfun([62.91, 0.5534],abs(wf_mirror_pulse(b,t)).^2)/2e3.*sin(2*pi*(b(2)+b(6).*t).*t+b(5)).*sign(wf_mirror_pulse(b,t));%sinc((t-b(1)/2)./b(4)).*b(3).*sin(2*pi*b(2)*t);%

wf_mirror_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_mirror_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%%
%%% 50:50 Beam Splitter pulse (halo 400 mus delay)
%--------------------------------------------------------------------------
dF_Bragg_1=237e3/2;%18e3;%17
dF_Bragg_2=237e3/2;%18e3;%17
f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
f2_Bragg_splitter=f0_AOM+dF_Bragg_2;

T_Bragg_splitter=200e-6;%10E-6;%35E-6;
t0_Bragg_splitter=nan;

sinc_scale_Bragg_splitter_1=6.6e-6;%4.6e-6
sinc_scale_Bragg_splitter_2=6.6e-6;%4.6e-6

Amp_sinc_Bragg_splitter_1=0.24;%0.5
Amp_sinc_Bragg_splitter_2=0.24;%0.5

Chirp_grad_splitter = 0;%7.6e6;

% wf_splitter_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);
wf_mirror_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%

wf_splitter_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_splitter_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

% dF_Bragg_1=0.100e6;%0.114e6;%0.064e6*2;%0.112e6;%0.114429e6;%0.114e6;
% dF_Bragg_2=0.100e6;%0.114e6;%0.064e6*2;%0.112e6;%0.114429e6;%0.114e6;
f1_Bragg_splitter_MZ=f0_AOM-dF_Bragg_1;
f2_Bragg_splitter_MZ=f0_AOM+dF_Bragg_2;

Amp_sinc_Bragg_splitter_1_MZ=0.133+0.006;%sqrt(13.0+a_shift);
Amp_sinc_Bragg_splitter_2_MZ=0.133+0.006;%sqrt(13.0+a_shift);

% %%% 50:50 Beam Splitter pulse
% %--------------------------------------------------------------------------
% dF_Bragg_1=17e3;%0.064e6;%0.097e6;%0.112e6;%0.064e6*2;%
% dF_Bragg_2=17e3;%0.064e6;%0.097e6;%0.112e6;%0.064e6*2;%0.112e6;
% f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
% f2_Bragg_splitter=f0_AOM+dF_Bragg_2;
% 
% dF_Bragg_1=0.100e6;%0.114e6;%0.064e6*2;%0.112e6;%0.114429e6;%0.114e6;
% dF_Bragg_2=0.100e6;%0.114e6;%0.064e6*2;%0.112e6;%0.114429e6;%0.114e6;
% %regular aom frequency
% f1_Bragg_splitter_MZ=f0_AOM-dF_Bragg_1;
% f2_Bragg_splitter_MZ=f0_AOM+dF_Bragg_2;
% 
% % f0_AOM_shift = 80e6-1.5e6; 
% % f1_Bragg_splitter_MZ=f0_AOM_shift-dF_Bragg_1;
% % f2_Bragg_splitter_MZ=f0_AOM_shift+dF_Bragg_2;
% 
% T_Bragg_splitter=40E-6;%60E-6;%10E-6;%35E-6;
% t0_Bragg_splitter=nan;
% 
% sinc_scale_Bragg_splitter_1=4e-6;%4e-6;%15e-6;%7.286e-6;%8.129e-6;%6.1e-6;
% sinc_scale_Bragg_splitter_2=4e-6;%4e-6;%15e-6;%7.286e-6;%6.65e-6;%6.1e-6;
% 
% a_shift = +0.0;
% 
% Amp_sinc_Bragg_splitter_1=0.145;%0.18;%0.18+0.02;%0.143;%sqrt(13.0+a_shift);
% Amp_sinc_Bragg_splitter_2=0.145;%0.18;%0.18+0.02;%0.143;%sqrt(13.0+a_shift);
% 
% Amp_sinc_Bragg_splitter_1_MZ=0.133+0.006;%0.18-0.03;%-0.02;%sqrt(13.0+a_shift);
% Amp_sinc_Bragg_splitter_2_MZ=0.133+0.006;%0.18-0.03;%-0.02;%sqrt(13.0+a_shift);
% 
% Chirp_grad_splitter = 0;%7.6e6;
% 
% % wf_splitter_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);
% wf_mirror_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%
% 
% % wf_splitter_pulse_1 = @(b,t) ampfun([136.5, 0.5234],abs(wf_splitter_pulse(b,t)).^2)/2e3.*sin(2*pi.*(b(2)-b(6).*t).*t+b(5)).*sign(wf_splitter_pulse(b,t));
% % wf_splitter_pulse_2 = @(b,t) ampfun([62.91, 0.5534],abs(wf_splitter_pulse(b,t)).^2)/2e3.*sin(2*pi.*(b(2)+b(6).*t).*t+b(5)).*sign(wf_splitter_pulse(b,t));
% 
% wf_splitter_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
% wf_splitter_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
% 
% % wf_splitter_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*(sin(2*pi.*(20e3+b(2)-b(6).*t).*t+b(5))+sin(2*pi.*(-20e3+b(2)-b(6).*t).*t+b(5)));
% % wf_splitter_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*(sin(2*pi.*(20e3+b(2)-b(6).*t).*t+b(5))+sin(2*pi.*(-20e3+b(2)-b(6).*t).*t+b(5)));
%%
%%% Momentum splitting
%--------------------------------------------------------------------------

%%% Bragg splitting: |k=0> |--> |k=+1K> + |k=0> + |k=-1K>
dF_Bragg_1= 220e3/2;%0*42.5e3/2%42.48e3/2;%19e3;%106e3;%21.55e3;%64e3-Ek;%64e3.*1.5;%0.096e6;%0.096e6;%0.099e6;%~0.093;% %%
dF_Bragg_2= 220e3/2;%0*42.5e3/2%42.48e3/2;%19e3;%0.096e6;%106e3;%21e3;%0.096e6;%21.55e3;%64e3-Ek;%64e3.*1.5;%0.096e6;%0.096e6;%0.099e6; %%


% %% Raman pulses 
% random_manual_delta = 0.0;%0.01;
% dF_Bragg_sym = (0.97-0.085+random_manual_delta)*1e6/4;
% % dF_Bragg_sym = (2.1-0.085+random_manual_delta)*1e6/4;
% % dF_Bragg_sym = (3.062)*1e6/4;
% dF_Bragg_1 = dF_Bragg_sym;
% dF_Bragg_2 = dF_Bragg_sym;

% dF_Bragg_1= 42.5e3/2; %stupid_shit;%0.845e6/2; %1*1.8635e6/2; % Raman %42.5e3/2
% dF_Bragg_2= 42.5e3/2; %stupid_shit;%0.845e6/2; %1*1.8635e6/2; % Raman
% dF_Bragg_1 = 0;
% dF_Bragg_2 = 0;

f1_Bragg_sym_f=f0_AOM-dF_Bragg_1;%sym;
f2_Bragg_sym_f=f0_AOM+dF_Bragg_2;%sym;

T_Bragg_sym_f=60e-6;%37e-6;%%70e-6;%31.5E-6; %
% P_Bragg_f = 26;%7.2;%20.5;%~7.8 %%
K_Bragg_sym_f_1= 0.5;%0.6%0.45;%0.4;%%0.296+0.0115;%-0.02;%0.23;%ampfun([136.5, 0.5234],P_Bragg_f)/2e3;%0.08; 1.25 1.35
%[60.117, 0.5
K_Bragg_sym_f_2= 0.5; %0.6;%0.45;%0.4;%0.35;%%0.35;%0.85;%0.212+0.0115;%-0.02;%0.21;%1.1*ampfun([62.91, 0.5534],P_Bragg_f)/2e3;%0.08;%[60.117, 0.5638] multiplier ~1.2 1.4 1.33 1.45
%[132.62, 0.5283]
% Gs_mod_Bragg_src_f_1=6.0;%0.9*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.63806937736142e-01);%~0.7 0.95
% Gs_mod_Bragg_src_f_2=6.0;%0.87*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.28341254744861e-01);
Gs_mod_Bragg_sym_f_1=2.6e-6;%3.3e-6;%5.5e-6;%3.05e-6;%3.7e-6;%2.8e-6;%7.843e-6;%3.9e-6;%6.1e-6;
Gs_mod_Bragg_sym_f_2=2.6e-6;%3.3e-6;%5.5e-6;%3.3e-6;%3.05e-6;%3.7e-6;3.7e-6;%3.05e-6;%4.6e-6;%7.15e-6;%3.9e-6;%6.1e-6;
t0_Bragg_sym_f=nan;

wf_bragg_sym_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);%.*cos(pi*(t-b(1)/2)/(b(1))).^2;
wf_bragg_sym_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*(sin(2*pi.*(b(2)-b(6)).*t+b(5))) %+ sin(2*pi.*(b(2)+b(6)).*t+b(5)));
wf_bragg_sym_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*(sin(2*pi.*(b(2)-b(6)).*t+b(5)))%+ sin(2*pi.*(b(2)+b(6)).*t+b(5)));


% %%% Bragg splitting: |k=0> |--> |k=0> + |k=-1K> + |k=-2K> 
% dF_Bragg_1=0.09e6;%21e3;%106e3;%21.55e3;%64e3-Ek;%64e3.*1.5;%0.096e6;%0.096e6;%0.099e6;%~0.093;% %%
% dF_Bragg_2=0.09e6;%21e3;%0.096e6;%106e3;%21e3;%0.096e6;%21.55e3;%64e3-Ek;%64e3.*1.5;%0.096e6;%0.096e6;%0.099e6; %%
% f1_Bragg_src_f=f0_AOM-dF_Bragg_1;
% f2_Bragg_src_f=f0_AOM+dF_Bragg_2;
% 
% T_Bragg_src_f=50e-6;%70e-6;%31.5E-6; %%
% % P_Bragg_f = 26;%7.2;%20.5;%~7.8 %%
% K_Bragg_src_f_1=0.65;%0.296+0.0115;%-0.02;%0.23;%ampfun([136.5, 0.5234],P_Bragg_f)/2e3;%0.08; 1.25 1.35
% %[60.117, 0.5638]
% K_Bragg_src_f_2=0.65;%0.85;%0.212+0.0115;%-0.02;%0.21;%1.1*ampfun([62.91, 0.5534],P_Bragg_f)/2e3;%0.08;%[60.117, 0.5638] multiplier ~1.2 1.4 1.33 1.45
% %[132.62, 0.5283]
% % Gs_mod_Bragg_src_f_1=6.0;%0.9*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.63806937736142e-01);%~0.7 0.95
% % Gs_mod_Bragg_src_f_2=6.0;%0.87*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.28341254744861e-01);
% Gs_mod_Bragg_src_f_1=2.8e-6;%2.8e-6;%7.843e-6;%3.9e-6;%6.1e-6;
% Gs_mod_Bragg_src_f_2=2.8e-6;%4.6e-6;%7.15e-6;%3.9e-6;%6.1e-6;
% t0_Bragg_src_f=nan;
% 
% wf_bragg_src_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);%.*cos(pi*(t-b(1)/2)/(b(1))).^2;
% wf_bragg_src_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
% wf_bragg_src_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));


%%
%%% Bragg splitting: |k=0> |--> |k=0> + |k=-1K>

% dF_Bragg_1=65e3;%0.06e6;
% dF_Bragg_2=65e3;%0.06e6;
% f1_Bragg_src_t=f0_AOM-dF_Bragg_1;
% f2_Bragg_src_t=f0_AOM+dF_Bragg_2;
% 
% T_Bragg_src_t=34e-6; %32
% P_Bragg_src = 5.8;%25;%4.8; %4.2; %power in mW 7 to 9 works well 4.93
% K_Bragg_src_1=ampfun([60.117, 0.5638],P_Bragg_src)/2e3;
% K_Bragg_src_2=ampfun([132.62, 0.5283],P_Bragg_src)/2e3;
% Gs_mod_Bragg_src_1=1.3*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01); %1.82
% Gs_mod_Bragg_src_2=1.3*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01); %1.81
% 
% t0_Bragg_src_t=nan;
dopp_shft = 11e3%232e3;%23e3 %146.6e3
dF_Bragg_1= (84.96e3+dopp_shft)/2; %64e3-Ek;%64e3;%0.06e6;
dF_Bragg_2= (84.96e3+dopp_shft)/2; %64e3-Ek;%64e3;%0.06e66;

f1_Bragg_src_t=f0_AOM-dF_Bragg_1;
f2_Bragg_src_t=f0_AOM+dF_Bragg_2;

T_Bragg_src_t=200e-6;%45e-6;%30e-6;%45e-6;%45 %32
% P_Bragg_src = 5.8;%25;%4.8; %4.2; %power in mW 7 to 9 works well 4.93
K_Bragg_src_1=0.62%0.28;%0.14%0.28;%0.43;%0.75;%0.6;%0.15;%0.0975;%0.4;%0.23;%0.2;%0.5;%0.8;%1;%0.8;%0.62;%0.2;%0.13+0.035;%ampfun([60.117, 0.5638],P_Bragg_src)/2e3;
K_Bragg_src_2=0.62%0.28;%0.14;%0.28;%0.43;%0.75;%0.6;%0.15;%0.0975;%0.4;%0.23;%0.2;%0.5;%0.08;%1;%0.8;%0.62;%0.8;%0.1+0.035;%ampfun([132.62, 0.5283],P_Bragg_src)/2e3; %change
Gs_mod_Bragg_src_1=3.8e-6%5.6e-6%6.4e-6;%6.5e-6;%9e-6;%8.129e-6;%1.3*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01); %1.82
Gs_mod_Bragg_src_2=3.8e-6%5.6e-6%6.4e-6;%6.5e-6;%8e-6;%9e-6;%6.65e-6;%1.3*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01); %1.81

t0_Bragg_src_t=nan;
wf_Bragg_t_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%
wf_Bragg_src_t_pulse_1 = @(b,t) wf_Bragg_t_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_Bragg_src_t_pulse_2 = @(b,t) wf_Bragg_t_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%%
%%% Bragg splitting: |k=0> |--> |k=0> + |k=+1K>

% dF_Bragg_1=65e3;%0.06e6;
% dF_Bragg_2=65e3;%0.06e6;
% f1_Bragg_src_t=f0_AOM-dF_Bragg_1;
% f2_Bragg_src_t=f0_AOM+dF_Bragg_2;
% 
% T_Bragg_src_t=34e-6; %32
% P_Bragg_src = 5.8;%25;%4.8; %4.2; %power in mW 7 to 9 works well 4.93
% K_Bragg_src_1=ampfun([60.117, 0.5638],P_Bragg_src)/2e3;
% K_Bragg_src_2=ampfun([132.62, 0.5283],P_Bragg_src)/2e3;
% Gs_mod_Bragg_src_1=1.3*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01); %1.82
% Gs_mod_Bragg_src_2=1.3*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01); %1.81
% 
% t0_Bragg_src_t=nan;
dopp_shft = 220e3;%23e3;%34e3
dF_Bragg_1= (84.96e3-dopp_shft)/2; %64e3-Ek;%64e3;%0.06e6;84.96e3-
dF_Bragg_2= (84.96e3-dopp_shft)/2;%64e3-Ek;%64e3;%0.06e66;
f1_Bragg_src_up=f0_AOM+dF_Bragg_1;
f2_Bragg_src_up=f0_AOM-dF_Bragg_2; 
         
T_Bragg_src_up= 60e-6;%45e-6;%30e-6;%45e-6;%45 %32
% P_Bragg_src = 5.8;%25;%4.8; %4.2; %power in mW 7 to 9 works well 4.93
K_Bragg_src_1_up= 0.3;%0.14;%0.28;%0.43;%0.75;%0.6;%0.15;%0.0975;%0.4;%0.23;%0.2;%0.5;%0.8;%1;%0.8;%0.62;%0.2;%0.13+0.035;%ampfun([60.117, 0.5638],P_Bragg_src)/2e3;
K_Bragg_src_2_up= 0.3;%0.14;%0.28;%0.43;%0.75;%0.6;%0.15;%0.0975;%0.4;%0.23;%0.2;%0.5;%0.08;%1;%0.8;%0.62;%0.8;%0.1+0.035;%ampfun([132.62, 0.5283],P_Bragg_src)/2e3; %change
Gs_mod_Bragg_src_1_up=5e-6;%6.5e-6;%9e-6;%8.129e-6;%1.3*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01); %1.82
Gs_mod_Bragg_src_2_up=5e-6;%6.5e-6;%8e-6;%9e-6;%6.65e-6;%1.3*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01); %1.81

t0_Bragg_src_up=nan;
wf_Bragg_up_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%
wf_Bragg_src_up_pulse_1 = @(b,t) wf_Bragg_up_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_Bragg_src_up_pulse_2 = @(b,t) wf_Bragg_up_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%%
%%% Bragg splitting: |k=0> |--> |k=-1K> + |k=-2K>
dopp_shft_1 = 38e3;
dF_Bragg_1=(84.96e3+dopp_shft_1)/2; %0.064e6;%0.0685e6;%0.085e6;%0.089e6;
dF_Bragg_2=(84.96e3+dopp_shft_1)/2; %0.064e6;%0.0685e6;%0.085e6;%0.089e6;
f1_Bragg_src_b=f0_AOM-dF_Bragg_1;
f2_Bragg_src_b=f0_AOM+dF_Bragg_2;

dopp_shft_2 = 208e3;
dF_Bragg_1=(84.96e3+dopp_shft_2)/2; %0.064e6+Ek;%0.0685e6;%0.085e6;%0.089e6;
dF_Bragg_2=(84.96e3+dopp_shft_2)/2; %0.064e6+Ek;%0.0685e6;%0.085e6;%0.089e6;
f1_Bragg_src_b_2=f0_AOM-dF_Bragg_1;
f2_Bragg_src_b_2=f0_AOM+dF_Bragg_2;

T_Bragg_src_b=50e-6;%22.5e-6;%34e-6;%3.9e-6;%4.949e-6;%10.2e-6;
% P_Bragg_src = 14.5;%15;%8.9;%4.0; %9
as_1=+0.01;
K_Bragg_src_b_1=0.64 %0.26+as_1;%ampfun([60.117, 0.5638],P_Bragg_src)/2e3;%[132.62, 0.5283]
K_Bragg_src_b_2=0.64 %0.24+as_1;%ampfun([132.62, 0.5283],P_Bragg_src)/2e3;%[60.117, 0.5638]

% P_Bragg_src = 4.55;%6.5; %9
as=0.04;
K_Bragg_src_b_3=0.55; %0.165+as;%ampfun([60.117, 0.5638],P_Bragg_src)/2e3;%ampfun([60.117, 0.5638],P_Bragg_src)/2e3;%0.0;%
K_Bragg_src_b_4=0.55;%0.36%0.135+as;%ampfun([132.62, 0.5283],P_Bragg_src)/2e3;%ampfun([132.62, 0.5283],P_Bragg_src)/2e3;%0.0;%

sg_1 = 0.25e-6;
Gs_mod_Bragg_src_b_1= 8.8e-6;%7.129e-6-sg_1;%3.242*1.72*T_Bragg_src_b/34e-6;%0.8*sqrt(2);%3.2*T_Bragg_src_b/10.2e-6;%1.5*T_Bragg_src_b/11.5e-6*sqrt(2)*sqrt(5.63806937736142e-01);
Gs_mod_Bragg_src_b_2= 8.8e-6;%5.65e-6-sg_1;%3.1393*1.81*T_Bragg_src_b/34e-6;%0.6*sqrt(2);%3.5*T_Bragg_src_b/10.2e-6;%1.5*T_Bragg_src_b/11.5e-6*sqrt(2)*sqrt(5.28341254744861e-01);

sg_2 = -0.2e-6;
Gs_mod_Bragg_src_b_3= 8.2e-6;%6.129e-6-sg_2;%3.242;%2.4*T_Bragg_src_b/11.5e-6;%1.5*T_Bragg_src_b/11.5e-6*sqrt(2)*sqrt(5.63806937736142e-01);
Gs_mod_Bragg_src_b_4= 8.2e-6;%4.65e-6-sg_2;%3.1393;%2.65*T_Bragg_src_b/11.5e-6;%1.5*T_Bragg_src_b/11.5e-6*sqrt(2)*sqrt(5.28341254744861e-01);

t0_Bragg_src_b=nan;%3.82e-6;%nan

wf_Bragg_src_b_pulse_1 = wf_Bragg_src_t_pulse_1;
wf_Bragg_src_b_pulse_2 = wf_Bragg_src_t_pulse_1;

%%
%%% Bragg splitting: |k=0> |--> |k=0> + |k=-2K>

dF_Bragg_1=0.119e6;%0.096e6;%0.099e6;%~0.093;% %%
dF_Bragg_2=0.119e6;%0.096e6;%0.099e6; %%
f1_Bragg_src_l=f0_AOM-dF_Bragg_1;
f2_Bragg_src_l=f0_AOM+dF_Bragg_2;

T_Bragg_src_l=45e-6;%31.5E-6; %%
K_Bragg_src_l_1=0.229;%0.343;%0.4;%0.314;%0.34;
K_Bragg_src_l_2=0.371;%0.286;%0.229;%0.314;%0.25;
Gs_mod_Bragg_src_l_1=8.129e-6;%7.843e-6;%3.9e-6;%6.1e-6;
Gs_mod_Bragg_src_l_2=6.65e-6;%7.15e-6;%3.9e-6;%6.1e-6;
t0_Bragg_src_f=nan;

wf_bragg_src_l_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);%.*cos(pi*(t-b(1)/2)/(b(1))).^2;
wf_bragg_src_l_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_bragg_src_l_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%% WAVEFORM BUILDER
%--------------------------------------------------------------------------
% Waveform Formats:
%{'sine',freq(hz),phase(rad),amplitude,Gauss mod,sample rate,duration}
%{'const',amplitude,sample rate,duration}
% {'double_sine',freq1(hz),freq2(hz),phase1(rad),phase2(rad),amplitude1,...
%     amplitude2,Gauss mod1,Gauss mod2,sample rate,duration}
%--------------------------------------------------------------------------
srate_all=srate_max;
ch1_raw={}; %waveform for chanel 1
ch2_raw={}; %waveform for chanel 2
for ii = 1:length(sequence) %run through each segment
    segment = sequence{ii}; %the current segment
    switch segment %add to the waveforms the desired segment
        case 'const'
            %%% A constant delay
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay}}%{'arb',nullfun, srate_all,T_delay, 0}}
                ];
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay}}%{{'arb',nullfun, srate_all,T_delay, 0}}
                ];
        case 'mag_transfer'
            %%% Magentic transfer
            ch2_raw=[ch2_raw(:)',...
{{'arb',nullfun, srate_all,T_Raman_mix+abs(T_pulse_del), 0}}
%                 {{'const',0, srate_all,T_Raman_mix+abs(T_pulse_del)}}
                ];
            ch1_raw=[ch1_raw(:)',...                                         
{{'double_sine',f1_Raman_mix,f2_Raman_mix,phi1_mix,phi2,K_R_mix,...
                K_R_mix,Gs_mod_R_mix,Gs_mod_R_mix,srate_all,T_Raman_mix,T_pulse_del}}
                ];

        case 'mag_k=-1'
            %%% Top Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix_1}},...
                {{'arb',  wf_Bragg_src_t_pulse_1,  srate_all,  T_Raman_mix,  f1_Raman_mix,  K_R_mix,  Gs_mod_R_mix, phi1_mix,0}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix_1}},...
                {{'arb',  wf_Bragg_src_t_pulse_2,  srate_all,  T_Raman_mix,  f2_Raman_mix,  K_R_mix,  Gs_mod_R_mix, phi2,0}}
                ];

        case 'mag_transfer_in_trap'
            %%% Magentic transfer
            ch1_raw=[ch1_raw(:)',...
                {{'double_sine',f1_Raman_mix_in_trap,f2_Raman_mix_in_trap,phi1_mix,phi2,K_R_mix,...
                K_R_mix,Gs_mod_R_mix,Gs_mod_R_mix,srate_all,T_Raman_mix,T_pulse_del}}
                ];
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_Raman_mix+abs(T_pulse_del)}}
                ];
            case 'k=+1,0,-1'
            %%% Full Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_bragg_sym_pulse_1,  srate_all,  T_Bragg_sym_f,  f1_Bragg_sym_f,  K_Bragg_sym_f_1,  Gs_mod_Bragg_sym_f_1, phi1_Bragg,0}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_bragg_sym_pulse_2,  srate_all,  T_Bragg_sym_f,  f2_Bragg_sym_f,  K_Bragg_sym_f_2,  Gs_mod_Bragg_sym_f_2, phi2_Bragg,0}}
                ];
        case 'k=0,-1,-2'
            %%% Full Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_bragg_src_pulse_1,  srate_all,  T_Bragg_src_f,  f1_Bragg_src_f,  K_Bragg_src_f_1,  Gs_mod_Bragg_src_f_1, phi1_Bragg,0}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_bragg_src_pulse_2,  srate_all,  T_Bragg_src_f,  f2_Bragg_src_f,  K_Bragg_src_f_2,  Gs_mod_Bragg_src_f_2, phi2_Bragg,0}}
                ];
        case 'k=0,-2'
            %%% Large Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_bragg_src_l_pulse_1,  srate_all,  T_Bragg_src_l,  f1_Bragg_src_l,  K_Bragg_src_l_1,  Gs_mod_Bragg_src_l_1, phi1_Bragg,0}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_bragg_src_l_pulse_2,  srate_all,  T_Bragg_src_l,  f2_Bragg_src_l,  K_Bragg_src_l_2,  Gs_mod_Bragg_src_l_2, phi2_Bragg,0}}
                ];
        case 'k=0,-1'
            %%% Top Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',  wf_Bragg_src_t_pulse_1,  srate_all,  T_Bragg_src_t,  f1_Bragg_src_t,  K_Bragg_src_1,  Gs_mod_Bragg_src_1, phi1_Bragg,0}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',  wf_Bragg_src_t_pulse_2,  srate_all,  T_Bragg_src_t,  f2_Bragg_src_t,  K_Bragg_src_2,  Gs_mod_Bragg_src_2, phi2_Bragg,0}}
                ];

        case 'k=+1,0'
            %%% Top Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_Bragg_src_up_pulse_1,  srate_all,  T_Bragg_src_up,  f1_Bragg_src_up,  K_Bragg_src_1_up,  Gs_mod_Bragg_src_1_up, phi1_Bragg,0}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_Bragg_src_up_pulse_2,  srate_all,  T_Bragg_src_up,  f2_Bragg_src_up,  K_Bragg_src_2_up,  Gs_mod_Bragg_src_2_up, phi2_Bragg,0}}
                ];
        case 'k=-1,-2'
            %%% Bottom Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_Bragg_src_b_pulse_1,  srate_all, T_Bragg_src_b,   f1_Bragg_src_b,          K_Bragg_src_b_1,       Gs_mod_Bragg_src_b_1, phi1_Bragg,0}},...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_Bragg_src_b_pulse_1,  srate_all, T_Bragg_src_b,   f1_Bragg_src_b_2,          K_Bragg_src_b_3,       Gs_mod_Bragg_src_b_3, phi1_Bragg,0}},...
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_Bragg_src_b_pulse_2,  srate_all,   T_Bragg_src_b    f2_Bragg_src_b,          K_Bragg_src_b_2,       Gs_mod_Bragg_src_b_2, phi2_Bragg, 0}},...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_Bragg_src_b_pulse_2,  srate_all,   T_Bragg_src_b    f2_Bragg_src_b_2,          K_Bragg_src_b_4,       Gs_mod_Bragg_src_b_4, phi2_Bragg, 0}},...
                ];
        case 'mirror'
            %%% Mirror pulse
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mirror}},...
                {{'arb',   wf_mirror_pulse_1,  srate_all,  T_Bragg_mirror,  f1_Bragg_mirror,  Amp_sinc_Bragg_mirror_1,  sinc_scale_Bragg_mirror_1, phi1_mirror,  Chirp_grad_mirror}}
                ];
            %             {'sine',    f1_Bragg_mirror       ,phi1,          K_Bragg_mirror_1,       Gs_mod_Bragg_mirror_1, srate_all,   T_Bragg_mirror, t0_Bragg_mirror}
            %
            %
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mirror}},...
                {{'arb',   wf_mirror_pulse_2,  srate_all,   T_Bragg_mirror, f2_Bragg_mirror,  Amp_sinc_Bragg_mirror_2,  sinc_scale_Bragg_mirror_2, phi2_mirror,  Chirp_grad_mirror}}
                ];
            %             {'sine',    f2_Bragg_mirror       ,phi2,          K_Bragg_mirror_2,       Gs_mod_Bragg_mirror_2, srate_all,   T_Bragg_mirror, t0_Bragg_mirror}
            %
            %
        case 'mirror2'
            %%% Mirror pulse
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mirror_2}},...
                {{'arb',   wf_mirror_pulse_1,  srate_all,  T_Bragg_mirror,  f1_Bragg_mirror,  Amp_sinc_Bragg_mirror_1,  sinc_scale_Bragg_mirror_1, phi1_mirror,  Chirp_grad_mirror}}
                ];
            %             {'sine',    f1_Bragg_mirror       ,phi1,          K_Bragg_mirror_1,       Gs_mod_Bragg_mirror_1, srate_all,   T_Bragg_mirror, t0_Bragg_mirror}
            %
            %
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mirror_2}},...
                {{'arb',   wf_mirror_pulse_2,  srate_all,   T_Bragg_mirror, f2_Bragg_mirror,  Amp_sinc_Bragg_mirror_2,  sinc_scale_Bragg_mirror_2, phi2_mirror,  Chirp_grad_mirror}}
                ];
            %             {'sine',    f2_Bragg_mirror       ,phi2,          K_Bragg_mirror_2,       Gs_mod_Bragg_mirror_2, srate_all,   T_Bragg_mirror, t0_Bragg_mirror}
            %
            %
        case 'splitter'
            %%% 50:50 Beam splitter pulse %                 {{'arb',ampfun, srate_all,T_delay_splitter, 0}}
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_splitter}},...
                {{'arb',   wf_splitter_pulse_1,  srate_all,  T_Bragg_splitter,  f1_Bragg_splitter,  Amp_sinc_Bragg_splitter_1,  sinc_scale_Bragg_splitter_1, phi1_splitter,Chirp_grad_splitter}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_splitter}},...
                {{'arb',   wf_splitter_pulse_2,  srate_all,   T_Bragg_splitter, f2_Bragg_splitter,  Amp_sinc_Bragg_splitter_2,  sinc_scale_Bragg_splitter_2, phi2_splitter,Chirp_grad_splitter}}
                ];
        case 'splitter_MZ'
            %%% 50:50 Beam splitter pulse
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_splitter_MZ}},...
                {{'arb',   wf_splitter_pulse_1,  srate_all,  T_Bragg_splitter,  f1_Bragg_splitter_MZ,  Amp_sinc_Bragg_splitter_1_MZ,  sinc_scale_Bragg_splitter_1, phi1_splitter_MZ,  Chirp_grad_splitter}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_splitter_MZ}},...
                {{'arb',   wf_splitter_pulse_2,  srate_all,   T_Bragg_splitter, f2_Bragg_splitter_MZ,  Amp_sinc_Bragg_splitter_2_MZ,  sinc_scale_Bragg_splitter_2, phi2_splitter_MZ,  Chirp_grad_splitter}}
                ]; 
        otherwise
            error('invalid sequence segment');
    end
end

%% Package and sent to instrument
chanels_dev1={ch_to_waveforms(ch1_raw),ch_to_waveforms(ch2_raw)};
%%
plot_segments(chanels_dev1,1);
%%
send_segments(chanels_dev1,1);