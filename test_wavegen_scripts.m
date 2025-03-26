%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interface between labview and matlab

%INPUT FROM LABVIEW
%   i          - integer,itteration number
%   auto_enable- boolean
%   file       - string : tells program where to look for setting files
%   file_exact - string
%   mloop      - boolean(legacy) do you want to interface with mloop or use matlab to scan over variables
%   control    - boolean(legacy)

%return
%exact_line- double (legacy just return a zero)
%new_path (the settings that labview will use for the next shot)

%% USER SETTINGS

update_keysight = 1; %do you want to update keysight
addaptive_evap = 1;%do you want to adaptively change the evap setpt to keep atom num constant
do_mag = 0; %mesure the magnetic transfer
do_num = 0; %measure the total atom number

num_points = 6;%6;
shots_per_point = 3;
phi_lim = pi;%2*pi;%

param_lim = phi_lim;%0.2;%12;%

mag_interval = 150;%how often to check the magnetic transfer percentage
num_interval = 300;% how often to check the total number

% if i == 1
%     shot_offset = 0;
%     next_shot_num = 0;
%     shot_mask = []; %mask for adaptive halo monitoring
% end
% 
% % marker = mod(floor((i-1-shot_offset)/shots_per_point),num_points)+1; %Counts from 1 to num shots before setpt update
% 
% if mod((i-1-shot_offset),shots_per_point) == 0 %alernative marker method where point is chosen randomly
%     marker = randi([1 num_points],1);
%     if (marker == 1 || marker == 5 || marker == 6) && (rand(1) > 0.45)
%         marker = randi([2 4],1);
%     end
% %     X = 1:num_points;
% %     P = [0.25, 0.05, 0.1, 0.05, 0.25, 0.05, 0.1, 0.05, 0.1];
% %     marker = randsample( X, 1, true, P );
% end

%main path for generating halos
new_path='c:\remote\settings202314Mar102319.xml';%'c:\remote\settings202315Mar143112.xml';

%path for measuring magnetic transfer efficency
mag_path='c:\remote\settings202127Jul091627.xml';%evap: ~0.842 MHz
%path for measuring the total atom number
num_path='c:\remote\settings202127Jul093144.xml';%evap: ~0.854 MHz

% 'c:\remote\settings202120Jul142856.xml';%evap: ~0.845 MHz
trap_type='quad_0_7_shunt_0_75';%'quad_3_0_shunt_0_9';


%% Sequence selection
% NOTE: choose one sequence from below
% NOTE: K is equal to the momentum recoil of the laser detuning
%--------------------------------------------------------------------------
% Segments:
%'mag_transfer':  Magnetic transfer
%'k=0,-1,-2':     delay + Bragg (|k=0> |--> |k=0> + |k=-1K> + |k=-2K>)
%'k=+1,0,-1':     delay + Bragg (|k=0> |--> |k=+1> + |k=0K> + |k=-1K>)
%'k=0,-1':        delay + Bragg (|k=0> |--> |k=0> + |k=-1K>)
%'k=+1,0':        delay + Bragg (|k=0> |--> |k=+1> + |k=0K>)
%'k=-1,-2':       delay + Bragg (|k=0> |--> |k=-1K> + |k=-2K>)
%'mirror':        delay + Mirror pulse
%'splitter':      delay + 50:50 Beam Splitter
%'const':         delay
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Common Sequences:
% Mach-Zender:    {'k=0,-1','splitter','mirror','splitter_MZ','mag_transfer'}
% Rarity-Tapster: {'k=+1,0,-1','mirror','splitter','const','mag_transfer'}

% sequence = {'k=0,-1,-2','mirror','splitter_MZ','mag_transfer'}; %Construct desired experimental sequenc from sequences above

sequence = {'k=+1,0,-1','mirror','splitter','const','mag_transfer'};%, RT interferometer sequence
% sequence = {'k=0,-1','mirror','const','mag_transfer'};% check mirror
% sequence = {'k=0,-1','splitter','const','mag_transfer'};% check splitter
% sequence = {'const','mag_transfer'};% check mag transfer

tag = 'main';
%% Adaptive number control
addpath('C:\Users\BEC Machine\cloudstor\PROJECTS\Momentum_Bells_test\dev')
top_target = 12.5;
update_interval = 5;
% if addaptive_evap
%     if i == 1 || ~exist('evap_setting','var')
%         %set intial evap height
%         evap_setting = 10;%15;%8;%10
%     end
%     [new_path, evap_setting, num_path]=evap_setting_update(evap_setting,i,update_interval,trap_type,top_target,shot_mask);
% end

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

%% General configs for Bragg
f0_AOM=80e6;            % [Hz]     Central AOMs frequency
Ek=84.9e3;              % beam geometry (90 deg)
ampfun = @(b,x) b(1).*x(:,1).^b(2);
ampfun_sin = @(b,x) asin(sqrt(x(:,1)./b(1)))./b(2);
nullfun = @(b,x) b(1).*b(2).*x(:,1).*0;

%% PULSE SETTINGS

%%% DELAYs between pulses
%--------------------------------------------------------------------------
T_delay_mix=0;         % Delay between the MAG and Bragg pulse
T_delay_mirror=210e-6; %Delay between the Bragg pulse and Mirror pulse
T_delay_mirror_2=5e-6; %Delay between the Bragg pulse and Mirror pulse
T_delay_splitter=210e-6;  %Delay between the Bragg pulse and Mirror pulse
T_delay_splitter_MZ=200e-6;  
T_delay=1420e-6-470e-6;%Const Delay


%%% phases
%--------------------------------------------------------------------------
% NOTE: separated out since we find zero gives fine results
phi1_Bragg=0;
phi2_Bragg=0;

phi1_mirror=0;
phi2_mirror=0;

phi1_splitter=1.571;%2.81;%-57.9415;
phi2_splitter=0;

phi1_splitter_MZ=0;%3*pi/4;
phi2_splitter_MZ=0;

%%% MAGNETIC TRANSFER pulse
%--------------------------------------------------------------------------
% quad 3.4, shunt 0.3
% B_trap_bottom=1.096e6+60e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;%
% B_trap_bottom_intrap=4.2e6;%1.025e6;%1.0603782e6;
% del = -4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
% del2 = 0e3;
% T_pulse_del = 0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% dF_Raman_in_trap=-(B_trap_bottom_intrap);     %[Hz]    Raman detuning
% T_Raman_mix=26e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
% Gs_mod_R_mix=1.35;%0.015;%0.8;%1.8;%2.69400894;
% phi1_mix=pi;
% phi2=0;
% K_R_mix=0.22;%0.4;%0.338;%0.338;%0.34063277;%0.31
%
% f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
% f2_Raman_mix=f0_AOM+dF_Raman/2-del2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"
%
% f1_Raman_mix_in_trap=f0_AOM-dF_Raman_in_trap/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
% f2_Raman_mix_in_trap=f0_AOM+dF_Raman_in_trap/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

% quad 3.0, shunt 0.9 at t=0
%--------------------------------------------------------------------------
% B_trap_bottom=1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;%
% B_trap_bottom_intrap=4.2e6;%1.025e6;%1.0603782e6;
% del = -4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
% del2 = 0e3;
% T_pulse_del = 0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% dF_Raman_in_trap=-(B_trap_bottom_intrap);     %[Hz]    Raman detuning
% T_Raman_mix=25e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
% Gs_mod_R_mix=0.95;%0.015;%0.8;%1.8;%2.69400894;
% phi1_mix=pi;
% phi2=0;
% K_R_mix=0.2;%0.4;%0.338;%0.338;%0.34063277;%0.31
%
% f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
% f2_Raman_mix=f0_AOM+dF_Raman/2-del2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"
%
% f1_Raman_mix_in_trap=f0_AOM-dF_Raman_in_trap/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
% f2_Raman_mix_in_trap=f0_AOM+dF_Raman_in_trap/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"
%

% quad 3.0, shunt 0.9 at t=360e-6 s
%--------------------------------------------------------------------------
B_trap_bottom=1.096e6-0e3;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;%
B_trap_bottom_intrap=4.2e6;%1.025e6;%1.0603782e6;
del = -7e3;%-4.8e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
del2 = 0e3;
T_pulse_del = 0.0e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
dF_Raman_in_trap=-(B_trap_bottom_intrap);     %[Hz]    Raman detuning
T_Raman_mix=25e-6;%25e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
Gs_mod_R_mix=1.1;%0.95;%0.015;%0.8;%1.8;%2.69400894;
phi1_mix=pi;
phi2=0;
K_R_mix=0.21;%0.2;%0.4;%0.338;%0.338;%0.34063277;%0.31

f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2-del2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

f1_Raman_mix_in_trap=f0_AOM-dF_Raman_in_trap/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix_in_trap=f0_AOM+dF_Raman_in_trap/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"


% quad 0.7, shunt 0.75 at t=1420e-6 s
%--------------------------------------------------------------------------
B_trap_bottom=0.929e6;%1.096e6+70e3-0e3;%1.025e6;%1.0603782e6; 1.0851e6;%0.8495e6;%
del = -11e3;%-10e3;%15e3;%-45e3;%-41.03543351e3; %detuning for second beam 3e3
del2 = 0e3;
T_pulse_del = 0.0e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
T_Raman_mix=45e-6;%25.5e-6;%24e-6;%23e-6;%19.20368817e-6;
Gs_mod_R_mix=0.7;%0.015;%0.8;%1.8;%2.69400894;
phi1_mix=pi;
phi2=0;
K_R_mix=0.12;%0.095;%%0.095;%0.4;%0.338;%0.338;%0.34063277;%0.31

f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2-del2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"


%%% MIRROR pulse (halo at 200 mus delay centered for k=+1,0,-1 halos)
%--------------------------------------------------------------------------
dF_Bragg_1=18e3;%
dF_Bragg_2=18e3;%

f1_Bragg_mirror=f0_AOM-dF_Bragg_1;
f2_Bragg_mirror=f0_AOM+dF_Bragg_2;

T_Bragg_mirror=40E-6;
t0_Bragg_mirror=nan;

sinc_scale_Bragg_mirror_1=4e-6;
sinc_scale_Bragg_mirror_2=4e-6;

Amp_sinc_Bragg_mirror_1=0.28;
Amp_sinc_Bragg_mirror_2=0.28;

Chirp_grad_mirror = 0;

% wf_mirror_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);
wf_mirror_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);

wf_mirror_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_mirror_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%%% MIRROR pulse (bec at 200 mus delay)
%--------------------------------------------------------------------------

% dF_Bragg_1=0.064e6;%17e3;
% dF_Bragg_2=0.064e6;%17e3;
% f1_Bragg_mirror=f0_AOM-dF_Bragg_1;
% f2_Bragg_mirror=f0_AOM+dF_Bragg_2;
%
% T_Bragg_mirror=40E-6;%10E-6;%35E-6;%15E-6;%32E-6;
% t0_Bragg_mirror=nan;%3.9895e-6;
%
% sinc_scale_Bragg_mirror_1=6e-6;%4e-6;
% sinc_scale_Bragg_mirror_2=6e-6;%4e-6;
%
%
% Amp_sinc_Bragg_mirror_1=0.17;%0.232;
% Amp_sinc_Bragg_mirror_2=0.17;%0.232;
%
% Chirp_grad_mirror = 0;%7.6e6;%
%
% % wf_mirror_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);
% wf_mirror_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);
%
% wf_mirror_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
% wf_mirror_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%%% 50:50 Beam Splitter pulse (halo 400 mus delay)
%--------------------------------------------------------------------------
dF_Bragg_1=18e3;%17
dF_Bragg_2=18e3;%17
f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
f2_Bragg_splitter=f0_AOM+dF_Bragg_2;

T_Bragg_splitter=40E-6;%10E-6;%35E-6;
t0_Bragg_splitter=nan;

sinc_scale_Bragg_splitter_1=4e-6;
sinc_scale_Bragg_splitter_2=4e-6;

Amp_sinc_Bragg_splitter_1=0.18;
Amp_sinc_Bragg_splitter_2=0.18;

Chirp_grad_splitter = 0;%7.6e6;

% wf_splitter_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);
wf_mirror_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%

wf_splitter_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_splitter_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

dF_Bragg_1=0.100e6;%0.114e6;%0.064e6*2;%0.112e6;%0.114429e6;%0.114e6;
dF_Bragg_2=0.100e6;%0.114e6;%0.064e6*2;%0.112e6;%0.114429e6;%0.114e6;
f1_Bragg_splitter_MZ=f0_AOM-dF_Bragg_1;
f2_Bragg_splitter_MZ=f0_AOM+dF_Bragg_2;

Amp_sinc_Bragg_splitter_1_MZ=0.133+0.006;%sqrt(13.0+a_shift);
Amp_sinc_Bragg_splitter_2_MZ=0.133+0.006;%sqrt(13.0+a_shift);

%%% 50:50 Beam Splitter pulse
%--------------------------------------------------------------------------
% dF_Bragg_1=0.064e6;%15e3;
% %0.097e6;%0.112e6;%0.064e6*2;%
% dF_Bragg_2=0.064e6;%15e3;
% %0.097e6;%0.112e6;%0.064e6*2;%0.112e6;
% f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
% f2_Bragg_splitter=f0_AOM+dF_Bragg_2;
%
% dF_Bragg_1=0.100e6;%0.114e6;%0.064e6*2;%0.112e6;%0.114429e6;%0.114e6;
% dF_Bragg_2=0.100e6;%0.114e6;%0.064e6*2;%0.112e6;%0.114429e6;%0.114e6;
% f1_Bragg_splitter_MZ=f0_AOM-dF_Bragg_1;
% f2_Bragg_splitter_MZ=f0_AOM+dF_Bragg_2;
%
% T_Bragg_splitter=40E-6;%10E-6;%35E-6;
% t0_Bragg_splitter=nan;
%
% sinc_scale_Bragg_splitter_1=6e-6;%4e-6;
% %7.286e-6;%8.129e-6;%6.1e-6;
% sinc_scale_Bragg_splitter_2=6e-6;%4e-6;
% %7.286e-6;%6.65e-6;%6.1e-6;
%
% a_shift = +0.0;
%
% Amp_sinc_Bragg_splitter_1=0.12;%0.157
% %0.16;%sqrt(13.0+a_shift);
% Amp_sinc_Bragg_splitter_2=0.12;%0.157
% %sqrt(13.0+a_shift);
%
% Amp_sinc_Bragg_splitter_1_MZ=0.133+0.006;%sqrt(13.0+a_shift);
% Amp_sinc_Bragg_splitter_2_MZ=0.133+0.006;%sqrt(13.0+a_shift);
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

%%% Momentum splitting
%--------------------------------------------------------------------------

%%% Bragg splitting: |k=0> |--> |k=+1K> + |k=0> + |k=-1K>
dF_Bragg_1=19e3;
dF_Bragg_2=19e3;
f1_Bragg_sym_f=f0_AOM-dF_Bragg_1;
f2_Bragg_sym_f=f0_AOM+dF_Bragg_2;

T_Bragg_sym_f=40e-6;
K_Bragg_sym_f_1=0.45;
K_Bragg_sym_f_2=0.45;

Gs_mod_Bragg_sym_f_1=3.3e-6;
Gs_mod_Bragg_sym_f_2=3.3e-6;
t0_Bragg_sym_f=nan;

wf_bragg_sym_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);
wf_bragg_sym_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_bragg_sym_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%%% Bragg splitting: |k=0> |--> |k=0> + |k=-1K> + |k=-2K>
dF_Bragg_1=0.09e6;
dF_Bragg_2=0.09e6;
f1_Bragg_src_f=f0_AOM-dF_Bragg_1;
f2_Bragg_src_f=f0_AOM+dF_Bragg_2;

T_Bragg_src_f=50e-6;
K_Bragg_src_f_1=0.65;
K_Bragg_src_f_2=0.65;

Gs_mod_Bragg_src_f_1=2.8e-6;
Gs_mod_Bragg_src_f_2=2.8e-6;
t0_Bragg_src_f=nan;

wf_bragg_src_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);%.*cos(pi*(t-b(1)/2)/(b(1))).^2;
wf_bragg_src_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_bragg_src_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%%% Bragg splitting: |k=0> |--> |k=0> + |k=-1K>
dF_Bragg_1=64e3;
dF_Bragg_2=64e3;
f1_Bragg_src_t=f0_AOM-dF_Bragg_1;
f2_Bragg_src_t=f0_AOM+dF_Bragg_2;

T_Bragg_src_t=40e-6;
K_Bragg_src_1=0.13;
K_Bragg_src_2=0.13;
Gs_mod_Bragg_src_1=5e-6;
Gs_mod_Bragg_src_2=5e-6;

t0_Bragg_src_t=nan;
wf_Bragg_t_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%
wf_Bragg_src_t_pulse_1 = @(b,t) wf_Bragg_t_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_Bragg_src_t_pulse_2 = @(b,t) wf_Bragg_t_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%%% Bragg splitting: |k=0> |--> |k=0> + |k=+1K>
dF_Bragg_1=64e3-Ek;
dF_Bragg_2=64e3-Ek;%
f1_Bragg_src_up=f0_AOM-dF_Bragg_1;
f2_Bragg_src_up=f0_AOM+dF_Bragg_2;

T_Bragg_src_up=40e-6;
K_Bragg_src_1_up=0.14;
K_Bragg_src_2_up=0.14;
Gs_mod_Bragg_src_1_up=5e-6;
Gs_mod_Bragg_src_2_up=5e-6;

t0_Bragg_src_up=nan;
wf_Bragg_up_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%
wf_Bragg_src_up_pulse_1 = @(b,t) wf_Bragg_up_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_Bragg_src_up_pulse_2 = @(b,t) wf_Bragg_up_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));


%%% Bragg splitting: |k=0> |--> |k=-1K> + |k=-2K>
dF_Bragg_1=0.064e6;
dF_Bragg_2=0.064e6;
f1_Bragg_src_b=f0_AOM-dF_Bragg_1;
f2_Bragg_src_b=f0_AOM+dF_Bragg_2;

dF_Bragg_1=0.064e6+Ek;
dF_Bragg_2=0.064e6+Ek;
f1_Bragg_src_b_2=f0_AOM-dF_Bragg_1;
f2_Bragg_src_b_2=f0_AOM+dF_Bragg_2;

T_Bragg_src_b=40e-6;

as_1=+0.01;
K_Bragg_src_b_1=0.26+as_1;
K_Bragg_src_b_2=0.24+as_1;

as=0.04;
K_Bragg_src_b_3=0.165+as;
K_Bragg_src_b_4=0.135+as;

sg_1 = 0.25e-6;
Gs_mod_Bragg_src_b_1=7.129e-6-sg_1;
Gs_mod_Bragg_src_b_2=5.65e-6-sg_1;

sg_2 = -0.2e-6;
Gs_mod_Bragg_src_b_3=6.129e-6-sg_2;
Gs_mod_Bragg_src_b_4=4.65e-6-sg_2;

t0_Bragg_src_b=nan;

wf_Bragg_src_b_pulse_1 = wf_Bragg_src_t_pulse_1;
wf_Bragg_src_b_pulse_2 = wf_Bragg_src_t_pulse_1;

%%% Bragg splitting: |k=0> |--> |k=0> + |k=-2K>

dF_Bragg_1=0.119e6;
dF_Bragg_2=0.119e6;
f1_Bragg_src_l=f0_AOM-dF_Bragg_1;
f2_Bragg_src_l=f0_AOM+dF_Bragg_2;

T_Bragg_src_l=45e-6;%
K_Bragg_src_l_1=0.229;%
K_Bragg_src_l_2=0.371;%;

Gs_mod_Bragg_src_l_1=8.129e-6;%
Gs_mod_Bragg_src_l_2=6.65e-6;%
t0_Bragg_src_f=nan;

wf_bragg_src_l_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);%.*cos(pi*(t-b(1)/2)/(b(1))).^2;
wf_bragg_src_l_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_bragg_src_l_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));

%% Iteration through parameters
% num_points_sub = sqrt(num_points);%num_points;%
% param_vec = linspace(0.2,0.4,num_points_sub);
% param_vec_2 = linspace(0.15,param_lim,num_points_sub);%[5,6.75,7.75];%

% %phi curve
% param_vec = [linspace(0,pi,5),1.17,2];%linspace(0.0,param_lim,num_points);%[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2;%[pi/4,pi/2];%linspace(0.0,param_lim,num_points);%[1.33,4.5];%[1.78,4.82,8,11.04];%
% if num_points == marker
%     Amp_sinc_Bragg_mirror_1 = 0;%param_vec(marker);%0.18+0.02;%0.154;%0.144;%0.329;%sqrt(25.0);
%     Amp_sinc_Bragg_mirror_2 = 0;%param_vec(marker);%0.18+0.02;%0.199;0.189;%0.307; %sqrt(25.0);5.714, 12.000
%
%     Amp_sinc_Bragg_splitter_1 = 0;%param_vec(marker);%0.143;%0.16;%sqrt(13.0+a_shift);
%     Amp_sinc_Bragg_splitter_2 = 0;%param_vec(marker);%0.143;%sqrt(13.0+a_shift);
%     phi1_splitter = -1;
%     param_vec(marker) = -1;
% else
% %     param_vec(marker) = 1;
%     phi1_splitter = param_vec(marker);
% end

%
% param_vec = [0,0.23,0.26,0.28,0.29];%[0,0.232];%[0,0.05,0.08,0.12,0.14,0.15,0.16,0.195];%[0,linspace(0.16,0.28,num_points-1)];%[0,0.388,0.5];%[0,0.52,0.55];%[0,linspace(0.35,0.65,num_points-1)];%
% 
% Amp_sinc_Bragg_mirror_1=param_vec(marker);%0.18+0.02;%0.154;%0.144;%0.329;%sqrt(25.0);
% Amp_sinc_Bragg_mirror_2=param_vec(marker);%0.18+0.02;%0.199;0.189;%0.307; %sqrt(25.0);5.714, 12.000

% %transmission curve
% param_vec = [0,0.1,0.125,0.155,0.17,0.19];%[0,5*pi/8,7*pi/8,9*pi/8,-1,pi/2];%[0,0.05,0.08,0.12,0.14,0.15,0.16,0.195];%[0,linspace(0.1,0.29,num_points-1)];%[0,linspace(0.16,0.28,num_points-1)];%[0,0.388,0.5];%[0,0.52,0.55];%[0,linspace(0.35,0.65,num_points-1)];%
% 
% %beam splitter settings for 5pi/16 pulse
% dF_Bragg_1=18e3;%17
% dF_Bragg_2=18e3;%17
% f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
% f2_Bragg_splitter=f0_AOM+dF_Bragg_2;
% 
% sinc_scale_Bragg_splitter_1=4e-6;
% sinc_scale_Bragg_splitter_2=4e-6;
% 
% %Amp_sinc_Bragg_splitter_1=0.155;
% %Amp_sinc_Bragg_splitter_2=0.155;
% 
% 
% Amp_sinc_Bragg_splitter_1=param_vec(marker);%0.143;%0.16;%sqrt(13.0+a_shift);
% Amp_sinc_Bragg_splitter_2=param_vec(marker);%0.143;%sqrt(13.0+a_shift);

param_vec = [0,5*pi/8,7*pi/8,9*pi/8,-1,pi/2];%

if marker == 1 %0
    Amp_sinc_Bragg_splitter_1 = 0;%param_vec(marker);%0.143;%0.16;%sqrt(13.0+a_shift);
    Amp_sinc_Bragg_splitter_2 = 0;%param_vec(marker);%0.143;%sqrt(13.0+a_shift);
elseif marker == 2 % 5pi/16 pulse
    dF_Bragg_1=18e3;%17
    dF_Bragg_2=18e3;%17
    f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
    f2_Bragg_splitter=f0_AOM+dF_Bragg_2;

    sinc_scale_Bragg_splitter_1=4e-6;
    sinc_scale_Bragg_splitter_2=4e-6;

    Amp_sinc_Bragg_splitter_1=0.16;
    Amp_sinc_Bragg_splitter_2=0.16;
elseif marker == 3 %cross pulse 7pi/8
    dF_Bragg_1=10e3;%17
    dF_Bragg_2=10e3;%17
    f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
    f2_Bragg_splitter=f0_AOM+dF_Bragg_2;

    sinc_scale_Bragg_splitter_1=7e-6;
    sinc_scale_Bragg_splitter_2=7e-6;

    Amp_sinc_Bragg_splitter_1=0.105;
    Amp_sinc_Bragg_splitter_2=0.105;
elseif marker == 4 %9pi/16 pulse
    dF_Bragg_1=19e3;%17
    dF_Bragg_2=19e3;%17
    f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
    f2_Bragg_splitter=f0_AOM+dF_Bragg_2;

    sinc_scale_Bragg_splitter_1=4e-6;
    sinc_scale_Bragg_splitter_2=4e-6;

    Amp_sinc_Bragg_splitter_1=0.085;
    Amp_sinc_Bragg_splitter_2=0.085;
elseif marker == 5
    Amp_sinc_Bragg_mirror_1 = 0;%param_vec(marker);%0.18+0.02;%0.154;%0.144;%0.329;%sqrt(25.0);
    Amp_sinc_Bragg_mirror_2 = 0;%param_vec(marker);%0.18+0.02;%0.199;0.189;%0.307; %sqrt(25.0);5.714, 12.000

    Amp_sinc_Bragg_splitter_1 = 0;%param_vec(marker);%0.143;%0.16;%sqrt(13.0+a_shift);
    Amp_sinc_Bragg_splitter_2 = 0;%param_vec(marker);%0.143;%sqrt(13.0+a_shift);

    param_vec(marker) = -1;
end


% param_vec = linspace(0,20,num_points).*1e3;
%
% f1_Bragg_mirror=f0_AOM-param_vec(marker);
% f2_Bragg_mirror=f0_AOM+param_vec(marker);

% T_delay_mirror=200e-6+40e-6;
% if mod(floor(i/100),3) == 0
%     T_delay_mirror=200e-6+40e-6;
% elseif mod(floor(i/100),3) == 1
%     T_delay_mirror=0e-6;
% else
%     T_delay_mirror=400e-6+80e-6-2/sqrt(pi)*95e-6;
% end

% phi1_splitter = param_vec(marker);
% phi2_splitter = param_vec(marker);
% phi1_mirror = param_vec(marker);

% param_vec = linspace(-0.1,0.1,num_points);
%
% a_shift = param_vec(marker);
%
% Amp_sinc_Bragg_splitter_1=sqrt(13.0+a_shift);
% Amp_sinc_Bragg_splitter_2=sqrt(13.0+a_shift);
%
% Amp_sinc_Bragg_splitter_1_MZ=sqrt(13.0+a_shift);
% Amp_sinc_Bragg_splitter_2_MZ=sqrt(13.0+a_shift);

% param_vec = linspace(113,115,num_points);
%
% dF_Bragg_1=param_vec(marker)*1e3;%0.114e6;
% dF_Bragg_2=param_vec(marker)*1e3;%0.114e6;
% f1_Bragg_splitter_MZ=f0_AOM-dF_Bragg_1;
% f2_Bragg_splitter_MZ=f0_AOM+dF_Bragg_2;

% param_vec = repmat(param_vec,1,num_points_sub);
% param_vec_2 = repelem(param_vec_2,num_points_sub);
%
% Amp_sinc_Bragg_mirror_1=param_vec(marker);
% Amp_sinc_Bragg_mirror_2=param_vec_2(marker);

% sinc_scale_Bragg_mirror_1=param_vec(marker).*1e-6;%3.9e-6;%6.1e-6;
% sinc_scale_Bragg_mirror_2=param_vec_2(marker).*1e-6;%3.9e-6;%6.1e-6;

% Gs_mod_Bragg_src_f_1=param_vec(marker).*1e-6;%7.843e-6;%3.9e-6;%6.1e-6;
% Gs_mod_Bragg_src_f_2=param_vec_2(marker).*1e-6;%7.15e-6;%3.9e-6;%6.1e-6;

% dF_Bragg_1=param_vec(marker).*1e3;%0.112e6;
% dF_Bragg_2=param_vec_2(marker).*1e3;%0.112e6;
% % f1_Bragg_mirror=f0_AOM-dF_Bragg_1;
% % f2_Bragg_mirror=f0_AOM+dF_Bragg_2;
% f1_Bragg_src_f=f0_AOM-dF_Bragg_1;
% f2_Bragg_src_f=f0_AOM+dF_Bragg_2;

% K_Bragg_src_f_1=param_vec(marker);
% K_Bragg_src_f_2=param_vec_2(marker);

% K_Bragg_src_l_1=param_vec(marker);
% K_Bragg_src_l_2=param_vec_2(marker);

%% Shot type
%add ability to switch between different shot types
if next_shot_num
    new_path = num_path;
    sequence = {'empty'}; %Construct desired experimental sequenc from sequences above
    tag = 'atom_num';
    update_keysight = 0;
    shot_offset = shot_offset + 1;
    next_shot_num = 0;
    shot_mask = [shot_mask,0];
else
    if do_mag && mod(i,mag_interval) == 0 && do_num && mod(i,num_interval) == 0
        new_path = mag_path;
        sequence = {'const','mag_transfer'}; %Construct desired experimental sequenc from sequences above
        tag = 'mag_transfer';
        shot_offset = shot_offset + 1;
        next_shot_num = 1;
        T_delay = T_delay_mirror+T_delay_splitter_MZ+125e-6;
        shot_mask = [shot_mask,0];
    elseif do_mag && mod(i,mag_interval) == 0
        new_path = mag_path;
        sequence = {'const','mag_transfer'}; %Construct desired experimental sequenc from sequences above
        tag = 'mag_transfer';
        shot_offset = shot_offset + 1;
        T_delay = T_delay_mirror+T_delay_splitter_MZ+125e-6;
        next_shot_num = 0;
        shot_mask = [shot_mask,0];
    elseif do_num && mod(i,num_interval) == 0
        new_path = num_path;
        sequence = {'empty'}; %Construct desired experimental sequenc from sequences above
        tag = 'atom_num';
        update_keysight = 0;
        next_shot_num = 0;
        shot_offset = shot_offset + 1;
        shot_mask = [shot_mask,0];
    else
        shot_mask = [shot_mask,1];
    end
end
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
            ch1_raw=[ch1_raw(:)',...
                {{'arb',nullfun, srate_all,T_Raman_mix+abs(T_pulse_del), 0}}
                %                 {{'const',0, srate_all,T_Raman_mix+abs(T_pulse_del)}}
                ];
            ch2_raw=[ch2_raw(:)',...
                {{'double_sine',f1_Raman_mix,f2_Raman_mix,phi1_mix,phi2,K_R_mix,...
                K_R_mix,Gs_mod_R_mix,Gs_mod_R_mix,srate_all,T_Raman_mix,T_pulse_del}}
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
                {{'arb',   wf_Bragg_src_t_pulse_1,  srate_all,  T_Bragg_src_t,  f1_Bragg_src_t,  K_Bragg_src_1,  Gs_mod_Bragg_src_1, phi1_Bragg,0}}
                ];

            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_Bragg_src_t_pulse_2,  srate_all,  T_Bragg_src_t,  f2_Bragg_src_t,  K_Bragg_src_2,  Gs_mod_Bragg_src_2, phi2_Bragg,0}}
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

%% Shot sequence settings
path_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_LabviewMatlab.txt';
path_param_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_KeysightMatlab.txt';
path_phi = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_Phi.txt';

%% convert waveforms to printables for logging
ch1_waveform_str = '';
ch2_waveform_str = '';
for waveforms = 1:numel(ch1_raw)
    if waveforms>1
        ch1_waveform_str = [ch1_waveform_str,', '];
    end
    ch1_waveform_str = [ch1_waveform_str,cell2str(array2str(ch1_raw{waveforms}))];
end
for waveforms = 1:numel(ch2_raw)
    if waveforms>1
        ch2_waveform_str = [ch2_waveform_str,', '];
    end
    ch2_waveform_str = [ch2_waveform_str,cell2str(array2str(ch2_raw{waveforms}))];
end

ch1_waveform_str = replace(ch1_waveform_str,"'",'');
ch2_waveform_str = replace(ch2_waveform_str,"'",'');

%% Interface

addpath('C:\Users\BEC Machine\cloudstor\PROJECTS\keysight-33600a')
% shot_info = shots.(shot_sequence{marker});
% new_path=shot_info.LVfile;
% Send waveforms
chanels_dev1={ch_to_waveforms(ch1_raw),ch_to_waveforms(ch2_raw)};
if update_keysight && mod((i-1-shot_offset),shots_per_point) == 0 %|| (mod(i,mag_interval) == 0 || mod(i,mag_interval) == 1))
    send_segments(chanels_dev1,1);
    plot_segments(chanels_dev1,1);
end
%write to log
%write to log
f_log=fopen(path_log,'a');  % append to log-file
nowdt=datetime('now');
fprintf(f_log,'shot num:%d, posixtime:%.3f, date:%s, matlab:interfacev8, labview settings:%s, %s\n',...
    i,posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),new_path,tag);
fclose(f_log);

f_log=fopen(path_param_log,'a');  % append to param-log-file
fprintf(f_log,'shot num:%d, posixtime:%.3f, ch1 waveform: %s, ch2 waveform: %s\n',...
    i,posixtime(nowdt),...
    ch1_waveform_str,ch2_waveform_str);
fclose(f_log);

if strcmp(tag,'main')
    f_log=fopen(path_phi,'a');  % append to param-log-file
    fprintf(f_log,'%d, %.3f, %.3f, %.3f\n',...
        i,posixtime(nowdt),...
        param_vec(marker),T_delay_mirror.*1e6);%param_vec(marker)
    fclose(f_log);
end

pause(0.1)
