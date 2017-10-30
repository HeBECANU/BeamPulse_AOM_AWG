% (c) Roman Khakimov, 2015
function instr_update(param)
% INSTR_UPDATE(PARAM)
%
% assign WaveformGenMain waveform parameters from input PARAM
% NOTE: this function must be configured by user in each scenario
% Remember to remove variable declaration in WaveformGenMain.m
%

global DIR_WATCH_ENBL

%     T_flight=param
%     K_R=param
%     Trigg_Delay=param
    %B_trap_bottom=param(1)
%     K_R=param(1)
%     %K_R=param(1)
%     T_Raman=param(2)
%     Gs_mod_R=param(3)

%% 2017-10-27: Bell v2, mixing beam characterisation for mf=0
% K_R_mix=param;

%% 2017-10-29: Bell v2
% par_value is a N-by-2 array of [MF, K_R_MIX]
% - spin-pol halos between mf with no rotation pulse to characterise source with noise
% - mf=1 only: more raman amp data
source_mf=param(1);
K_R_mix=param(2);

%% generate waveform
%     WaveformGenMainv4;
WaveformGenMain;
end