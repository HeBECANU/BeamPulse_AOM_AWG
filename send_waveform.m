function send_waveform(ch1_raw,ch2_raw,plot_opt)
%takes in a waveform and sends it to the keysight
chanels_dev1={ch_to_waveforms(ch1_raw),ch_to_waveforms(ch2_raw)};
if plot_opt
    plot_segments(chanels_dev1,1)
end
send_segments(chanels_dev1,1)
end