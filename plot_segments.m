function plot_segments(chanels,zeros)
%this plots the wavefroms defined by channels
%zeros=false;
%chanels=chanels_dev1;
%to do
    %add a time accululator to handle variable sample rates (if the gen
    %supports this)
    %titles ect


    figure(10)
	for n=1:2
    	chan_amp=[]; %point list
        chan_time=[];
        	for m=1:size(chanels{n},2)
            	sub_wf=chanels{n}.waveform;
                if zeros
                    chan_amp=[chan_amp,repmat(chanels{n}(m).waveform,1,chanels{n}(m).repeats)];
                else
                    chan_amp=[chan_amp,chanels{n}(m).waveform];
                end
           end
        fprintf('total points',size(chan_amp,2))   
        subplot(2,1,n)
        plot(chan_amp)

    end
    clear chan_amp
    pause(0.1)
end
