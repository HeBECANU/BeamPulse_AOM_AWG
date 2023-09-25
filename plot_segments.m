function plot_segments(chanels,zeros)
%this plots the wavefroms defined by channels
%zeros=false;
%chanels=chanels_dev1;
%to do
%add a time accululator to handle variable sample rates (if the gen
%supports this)
%titles ect

%dt = 1/1e9;
% figure(10)
% for n=1:size(chanels,2)
%     chan_amp=[]; %point list
%     chan_time=[];
%     for m=1:size(chanels{n},2)
%         sub_wf=chanels{n}.waveform;
%         sr=chanels{n}.sr;
%         dt = 1/sr;
%         if zeros
%             chan_amp=[chan_amp,repmat(chanels{n}(m).waveform,1,chanels{n}(m).repeats)];
%         else
%             chan_amp=[chan_amp,chanels{n}(m).waveform];
%         end
%     end
%     fprintf('total points',size(chan_amp,2))
%     subplot(2,1,n)
%     t=(0:(length(chan_amp)-1)).*dt*1e6;
%     plot(t,chan_amp)
%     xlabel('time (\(\mu\)s)','interpreter','latex')
%     ylabel('V','interpreter','latex')
%     title(['Channel ',num2str(n)])
%     
% end
% clear chan_amp
% pause(0.1)

figure(10)
for n=1:size(chanels,2)
    chan_amp=[]; %point list
    chan_time=[];
    for m=1:size(chanels{n},2)
        sub_wf=chanels{n}.waveform;
        sr=chanels{n}.sr;
        dt = 1/sr;
        if zeros
            chan_amp=[chan_amp,repmat(chanels{n}(m).waveform,1,chanels{n}(m).repeats)];
        else
            chan_amp=[chan_amp,chanels{n}(m).waveform];
        end
    end
%     fprintf('total points',size(chan_amp,2))
    subplot(2,1,n)
    t=(0:(length(chan_amp)-1)).*dt*1e6;
    plot(t,chan_amp)
    xlabel('time (\(\mu\)s)','interpreter','latex')
    ylabel('V','interpreter','latex')
    title(['Channel ',num2str(n)])
    
end
clear chan_amp
pause(0.1)
% 
% g1=chanels{1}(4).waveform;
% g2=chanels{2}(5).waveform;
% sinc2=chanels{2}(8).waveform;
% sinc1=chanels{1}(7).waveform;

end
