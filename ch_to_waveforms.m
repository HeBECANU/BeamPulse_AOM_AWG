function ch_processed=ch_to_waveforms(ch_raw)
%this function takes an array of waveform specifications and returns a (cell array/or structure) of wavefroms
%this code deals with constants by looping the minimum number of points (32) in order to conserve
%wavefrom memory and allow for longer sequences
%[type('sine','const'),freq(hz),phase(rad),duration,amplitude,sample rate,Gauss mod]
global max_points  points_min repeats_max

points_min=double(256);%double(16);%double(32);%double(32);

%to improve
    %if the number of repeats exceeds the max then make the sequence larger
    %so it fits

n_points_tot=0;
seg_index=1; %need a semprate wf index because the constant parts will generate 1 or 2
zero_at_end=1;
ch_processed=[];

for n=1:size(ch_raw,2)
    
    seg_raw=ch_raw{n};
    
    seg_type=seg_raw{1};
    
 
    
    if strcmp(seg_type,'sine')
        seg_freq=seg_raw{2};
        seg_phase=seg_raw{3};
        seg_amp=seg_raw{4};
        seg_gmod=seg_raw{5};
        seg_sr=seg_raw{6}; %sample rate
        seg_duration=seg_raw{7};
        dt=1/seg_sr;
        n_points=round(seg_duration/dt);
        
        L = (n_points)-1;  % codegeneration special case: Ensure order is scalar
        
        if length(seg_raw)>7 && ~isnan(seg_raw{8})
            seg_center = seg_raw{8}/dt;
        else
            seg_center = L/2;
        end
        
        if n_points>max_points
            error('Requested segment length is too long\n')
        end
        
        t=linspace(0, n_points*dt, n_points);   
        if ~isnan(seg_gmod)
%             envelope=seg_amp*gausswin(n_points,seg_gmod).';
%             envelope=seg_amp*chebwin(n_points,seg_gmod*50).';
        npts = (0:(n_points-1))'-seg_center;
        envelope=seg_amp.* exp(-(1/2).*(seg_gmod.*npts./(L/2)).^2);
        else
            envelope = seg_amp(ones(1,n_points));
        end
        wf_out=sin(2*pi*seg_freq*t+seg_phase);    
        wf_out=wf_out.*envelope';
        PA=trapz(t, envelope);
        if zero_at_end & n==size(ch_raw,2)
            wf_out=[wf_out,0];
        end
        
		fprintf('Pulse_area=%2.3e for duration %3.1f µs\n', PA, seg_duration*1e6)
        
        ch_processed(seg_index).waveform=wf_out;
        
        %due to the intrinsic asymetery that can happen when the waveform
        %is enveloped i add 2% to the pk-pk voltage
        
        ch_processed(seg_index).vpp=(max(wf_out)-min(wf_out))*1.02;
        ch_processed(seg_index).sr=seg_sr;
        ch_processed(seg_index).repeats=1;
        seg_index=seg_index+1;
        n_points_tot=n_points_tot+n_points;
    elseif strcmp(seg_type,'sweep_sine')
        %freq start, freq end, phase, amp, sample rate, duration
        seg_freq_start=seg_raw{2};
        seg_freq_end=seg_raw{3};
        seg_phase=seg_raw{4};
        seg_amp=seg_raw{5};
        seg_sr=seg_raw{6}; %sample rate
        seg_duration=seg_raw{7};
        dt=1/seg_sr;
        n_points=round(seg_duration/dt);
        
        if n_points>max_points
            error('Requested segment length is too long\n')
        end
        
        t=linspace(0, n_points*dt, n_points); 
        freq_grad=(seg_freq_end-seg_freq_start)/seg_duration;
        wf_out=seg_amp*sin(2*pi*(seg_freq_start+freq_grad*t).*t+seg_phase);    
        if zero_at_end & n==size(ch_raw,2)
            wf_out=[wf_out,0];
        end
        ch_processed(seg_index).waveform=wf_out;
        
        %due to the intrinsic asymetery that can happen when the waveform
        %is enveloped i add 2% to the pk-pk voltage
        
        ch_processed(seg_index).vpp=(max(wf_out)-min(wf_out))*1.02;
        ch_processed(seg_index).sr=seg_sr;
        ch_processed(seg_index).repeats=1;
        seg_index=seg_index+1;
        n_points_tot=n_points_tot+n_points;
    elseif strcmp(seg_type,'double_sine')
        seg_freq1=seg_raw{2};
        seg_freq2=seg_raw{3};
        seg_phase1=seg_raw{4};
        seg_phase2=seg_raw{5};
        seg_amp1=seg_raw{6};
        seg_amp2=seg_raw{7};
        seg_gmod1=seg_raw{8};
        seg_gmod2=seg_raw{9};
        seg_sr=seg_raw{10}; %sample rate
        seg_duration=seg_raw{11};
        seg_del=seg_raw{12};%delay between pulses
        
        
        
        dt=1/seg_sr;
        n_points=round(seg_duration/dt);
        n_delayed=abs(round(seg_del/dt));

        L = (n_points)-1;  % codegeneration special case: Ensure order is scalar
        
        if length(seg_raw)>12
            seg_center = seg_raw{13}/dt;
        else
            seg_center = L/2;
        end
        
        if n_points>max_points
            error('Requested segment length is too long\n')
        end
        
         
%         envelope1=seg_amp1*gausswin(n_points,seg_gmod1).';
%         envelope2=seg_amp2*gausswin(n_points,seg_gmod2).';
        
        n1 = (0:(n_delayed+n_points-1))'-seg_center;
        n2 = (0:(n_delayed+n_points-1))'-seg_center-n_delayed;
        envelope1=seg_amp1.* exp(-(1/2).*(seg_gmod1.*n1./(L/2)).^2);
        envelope2=seg_amp2.* exp(-(1/2).*(seg_gmod1.*n2./(L/2)).^2);
%         envelope1=seg_amp1*chebwin(n_points,seg_gmod1*50).';
%         envelope2=seg_amp2*chebwin(n_points,seg_gmod2*50).';
%         envelope1=seg_amp1*tukeywin(n_points,1/seg_gmod1).';
%         envelope2=seg_amp2*tukeywin(n_points,1/seg_gmod2).';
        t=linspace(0, (n_delayed+n_points)*dt, (n_delayed+n_points));
%         t=linspace(0, (n_delayed+n_points)*dt, (n_delayed+n_points));  
        wf_out1=sin(2*pi*seg_freq1*t+seg_phase1);    
        wf_out1=wf_out1.*envelope1';
%         if seg_del>0
%             wf_out1=[wf_out1,zeros(1,n_delayed)];
%         else
%             wf_out1=[zeros(1,n_delayed),wf_out1];
%         end
        
%         t=linspace(0, (n_points)*dt, n_points);  
        wf_out2=sin(2*pi*seg_freq2*t+seg_phase2);    
        wf_out2=wf_out2.*envelope2';
%         if seg_del>0
%             wf_out2=[zeros(1,n_delayed),wf_out2];
%         else
%             wf_out2=[wf_out2,zeros(1,n_delayed)];
%         end
        
        wf_out=wf_out1+wf_out2;
        PA=trapz(t, envelope1)+trapz(t, envelope2);
        if zero_at_end & n==size(ch_raw,2)
            wf_out=[wf_out,0];
        end
        
		fprintf('Pulse_area=%2.3e for duration %3.1f µs\n', PA, seg_duration*1e6)
        
        ch_processed(seg_index).waveform=wf_out;
        
        %due to the intrinsic asymetery that can happen when the waveform
        %is enveloped i add 2% to the pk-pk voltage
        
        ch_processed(seg_index).vpp=(max(wf_out)-min(wf_out))*1.02;
        ch_processed(seg_index).sr=seg_sr;
        ch_processed(seg_index).repeats=1;
        seg_index=seg_index+1;
        n_points_tot=n_points_tot+n_points; 
        

    elseif strcmp(seg_type,'const') 
        seg_amp=seg_raw{2};
        seg_sr=seg_raw{3}; %sample rate
        seg_duration=seg_raw{4};
        

        dt=1/seg_sr;
        n_points=round(seg_duration/dt);
        %here we break the contant spacing into as many points_min blocks
        %as possible with one points_min+excess
        seg_repeats=floor(n_points/points_min-1); 
        points_rem=n_points-(seg_repeats-1)*points_min;
  
        if seg_repeats>repeats_max
            error('Requested repeats for constant exceeds %i \n',repeats_max)
        elseif seg_repeats>0
            ch_processed(seg_index).waveform=ones(1,points_min)*seg_amp;
            ch_processed(seg_index).sr=seg_sr;
            ch_processed(seg_index).repeats=seg_repeats;
            ch_processed(seg_index).vpp=0.001; %set range to min
            seg_index=seg_index+1;
            n_points_tot=n_points_tot+points_min;
        end
        
        if points_rem>0
            ch_processed(seg_index).waveform=ones(1,points_rem)*seg_amp;
            ch_processed(seg_index).sr=seg_sr;
            ch_processed(seg_index).repeats=1;
            ch_processed(seg_index).vpp=0.001;
            seg_index=seg_index+1;
            n_points_tot=n_points_tot+points_rem;
        end    
        
        
        if zero_at_end & n==size(ch_raw,2) %if this is the last segment add a zero at the end
            ch_processed(seg_index-1).waveform=[ch_processed(seg_index-1).waveform,0]
        end
    elseif strcmp(seg_type,'arb')
        %arbitrary waveform
        seg_sr=seg_raw{3};
        seg_duration=seg_raw{4};
        
        params = [seg_raw{4:end}];
        

        dt=1/seg_sr;
        n_points=round(seg_duration/dt);
        wf_func = seg_raw{2};
        
        t=linspace(0, n_points*dt, n_points); 
        wf_out=wf_func(params,t')';
        
        PA=trapz(t, wf_out);
        if zero_at_end && n==size(ch_raw,2)
            wf_out=[wf_out,0];
        end
        
		fprintf('Pulse_area=%2.3e for duration %3.1f µs\n', PA, seg_duration*1e6)
        
        ch_processed(seg_index).waveform=wf_out;
        
        %due to the intrinsic asymetery that can happen when the waveform
        %is enveloped i add 2% to the pk-pk voltage
        
        ch_processed(seg_index).vpp=(max(wf_out)-min(wf_out))*1.02;
        ch_processed(seg_index).sr=seg_sr;
        ch_processed(seg_index).repeats=1;
        seg_index=seg_index+1;
        n_points_tot=n_points_tot+n_points;
    else
        error('Unknown Waveform type\n')    
    end
       
end

%check that the total length is resonable
if n_points_tot>max_points
    error('Requested total waveform length is too long\n')
end

fprintf('done generating waveforms\n')


end
