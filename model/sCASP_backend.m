function [ output] = sCASP_backend(int_rep_speech, int_rep_mix, fs, cf_aud, cf_mod)

 if size(int_rep_speech)~=size(int_rep_mix)
    error('internal representations must have the same size');
end

[Nsamp Naud_ch Nmod_ch] = size(int_rep_mix);

   
WinDurs = 1./cf_mod; % The window duration is the inverse of the centerfrequency of the modulation channel
WinDurs(1) = 1/2.5;

WinLengths = floor(WinDurs * fs);
Nsegments = floor(Nsamp./WinLengths)+ ones(1,length(cf_mod)); % The total number of segments is Nframes plus any additional "leftover"
    
      
if find(WinLengths == Nsamp)% If the duration of the stimulus is exactly equal to the window duration
        segIdx = find(WinLengths == Nsamp);
        Nsegments(segIdx) =  Nsegments(segIdx)-1;
end
    

for m=1:length(cf_mod)
    
    rule4th=find(cf_aud > 4*cf_mod(m)); % Apply rule of cf_mod < 1/4 cf_aud
   
    % Build the modulation 2D matrices:                                       %
    speech=squeeze(int_rep_speech(:, :, m));
    mix=squeeze(int_rep_mix(:, :, m));                                  
                                      
    % Delete discarded bands
     speech=speech(:, rule4th);
     mix=mix(:, rule4th);    
    
    tmp_ssnn = zeros(WinLengths(m), size(mix, 2), Nsegments(m)); % Allocate memory for multi-resolution
    tmp_ss = tmp_ssnn;
        
    segLengths = zeros(1,Nsegments(m)) ;
    
    % Find starting and ending points of the segments:
    
     for i = 1:Nsegments(m) % For each temporal segment of the signal
                               % find the start and end index of the frame
            if i > (Nsegments(m)-1)
                startIdx = 1 + (i-1)*WinLengths(m);
                endIdx = Nsamp;
            else
                startIdx = 1 + (i-1)*WinLengths(m);
                endIdx = startIdx + WinLengths(m)-1;
            end
            
            segment = startIdx:endIdx;
            segLengths(i) = length(segment);
            
         % internal representation of the temporal segments (samplesPerSegment x all bands x  number of segments)
         
         tmp_ss(1:segLengths(i), :, i) = speech(segment,:);
         tmp_ssnn(1:segLengths(i), :, i) = mix(segment,:);
     
        
         dint_temp(i, m) = corr2(tmp_ss(1:segLengths(i), :, i), tmp_ssnn(1:segLengths(i),:, i));
         dint_temp(dint_temp<0)=0; % Remove negative correlations
  
   end
        
 dmod(m)=sum(dint_temp(1:Nsegments(m), m))/Nsegments(m);
     
end


output.dint=dmod;
output.dsegments=dint_temp;
output.dfinal=mean(dmod);

end