function sig_with_ramp = aux_add_ramps(sig_in, fs, dur_ramp, type_ramp)
% sig_with_ramp = aux_add_ramps(sig_in, fs, dur_ramp, type_ramp)
%
% This function applies ramps of the type defined by 'type_ramp':
% 'blackman', 'cosine' or 'hanning' to the on onset and offset of the
% given signal 'sig_in'. The ramp duration is defined by 'dur_ramp'.

%% Add onset and offset ramps of dur_ramp ms to a signal

% Ramp in samples
ramp_samples        = round(dur_ramp*fs);
stim_samples        = length(sig_in);

% Select ramp type
switch type_ramp
    case 'blackman'
        ramp_full   = blackman(2*ramp_samples);
        ramp_full   = ramp_full(:)';
        ramp_asc    = ramp_full(1:ramp_samples);
        ramp_desc   = ramp_full(ramp_samples+1:end);
        
    case 'cosine'
        ramp_full   = sin(linspace(0, pi, 2*ramp_samples));
        ramp_full   = ramp_full(:)';
        ramp_asc    = ramp_full(1:ramp_samples);
        ramp_desc   = ramp_full(ramp_samples+1:end);
        
       
    case 'hanning'
        ramp_full   = hanning(2*ramp_samples);
        ramp_full   = ramp_full(:)';
        ramp_asc    = ramp_full(1:ramp_samples);
        ramp_desc   = ramp_full(ramp_samples+1:end);
end

% Generate ramp shape
n_ones              = stim_samples - length(ramp_full);

ramp_shape          = [ramp_asc, ones(1, n_ones), ramp_desc];
ramp_shape          = ramp_shape(:)';

% Apply ramp to the input signal
sig_with_ramp       = sig_in .* ramp_shape;

end