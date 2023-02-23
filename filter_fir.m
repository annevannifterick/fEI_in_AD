function [Data_filtered] = filter_fir(Data,hp,lp,fs,fir_order);
%
% Modified klaus.linkenkaer@cncr.vu.nl, 070521.
% Modified S. Poil, 12-Feb-2008, added warning for wrong filter order
%
%% NEXT version
% #1  time should be along 1st dimension to improve memory access. 
% 
%******************************************************************************************************************
% Purpose...
%
% Causal (feedforward) bandpass filter the signal 'F' with a Hamming window.
%
%
%******************************************************************************************************************
% Input parameters...
%
% Data          : data matrix (or vector), time along the 2nd dimension!
% hp            : highpass corner (e.g., 8 Hz).
% lp            : lowpass corner (e.g., 13 Hz).
% fs            : sampling frequency.
% fir_order     : filter order in units of seconds (to prevent an influence of sampling frequency!)
%
%******************************************************************************************************************
% Default parameters (can be changed)...

% time window included in the filtering process.
% Filter orders suitable for alpha and beta oscillations and based on:
% Nikulin. 2005. Neurosci. Long-range temporal correlations in EEG oscillations...

%fir_order = 0.22;      % Use for high time resolution and low frequency resolution.
%fir_order = 0.38       % Use for low time resolution and high frequency resolution.

%**************************************************************************
% Warn if Data might have wrong dimensions
if(size(Data,2) > size(Data,1))
    warning('Filter_fir:Dimension','Filter_fir::Input migth have wrong dimensions. Time should be along the 1st dimension!')
end
% Warn if filter order is too low
if(2/hp > fir_order)
    warning('Filter_fir:FilterOrder','Filter_fir:: The filter order is too low for the given high-pass frequency.')
    disp('Use minimum')
    disp(2/hp)
end

%******************************************************************************************************************
% Define filter characteristics:
%%Finite impulse response filter of fir_order*fs order, and [2hp/fs 2lp/fs] Hamming window
b = fir1(floor(fir_order*fs),[hp lp]/(fs/2));

%******************************************************************************************************************
% Filter the vector 'F' using the filter characteristics from above:

Data_filtered = zeros(size(Data));
for ii = 1:size(Data,2)
    Data_filtered(:,ii) = filter(b,1,Data(:,ii));
end