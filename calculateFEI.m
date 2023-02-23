% Originally created by Richard Hardstone (2020), rhardstone@gmail.com
% Richard Hardstone and Klaus Linkenkaer-Hansen have filed the patent claim (PCT/NL2019/050167) “Method of determining brain activity”; with priority date 16 March 2018.

function [fEI, fEI_outliers_removed, wAmp, wDNF] = calculateFEI(Signal, windowSize, windowOverlap, DFAexponent)

% Signal dimensions (numSamples,numChannels)
% windowSize in samples
% windowOverlap is fraction of overlap between windows (0-1)
% DFA exponent value for all channels of the signal

lengthSignal = size(Signal,1);
numChannels = size(Signal,2);

windowOffset = floor(windowSize * (1-windowOverlap));
allWindowIndex = createWindowIndices(lengthSignal, windowSize, windowOffset);
numWindows = size(allWindowIndex,1);

fEI = zeros(numChannels,1);
wAmp = zeros(numChannels,numWindows);
wDNF = zeros(numChannels,numWindows);

for i_channel = 1:numChannels
    originalAmplitude = Signal(:,i_channel);    
    signalProfile = cumsum(originalAmplitude - mean(originalAmplitude));            %% Step ii -> iii       
    w_originalAmplitude = mean(originalAmplitude(allWindowIndex),2); %Calculate mean amplitude for each window
    xAmp = repmat(w_originalAmplitude,[1 windowSize]); 
    xSignal = signalProfile(allWindowIndex); %Arrange Signals into windows
    xSignal = (xSignal ./ xAmp)';                                                   %% Step iii -> iv       
    dSignal = detrend(xSignal);                                                     %% Step iv -> v       
    w_detrendedNormalizedFluctuations = std(dSignal,1)';                            %% Step v -> vi       
    fEI(i_channel)  = 1-corr(w_detrendedNormalizedFluctuations,w_originalAmplitude); %% Step vi -> vii

    %outlier removal
    amp_outliers = isoutlier(w_originalAmplitude,'gesd');
    fluctuations_outliers = isoutlier(w_detrendedNormalizedFluctuations,'gesd');
    not_outliers = ~amp_outliers & ~fluctuations_outliers;
    fEI_outliers_removed(i_channel) = 1-corr(w_detrendedNormalizedFluctuations(not_outliers),w_originalAmplitude(not_outliers));

    wAmp(i_channel,:) = w_originalAmplitude;
    wDNF(i_channel,:) = w_detrendedNormalizedFluctuations;
end

% set values to nan for all channels that do not have long-range temporal
% correlations (DFA<=0.6)
fEI(DFAexponent<=0.6) = nan;
fEI_outliers_removed(DFAexponent<=0.6) = nan;

function allWindowIndex = createWindowIndices(lengthSignal, lengthWindow, windowOffset)

windowStarts = (1:windowOffset:lengthSignal-lengthWindow+1)-1;
numWindows = length(windowStarts);

oneWindowIndex = 1:lengthWindow;
allWindowIndex = repmat(oneWindowIndex,[numWindows,1]);

allWindowIndex = allWindowIndex + repmat(windowStarts',[1 lengthWindow]);