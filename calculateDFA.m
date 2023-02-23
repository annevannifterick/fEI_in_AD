% Originally created by Richard Hardstone (2020), rhardstone@gmail.com

function [DFAExponent,DFAIntercept,meanDF,windowSizes] = calculateDFA(Signal, ComputeInterval, FitInterval, sampling_frequency, windowOverlap)

% Signal dimensions (numSamples,numChannels)
% CalcInterval in seconds
% windowOverlap in samples

% check that the intervals are properly set
assert(FitInterval(1)>=ComputeInterval(1) & FitInterval(2)<=ComputeInterval(2),'CalcInterval should be included in ComputeInterval');
assert(ComputeInterval(1)>=0.1 & ComputeInterval(2)<=1000,'ComputeInterval should be between 0.1 and 1000 seconds');

%compute DFA window sizes for the given CalcInterval
windowSizes = floor(logspace(-1,3,40) * sampling_frequency);%logspace from 0.1 seccond (10^-1) to 1000 (10^3) seconds

windowSizes = windowSizes(windowSizes >= ComputeInterval(1)*sampling_frequency & ...
    windowSizes <= ComputeInterval(2)*sampling_frequency);

fitIntervalFirstWindow = find(windowSizes>=FitInterval(1)*sampling_frequency,1,'first');
fitIntervalLastWindow = find(windowSizes<=FitInterval(2)*sampling_frequency,1,'last');
fitIntervalLength = fitIntervalLastWindow - fitIntervalFirstWindow + 1;

lengthSignal = size(Signal,1);
numChannels = size(Signal,2);
meanDF = zeros(numChannels,length(windowSizes));
DFAExponent = zeros(numChannels,1);
DFAIntercept = zeros(numChannels,1);
windowSizes = windowSizes(:);

for i_channel = 1:size(Signal,2)
    for i_windowSize = 1:length(windowSizes)
        windowOffset = floor(windowSizes(i_windowSize) * (1-windowOverlap));
        allWindowIndex = createWindowIndices(lengthSignal, windowSizes(i_windowSize), windowOffset);
        originalAmplitude = Signal(:,i_channel);
        signalProfile = cumsum(originalAmplitude - mean(originalAmplitude));    
        xSignal = signalProfile(allWindowIndex);
        dSignal = detrend(xSignal');
        w_detrendedFluctuations = std(dSignal,1); %std(dSignal,1) means normalized by N instead of N-1
        meanDF(i_channel,i_windowSize) = mean(w_detrendedFluctuations); 
    end
   
    
    if fitIntervalLength > 1
        X = [ones(fitIntervalLength,1) log10(windowSizes(fitIntervalFirstWindow:fitIntervalLastWindow))];
        Y = log10(meanDF(i_channel,fitIntervalFirstWindow:fitIntervalLastWindow))';
        regressOutput = regress(Y,X);
        DFAExponent(i_channel) = regressOutput(2,1);
        DFAIntercept(i_channel) = regressOutput(1,1);
    else
        DFAExponent(i_channel) = nan;
        DFAIntercept(i_channel) = nan;
    end
end

function allWindowIndex = createWindowIndices(lengthSignal, lengthWindow, windowOffset)

windowStarts = (1:windowOffset:lengthSignal-lengthWindow+1)-1;
numWindows = length(windowStarts);

oneWindowIndex = 1:lengthWindow;
allWindowIndex = repmat(oneWindowIndex,[numWindows,1]);

allWindowIndex = allWindowIndex + repmat(windowStarts',[1 lengthWindow]);
