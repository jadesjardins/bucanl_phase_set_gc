function EEG=eegPSGC(EEG,dattype,datind,evttype,winSize,refractPeriod,PSThresh,freqCenter,freqWidth,pntLikely,MConsid,fnameout)
% % Global parameters for the method

datind=str2num(datind);

global Data
switch dattype;
    case 1
        Data=double(EEG.data(datind,:));
    case 2
        Data=EEG.icaact(datind,:);
end

global boundInd
boundInd=ceil([EEG.event(strmatch('boundary',{EEG.event.type},'exact')).latency]);

global NumChan
NumChan=size(Data,1);

%W is the size of the windows in the regression
%RW is the size of the refractionary period between shifts
%2ms per point, 50 points is 100ms window
global W RW
W=winSize;%50;
RW=refractPeriod;%5;

%Phase shifting threshold
global c
c = PSThresh;%0.02;

%frequency band of interest 8-10Hz (w +/- width)
global w width
w=freqCenter;%9;
width=freqWidth;%1;

% Number of points (deltaP) and actual time (delta) in 
% likelihood approximation
global delta deltaP
deltaP=pntLikely;%2;
delta=pntLikely/EEG.srate;%2/500;

% Sampling Rate  
global T
T = EEG.srate;%500;


% The values of M that are considered
Mvec=MConsid;%3:6;

%SaveString is for the output, same as the base of the Batchfile
SaveString=fnameout;


PSGCplugin_FindModel();
PSGCplugin_FindResults();
