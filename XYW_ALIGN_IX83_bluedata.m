
function XYW_ALIGN_IX83(CS)
%function XYW_ALIGN_IX83(CS)
%
%for n=1:12, XYW_ALIGN_IX83(n); end

clear global;
global RED
RootName = ['CS' num2str(CS) '_'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALIGNMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    load([RootName '_dRXY.mat']);
    disp([RootName ' dRXY alignment previously done']);
catch
    D{1} = dir([RootName '*GFP*.vsi']);   %green data
    D{2} = dir([RootName '*DAPI*.vsi']); %blue data, ddHT
    %ENSURE THAT FILES ARE IN PROPER ORDER
    for q=1:2
        CaptureNum=[];
        for k=1:length(D{q})
            tmp = NumbersInString(D{q}(k).name);
            CaptureNum(k) = tmp(2);
        end
        [dum,I] = sort(CaptureNum);
        D{q} = D{q}(I);
    end
    %CHECK
    if length(D{2}) == 2*length(D{1})
        disp('Aha, you used the new 2 x mTRITC protocol. Nice!');
        %9/2018 -- two mTRITC exposures (5ms=new and 20ms=old)
        D{3} = D{2}(1:2:end);
        D{2} = D{2}(2:2:end); %original exposure time
    end
    if length(D{1}) ~= length(D{2})
        error('mismatch in number of mGFP vs mTRITC files');
    end
    r=0;
    for k=1:length(D{1})
        disp(D{2}(k).name);
        Rtmp = bfopenIX83(D{2}(k).name); %load the next mTRITC data
        if k==1
            NumRuns = length(D{1}); %number of runs
            
            RFPReps = size(Rtmp,3); %number of repeitions for a single run
            nTR = RFPReps * NumRuns; %total # timepoints
            RED = zeros(size(Rtmp,1),size(Rtmp,2),nTR, class(Rtmp)); %allocate memory
        end
        RED(:,:,r+1:r+RFPReps) = Rtmp;  %copy
        r = r + RFPReps;
    end
    clear Rtmp;  %clean up
    
    
    h = AlignImageGUI_IX83(0);
    waitfor(h);
    dRXY = evalin('base','AIH.dRXY');
    evalin('base','clear AIH');
    save([RootName '_dRXY.mat'], 'dRXY');
end

