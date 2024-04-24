
%NOTE: flumazenil-DART2.0 seems very left-shifted.  
%but data is very noisy, so I'm just guessing
%  HTneg ~1-3 uM vs  (dose groups are all different, but can see trends) 
%  HTpos ~0.01 uM    (very poorly constrained, given no believable baseline)
% so, DW ~ 100 to 300.  
% this is reason enough to favor diazepam-DART over flumazenil
% (not to mention that flumazenil-DART has altered funciton vs flumazenil)


if 0
    clear;
end

try
    CELLS;
catch
    PathList = { ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_08_GABA_PAM_Darts\batch02' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_08_GABA_PAM_Darts\batch03' ...
        ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_15_GABA_PAM_Darts\batch_01' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_15_GABA_PAM_Darts\batch_02' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_16_GABA_PAM_Darts\batch_03' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_16_GABA_PAM_Darts\batch_04' ...
        ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_22_GABA_PAM_Darts\batch_01' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_22_GABA_PAM_Darts\batch_02' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_23_GABA_PAM_Darts_ver2\batch_03' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\DIAZEPAM_DART2.0\2019_05_23_GABA_PAM_Darts_ver2\batch_04' ...
        };
    
    if isempty(dir(PathList{1}))
        %running on laptop -- change root folder
        for k=1:length(PathList)
            PathList{k} = strrep(PathList{k},  'E:\Dropbox (TadrossLab)\BIG_DATA',   'C:\Users\mrt31\Dropbox (TadrossLab)\BIG_DATA');
        end
    end    
    CELLS = XYW_LoadCells(PathList);
    for k=1:length(CELLS)
        %MQ207-DART is inactive, so use these as a proxy for CONTROL
        if strcmpi(CELLS(k).SeriesName, 'MQ207-DART2.0')
            CELLS(k).SeriesName = 'CONTROL';
        end
    end
end


%For subplots -- choose one of these
OPTIONS.REFBINS = 7;
OPTIONS.REFTXT = 'RESPONSEcv';   OPTIONS.REFLIM = [ 0  0.5 ];
OPTIONS.REFTXT = 'log10(delGstd)';   OPTIONS.REFLIM = [ 1  4 ];
OPTIONS.REFTXT = 'mean(grnMin)'; OPTIONS.REFLIM = [ 0  0.7];
OPTIONS.REFTXT = 'mean(grnRaw)'; OPTIONS.REFLIM = [ 0  10000];
OPTIONS.REFTXT = 'mean(grnNorm)'; OPTIONS.REFLIM = [ 0.2  0.5];
OPTIONS.REFTXT = 'FilledArea';  OPTIONS.REFLIM = [ 0 1000 ];
OPTIONS.REFTXT = 'RESP(5)'; OPTIONS.REFLIM = [ 0  1 ];
OPTIONS.REFTXT = 'RESP(1)';     OPTIONS.REFLIM = [ 0.5  1.5 ];
OPTIONS.REFTXT = 'mean(RESPONSEstd)';   OPTIONS.REFLIM = [0  0.5 ];
OPTIONS.REFTXT = 'SHAPE';   OPTIONS.REFLIM = [ 0.95  1 ];

OPTIONS.REFTXT = 'grnMinWavyness'; OPTIONS.REFLIM = [ 0  0.06 ];
OPTIONS.REFTXT = 'RESPONSEstd';   OPTIONS.REFLIM = [ 0  0.5 ];
OPTIONS.REFTXT = 'REF0';        OPTIONS.REFLIM = [ 0  0.5 ];
OPTIONS.REFTXT = 'REF0/REF1';   OPTIONS.REFLIM = [ 0.  1.2 ]; %**THIS IS THE CRITICAL GUIDANCE FOR TUNING OPTIONS.dwREF0
OPTIONS.REFTXT = 'RESP(1)-RESP(6)'; OPTIONS.REFLIM = [ -1  1 ]; %THIS WAS USED TO check that ALLOWABLE RANGE OF RESP(1)-RESP(6) is al OK.  choice of max delta -1:+1 was arbitrary
OPTIONS.REFTXT = 'log10(delG)'; OPTIONS.REFLIM = [3.0  4.5];
OPTIONS.REFTXT = 'REF1-REF0'; OPTIONS.REFLIM = [0  0.6];


OPTIONS.DoseVectPlot = [ 0 3 10 30  100  300  1000  3000  10000  Inf  30000 ]; %additional figures for given doses

OPTIONS.YLIM = [-0.25 1.4];
OPTIONS.bScatterLog = 1;  %for dTomato
OPTIONS.bShowScatterPos = 1;
OPTIONS.bShowScatterNeg = 1;
OPTIONS.bShowLine = 0;
OPTIONS.LegendLocation = 'southwest';

OPTIONS.bUseMedian=0;
OPTIONS.bUseMax = 0;  %use GCaMP max, instead of mean
OPTIONS.MAXGrnClip = inf; %65000;

OPTIONS.NumRepsPerDose = 6;    %THIS IS HARD CODED
OPTIONS.RepVect = 1:6; %use all reps
OPTIONS.MaxNperCS = inf;  %assume every cell is independent, even if on the same CS

OPTIONS.bLinearInterpolateFloor = 2; %use (higher-order) interpolation for the floor

OPTIONS.RLOW = 3.0;
OPTIONS.RHIGH = [3.5  inf];
OPTIONS.RBIN = [-inf 2.5:0.5:4.5 inf];

OPTIONS.nREF1 = [1 6];    %one-normalization
OPTIONS.bFixRESP16 = 1;  %flag to correct for rundown-dependence of middle doses

OPTIONS.nREF0 = [5];  %zero-subtractions
%note: plots of CONTROL Series(1) data as a function of REF0/REF1 showed a strong positive slope
%note: this slope was not present for dose(1) or dose(6), but became steep for intermediate doses
%      I'm talking aobut the control, so there was no drug on board. It's a rundown phenomenon. 
%      Something about rundown is related to impact of diazepam during dose(5)
%note2: (X-REF0)/(REF1-REF0) produced a very positive slope
%       (X-REF0)/REF1        produced a very negative slope
%       (X-REF0)/(REF1-dw*REF0) could be made flat for a chosen dw.
%so, to account for this, I added the following:
if 1
    %fix the issue
    OPTIONS.dwREF0         = [1      0.9   0.7   0.3        1 1];  %the smaller the number, the more it counteracts a positive slope
OPTIONS.RundownCorrection = [1.1        1  0.85  0.73       0.73  0.9]; %[1 1 1 1 1 1];
    OPTIONS.ScalingPostREF = [1      0.88  0.85  0.6       1 1]; %used in concert with dwREF0
elseif 1
    %fix the issue
    OPTIONS.dwREF0         = 1;
OPTIONS.RundownCorrection  = [1.1       1  1  1                 0.73  0.9]; %[1 1 1 1 1 1];
    OPTIONS.ScalingPostREF = [1        0.95  0.8   0.5          1    1]; %used in concert with dwREF0
else
    %default -- useful to adjust RundownCorrection
    OPTIONS.dwREF0         = 1;
    OPTIONS.ScalingPostREF = 1;
OPTIONS.RundownCorrection = [1.1    1  0.88  0.73    0.73  0.9]; %[1 1 1 1 1 1];
end

MM=[-0.5  2.5];
OPTIONS.RespMinMax = { ...
    [1 MM] [2 MM] [3 MM] [4 MM] [5 MM] [6 MM] ...
    [1 -6  -1  1]  ...
    {'delG' 10.^[3.8 5]} ...
    {'REF1-REF0' 0 1} ...
    };  %avoid excessive runup/down (factor of 2)

%    {'delG' 10.^[3.5 5]} {'REF0/REF1' 0 1} ...
%    {'RESPONSEstd' 0 0.25}  {'SHAPE' 0 1.0}...
%    {'grnMinWavyness' 0 0.02} ...

[Series, CELLS] = XYW_DoseResponsePlot2(CELLS, OPTIONS, [1 2 3 ]);

figure(1); try delete(h); catch, end
x = 10.^[-3.5:.01:2.8];
m=1.0; z=x.^m;
n=1.2; y=x.^n;
h=plot(...
    x, 0*x+1, 'k', ...
    x, 1-0.86*z./(z + 0.05^m), 'b', ...
    x, 1-0.86*z./(z + 38^m),   'r', ...
    x, 1-0.81*y./(y + 0.05^n), 'r', ...
    'linewidth', 3);


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check a specific dose
CHECK_SERIES = 2; %length(Series);
CHECK_DOSEGRP = 2;
SORT_DOSE_ =  3;
CHECK_FLD     = 'Ineg';
CELLCHECK = Series(CHECK_SERIES).CELLS(Series(CHECK_SERIES).(CHECK_FLD));
TMP = vertcat(CELLCHECK.RESPONSE);
CELLCHECK = CELLCHECK(~isnan(TMP(:,CHECK_DOSEGRP+2)));
if SORT_DOSE_>0
    TMP = vertcat(CELLCHECK.RESPONSE_);
    [~,I] = sort(TMP(:,SORT_DOSE_), 'descend'); CELLCHECK = CELLCHECK(I);
end
for n=1:length(CELLCHECK)
    figure(21); clf; imagesc(CELLCHECK(n).SHAPECOV); axis image; hold on; set(gca,'ydir','normal', 'xtick', 0:6:36);
    plot(0.5+(1:length(CELLCHECK(n).grnMin))/16, 1+5*CELLCHECK(n).grnMin, '.r', 'linewidth', 2);
    plot(0.5+(1:length(CELLCHECK(n).grnMin))/16, 1+0*CELLCHECK(n).grnMin, '.r', 'linewidth', 2);
    plot(0.5+(1:length(CELLCHECK(n).grnNorm))/16, 1+5*CELLCHECK(n).grnNorm, '.k', 'linewidth', 2);
    title([num2str(n) ' of ' num2str(length(CELLCHECK))  ] ); 
    colorbar; caxis([-1.1  1.1]);
    disp(CELLCHECK(n))
    pause;
end




