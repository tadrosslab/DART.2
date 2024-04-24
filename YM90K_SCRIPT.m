
if 0
    clear;
end

try
    CELLS;
catch
    PathList = { ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_05_YM90Kcomparison_NewProt\Plate1_2018_09_05' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_05_YM90Kcomparison_NewProt\plate2_2018_09_05' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_05_YM90Kcomparison_NewProt\plate3_2018_09_05' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_05_YM90Kcomparison_NewProt\plate4_2018_09_06' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_05_YM90Kcomparison_NewProt\plate5_2018_09_06' ...
        ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_26_YM90K\plate1_2018_09_26' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_26_YM90K\plate2_2018_09_26' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_26_YM90K\plate4_2018_09_28' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_26_YM90K\plate5_2018_09_28' ...
        };
    
    %one extra folder with a camera glitch in some doses
    %PathList{end+1} = 'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\YM90K_DART2.0\2018_09_26_YM90K\plate3_2018_09_26';
    
    if length(dir(PathList{1}))==0
        %running on laptop -- change root folder
        for k=1:length(PathList)
            PathList{k} = strrep(PathList{k},  'E:\Dropbox (TadrossLab)\BIG_DATA',   'C:\Users\mrt31\Dropbox (TadrossLab)\BIG_DATA');
        end
    end
    CELLS = XYW_LoadCells(PathList);
end

  
%For subplots -- choose one of these
OPTIONS.REFBINS = 7;
OPTIONS.REFTXT = 'red1Std';   OPTIONS.REFLIM = [ 0  100 ];
OPTIONS.REFTXT = 'log10(delGstd)';   OPTIONS.REFLIM = [ 1  4 ];
OPTIONS.REFTXT = 'max(RESP)';   OPTIONS.REFLIM = [ 0.9 4 ];
OPTIONS.REFTXT = 'min(RESP)';   OPTIONS.REFLIM = [ -0.25 0.05 ];
OPTIONS.REFTXT = 'RESPONSEcv';   OPTIONS.REFLIM = [ -1  2 ];
OPTIONS.REFTXT = 'REF0/REF1';   OPTIONS.REFLIM = [ 0  1 ];
OPTIONS.REFTXT = 'REF1-REF0'; OPTIONS.REFLIM = [ 0  1 ];
OPTIONS.REFTXT = 'grnMinWavyness'; OPTIONS.REFLIM = [ 0  0.5 ];
OPTIONS.REFTXT = 'RESP(1)';     OPTIONS.REFLIM = [ -0.08  0.16 ];
OPTIONS.REFTXT = 'REF0';        OPTIONS.REFLIM = [ -0.2  1.4 ];
OPTIONS.REFTXT = 'REF1';        OPTIONS.REFLIM = [ 0.25  0.65 ];
OPTIONS.REFTXT = 'min(RESPONSEstd)';   OPTIONS.REFLIM = [0  0.2 ];
OPTIONS.REFTXT = 'max(RESPONSEstd)';   OPTIONS.REFLIM = [0  0.2 ];
OPTIONS.REFTXT = 'mean(RESPONSEstd)';   OPTIONS.REFLIM = [0  0.2 ];
OPTIONS.REFTXT = 'RESPONSEstd';   OPTIONS.REFLIM = [ 0  0.2 ];
OPTIONS.REFTXT = 'RESP(1)-RESP(2)'; OPTIONS.REFLIM = [ -0.5  0.5 ];
OPTIONS.REFTXT = 'mean(grnMin)'; OPTIONS.REFLIM = [ -0.1  1.1];
OPTIONS.REFTXT = 'FilledArea';  OPTIONS.REFLIM = [ 0 1000 ];
OPTIONS.REFTXT = 'log10(delG)'; OPTIONS.REFLIM = [3 4];
OPTIONS.REFTXT = 'SHAPE';   OPTIONS.REFLIM = [ 0  1 ];
OPTIONS.REFTXT = 'RESP(4)'; OPTIONS.REFLIM = [ -0.2  1.6 ];
OPTIONS.REFTXT = 'FloorDelta'; OPTIONS.REFLIM = [ 0  0.5 ];

OPTIONS.DoseVectPlot = [ 1  10  30  100  300  1000 3000  10000  inf ]; %additional figures for given doses

OPTIONS.YLIM = [-0.25 1.4];
OPTIONS.bScatterLog = 1;  %for dTomato
OPTIONS.bShowScatterPos = 1;
OPTIONS.bShowScatterNeg = 1;
OPTIONS.LegendLocation = 'southwest';
OPTIONS.bShowLine = 0;

OPTIONS.bUseMedian=0;
OPTIONS.bUseMax = 0;  %use GCaMP max, instead of mean
OPTIONS.MAXGrnClip = inf; %65000;

OPTIONS.NumRepsPerDose = 6;    %THIS IS HARD CODED
%OPTIONS.RepVect =  1:6; %use all
OPTIONS.RepVect =  4:6; %use later reps only
OPTIONS.MaxNperCS = 1;  %data on a given coverslip is highly correlated, so N=1 per coverslip

OPTIONS.bLinearInterpolateFloor = 1; %use linear (or higher-order) interpolation for the floor
OPTIONS.RundownCorrection = [1 1  0.92   0.84  0.82  0.79];

OPTIONS.RLOW = 2.5;
%OPTIONS.RLOW_REF_RANGE = []; %(don't use) require RESP(6) to be something for this to be HT-?

OPTIONS.RHIGH = 4.5;
%OPTIONS.RHIGH_REF_RANGE = [-inf 0.2]; %requite RESP(6) to be close to 0 for this to be HT+
%OPTIONS.RHIGH_REF_RANGE = []; %don't use
OPTIONS.RBIN = [-inf 2.25:.25:4.75 inf];

OPTIONS.nREF0 = [];  %zero-subtraction
OPTIONS.nREF1 = [1 2]; %one-normalization

MM = [-0.5 1.25];
OPTIONS.RespMinMax = {...
    [1 -2 -0.5 0.5]  [3 0.2 MM(2)]  [4 MM] [5 MM] [6 MM] ...
    {'delG' 1e3 1e4}};  %avoid excessive runup/down (factor of 2)


OPTIONS.RespMinMax = {...
    [1 -2 -0.5 0.5]  [3 0.2 MM(2)]  [4 MM] [5 MM] [6 MM] ...
    {'delG' 1e3 1e4}};  %avoid excessive runup/down (factor of 2)

[Series, CELLS] = XYW_DoseResponsePlot2(CELLS, OPTIONS, [1 4]);

figure(1); try, delete(h); catch, end;
x = 10.^[-3.5:.01:1.5]; 
m=1.3; z=x.^m; 
n=1.6; y=x.^n; 
h=plot(...
    x, 1-1.0*z./(z + 6^m), 'b', ...
    x, 1-1.0*z./(z + 5^m), 'r', ...
    x, 1-0.9*y./(y + 0.19^n), 'b', ...
    x, 1-0.9*y./(y + 0.017^n), 'r', ...
    'linewidth', 3);

%    x, 1-1.0*z./(z + 4^m), 'b', ...
%    x, 1-0.9*y./(y + 0.035^n), 'b', ...

return;

%% Code For YM90K.2 Exemplars
bCollapse = 1;
q = 4;  %weigh the largest deviations from the mean the most
for s=1:2  %series, 1=YM90K.1DART.1, 2=YM90K.1DART.2
    
    figure(20+s); clf;
    %HTpos
    TMP=vertcat(Series(s).CELLS(Series(s).Ipos).RESPONSE);  %imagesc(TMP); set(gca,'xtick',1:length(Series(s).DOSE),'xticklabel',Series(s).DOSE);
    [~,I]=sort(MeanNoNaN(abs(TMP - MeanNoNaN(TMP,1)).^q, 2))  %sort cells by similarity to the mean
    TMP = TMP(I,:);  %sorted version
    for group=1:3
        P(group).I = I(~isnan(TMP(:,2+group)));  %divvy up by doisng group 1,2,3
        P(group).I = Series(s).Ipos(P(group).I)';
    end
    %HTneg
    TMP=vertcat(Series(s).CELLS(Series(s).Ineg).RESPONSE);
    [~,I]=sort(MeanNoNaN(abs(TMP - MeanNoNaN(TMP,1)).^q, 2))  %sort cells by similarity to the mean
    TMP = TMP(I,:);  %sorted version
    for group=1:3
        G(group).I = I(~isnan(TMP(:,2+group)));  %divvy up by doisng group 1,2,3
        G(group).I = Series(s).Ineg(G(group).I)';
    end
    
    if s==1
        %%create the plot                             Positive exemplar cells          neGative exemplar cells
        PlotCSExemplars_PUBLICATION2(Series(s).CELLS, [P(1).I(5) P(2).I(1) P(3).I(1)], [G(1).I(1) G(2).I(5) G(3).I([7])], bCollapse)
    else
        %%create the plot                             Positive exemplar cells          neGative exemplar cells
        PlotCSExemplars_PUBLICATION2(Series(s).CELLS, [P(1).I(1) P(2).I(8) P(3).I(1)], [G(1).I(2) G(2).I(27) G(3).I(10)], bCollapse)
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% check a specific dose
CHECK_SERIES = length(Series);
CHECK_DOSEGRP = 1;
SORT_DOSE_ =  4;
CHECK_FLD     = 'Ipos';
CELLCHECK = Series(CHECK_SERIES).CELLS(Series(CHECK_SERIES).(CHECK_FLD));
TMP = vertcat(CELLCHECK.RESPONSE);
CELLCHECK = CELLCHECK(~isnan(TMP(:,CHECK_DOSEGRP+2)));
if SORT_DOSE_>0
    TMP = vertcat(CELLCHECK.RESPONSE_);
    [~,I] = sort(TMP(:,SORT_DOSE_), 'descend'); 
    CELLCHECK = CELLCHECK(I);
end
for n=1:length(CELLCHECK)
    figure(21); clf; imagesc(CELLCHECK(n).SHAPECOV); axis image; hold on; set(gca,'ydir','normal', 'xtick', 0:6:36);
    plot(0.5+(1:length(CELLCHECK(n).grnMin))/16, 1+5*CELLCHECK(n).grnMin, '.r', 'linewidth', 2);
    plot(0.5+(1:length(CELLCHECK(n).grnMin))/16, 1+0*CELLCHECK(n).grnMin, '.r', 'linewidth', 2);
    plot(0.5+(1:length(CELLCHECK(n).grnNorm))/16, 1+5*CELLCHECK(n).grnNorm, '.k', 'linewidth', 2);
    title([num2str(n) ' of ' num2str(length(CELLCHECK))] ); 
    colorbar; caxis([-1.1  1.1]);
    disp(CELLCHECK(n))
    pause;
end




