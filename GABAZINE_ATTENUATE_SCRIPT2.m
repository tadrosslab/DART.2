
if 0
    clear;
end

try
    CELLS;
catch
    PathList1 = { ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2018_10_31_gabazineDART_comparison\plate2_2018_11_01' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2018_10_31_gabazineDART_comparison\plate3_2018_11_02' ...
        ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2019_02_13_gabazineDART_comparison\plate01_2019_02_13' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2019_02_13_gabazineDART_comparison\plate02_2019_02_14' ...
        ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2019_02_20_gabazine_comparison_in_7Ktyrode\plate01_2019_02_20' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2019_02_20_gabazine_comparison_in_7Ktyrode\plate02_2019_2_21' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2019_02_20_gabazine_comparison_in_7Ktyrode\plate03_2019_02_21' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2019_02_20_gabazine_comparison_in_7Ktyrode\plate04_2019_02_22' ...
        };
    PathList2 = { ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_02_20_GH_ddHT_vs_HT\plate_01' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_02_20_GH_ddHT_vs_HT\plate_02' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_02_20_GH_ddHT_vs_HT\plate_03' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_02_26_GH_HT_vs_ddHT\plate01' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_02_26_GH_HT_vs_ddHT\plate02' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_02_26_GH_HT_vs_ddHT\plate03' ...
        };
    PathList3 = { ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_11_11_XYW_new_gabazine_variants\batch_01' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_11_11_XYW_new_gabazine_variants\batch_02' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_11_11_XYW_new_gabazine_variants\batch_03' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_11_18_new_gabazine_variants\batch_01' ...
        'E:\Dropbox (TadrossLab)\BIG_DATA\IX83_NEURONAL\BY_PROJECT\GABAZINE_DART2.0\2020_11_18_new_gabazine_variants\batch_02' ...
        };
    
    PathList = [PathList1  PathList2 PathList3]; 
    %PathList = [ PathList3]; 
    
    if isempty(dir(PathList{1}))
        %running on laptop -- change root folder
        for k=1:length(PathList)
            PathList{k} = strrep(PathList{k},  'E:\Dropbox (TadrossLab)\BIG_DATA',   'C:\Users\mrt31\Dropbox (TadrossLab)\BIG_DATA');
        end
    end
    CELLS = XYW_LoadCells(PathList);
    for k=1:length(CELLS)
        if strcmpi(CELLS(k).SeriesName, 'NsC00/gabazine') || strcmpi(CELLS(k).SeriesName, 'NsC07s/gabazine')
            CELLS(k).SeriesName = 'gabazine';
        elseif strcmpi(CELLS(k).SeriesName, 'NsC07s/gabazine-DART2.0') 
            CELLS(k).SeriesName = 'gabazine-DART2.0';
        end
    end
end
if 0 % 1
    %conbine three doses in the DW for better N on the activity vs red plot
    for k=1:length(CELLS)
        for j=1:length(CELLS(k).DOSE)
            if ~isempty(intersect(CELLS(k).DOSE{j}, [3000 300 1000]))
                CELLS(k).DOSE{j} = 300;
            end
        end
    end
end


%For subplots -- choose one of these
OPTIONS.REFBINS = 7;
OPTIONS.REFTXT = 'SCOREGOOD';   OPTIONS.REFLIM = [ 0  1 ];
OPTIONS.REFTXT = 'mean(grnMin)'; OPTIONS.REFLIM = [ 0  8];
OPTIONS.REFTXT = 'red1Std';   OPTIONS.REFLIM = [ 0  100 ];
OPTIONS.REFTXT = 'log10(delGstd)';   OPTIONS.REFLIM = [ 1  4 ];
OPTIONS.REFTXT = 'max(RESP)';   OPTIONS.REFLIM = [ 0.9 4 ];
OPTIONS.REFTXT = 'min(RESP)';   OPTIONS.REFLIM = [ -0.25 0.05 ];
OPTIONS.REFTXT = 'SHAPE';   OPTIONS.REFLIM = [ 0  1 ];
OPTIONS.REFTXT = 'RESPONSEcv';   OPTIONS.REFLIM = [ -1  2 ];
OPTIONS.REFTXT = 'min(RESPONSEstd)';   OPTIONS.REFLIM = [0  0.05 ];
OPTIONS.REFTXT = 'max(RESPONSEstd)';   OPTIONS.REFLIM = [0  0.2 ];
OPTIONS.REFTXT = 'mean(RESPONSEstd)';   OPTIONS.REFLIM = [0  0.07 ];
OPTIONS.REFTXT = 'REF0/REF1';   OPTIONS.REFLIM = [ 0  1 ];
OPTIONS.REFTXT = 'REF1-REF0'; OPTIONS.REFLIM = [ 0  1 ];
OPTIONS.REFTXT = 'grnMinWavyness'; OPTIONS.REFLIM = [ 0  0.5 ];
OPTIONS.REFTXT = 'FilledArea';  OPTIONS.REFLIM = [ 0 1000 ];
OPTIONS.REFTXT = 'RESP(1)';     OPTIONS.REFLIM = [ -0.08  0.16 ];
OPTIONS.REFTXT = 'RESPONSEstd';   OPTIONS.REFLIM = [ 0  0.1 ];
OPTIONS.REFTXT = 'REF0';        OPTIONS.REFLIM = [ -0.2  1.4 ];
OPTIONS.REFTXT = 'REF1';        OPTIONS.REFLIM = [ 0.25  0.65 ];
OPTIONS.REFTXT = 'RESP(1)-RESP(2)'; OPTIONS.REFLIM = [ -0.2  0.7 ];
OPTIONS.REFTXT = 'log10(delG)'; OPTIONS.REFLIM = [2 5];

OPTIONS.DoseVectPlot = [0  0.3  1 3 10  30  100  300  1000  3000  10000  30000  100000  300000]; %additional figures for given doses

if 0 %IDENTICAL TO GABAZINE_SCRIPT
    OPTIONS.YLIM = [-0.25 1.4];
    OPTIONS.bScatterLog = 1;  %for dTomato
    OPTIONS.bShowScatterPos = 1;
    OPTIONS.bShowScatterNeg = 1;
    OPTIONS.bShowLine = 0;
    
    OPTIONS.bUseMedian=0;
    OPTIONS.bUseMax = 0;  %use GCaMP max, instead of mean
    OPTIONS.MAXGrnClip = inf; %65000;
    
    OPTIONS.NumRepsPerDose = 4;    %THIS IS HARD CODED
    OPTIONS.RepVect = 1:4; %use all reps
    OPTIONS.MaxNperCS = inf;  %assume every cell is independent, even if on the same CS
    
    OPTIONS.bLinearInterpolateFloor = 2; %use (higher-order) interpolation for the floor
    OPTIONS.RundownCorrection = [1.1  1.1  1.1  1.25  1  1];
    
    OPTIONS.RLOW = 2.5;
    OPTIONS.RHIGH = 4.25;
    OPTIONS.RBIN = [-inf 2.25:.25:4.75 inf];
    
    OPTIONS.nREF0 = [1 2];  %zero-subtractions
    OPTIONS.nREF1 = [6];    %one-normalization
    
    MM=[-0.2  2];
    OPTIONS.RespMinMax = {...
        [1 -2 -0.1  0.25] [3 -0.1  2] [4 MM] [5 MM] ...
        {'delG' 1e3 1e4}};  %avoid excessive runup/down (factor of 2)
    
elseif 0  %lowest quality control, most N
    OPTIONS.YLIM = [-0.2 0.8];
    OPTIONS.bScatterLog = 1;  %for dTomato
    OPTIONS.bShowScatterPos = 1;
    OPTIONS.bShowScatterNeg = 1;
    OPTIONS.bShowLine = 0;
    
    OPTIONS.bUseMedian=0;
    OPTIONS.bUseMax = 0;  %use GCaMP max, instead of mean
    OPTIONS.MAXGrnClip = inf; %65000;
    
    OPTIONS.NumRepsPerDose = 4;    %THIS IS HARD CODED
    OPTIONS.RepVect = 1:4; %use all reps
    OPTIONS.MaxNperCS = inf;  %assume every cell is independent, even if on the same CS
    
    OPTIONS.bLinearInterpolateFloor = 2; %use (higher-order) interpolation for the floor
    OPTIONS.RundownCorrection = [];  %[1.5  1  0.95  0.9  0.85  1.1];
    
    OPTIONS.RLOW =  2.5;
    OPTIONS.RHIGH = 4.0; %4.25;
    %OPTIONS.RBIN = [-inf 2.25:.25:4.75 inf];
    OPTIONS.RBIN = [-inf 2.5:0.5:4.5 inf];
    
    OPTIONS.nREF0 = [2];  %zero-subtractions
    OPTIONS.nREF1 = [];    %one-normalization
    
    MM=[-0.2  2];
    OPTIONS.RespMinMax = {...
        [1 MM] [2 MM] [3 MM] [4 MM] [5 MM] [6 MM]...
        {'delG' 1e3 1e4}};  %avoid excessive runup/down (factor of 2)
   
elseif 1  %pseudo-optimized for this data set
    OPTIONS.YLIM = [-0.2 1.8];
    OPTIONS.bScatterLog = 1;  %for dTomato
    OPTIONS.bShowScatterPos = 1;
    OPTIONS.bShowScatterNeg = 1;
    OPTIONS.bShowLine = 0;
    
    OPTIONS.bUseMedian=0;
    OPTIONS.bUseMax = 0;  %use GCaMP max, instead of mean
    OPTIONS.MAXGrnClip = inf; %65000;
    
    OPTIONS.NumRepsPerDose = 4;    %THIS IS HARD CODED
    OPTIONS.RepVect = 1:4; %use all reps
    OPTIONS.MaxNperCS = inf;  %assume every cell is independent, even if on the same CS
    
    OPTIONS.bLinearInterpolateFloor = 2; %use (higher-order) interpolation for the floor
    OPTIONS.RundownCorrection = [1  1  0.95  0.95  0.8  0.8];
    
    OPTIONS.RLOW =  2.5;
    OPTIONS.RHIGH = 4.1; %4.25;
    %OPTIONS.RBIN = [-inf 2.25:.25:4.75 inf];
    OPTIONS.RBIN = [-inf 2.25:0.25:4.5 inf];
    
    OPTIONS.nREF0 = [];  %zero-subtractions
    OPTIONS.nREF1 = [6];    %one-normalization
    
    MM=[-0.2  2];
    OPTIONS.RespMinMax = {...
        [1 -0.2 0.4]  [1 -2 -0.1  0.25] [3 MM] [4 MM] [5 MM] ...
        {'delG' 10^2.5  inf}};  %avoid excessive runup/down (factor of 2)
end

[Series, CELLS] = XYW_DoseResponsePlot2(CELLS, OPTIONS, 5+[[1  [ 5   7 ]-1]]);

figure(1); try delete(h); catch, end
x = 10.^[-3.5:.01:3];
m=1; z=x.^m;
n=1.3; y=x.^n;
h=plot(...
    x, 1.02*z./(z + 13^m), '--k', ...
    x, 1.01*y./(y + 0.025^n), 'k', ...
    ...
    x, 1.00*z./(z + 70^m), '--b', ...
    x, 0.90*y./(y + 0.035^n), 'b', ...
    ...
    x, 1.00*z./(z + 180^m), '--r', ...
    x, 0.77*y./(y + 0.040^n), 'r', ...
    'linewidth', 3);


return;


%% Code For YM90K.2 Exemplars
bCollapse = 1;
q = 4;  %weigh the largest deviations from the mean the most
for s=3  %series, 1=YM90K.1DART.1, 2=YM90K.1DART.2
    
    figure(20+s); clf;
    %HTpos
    TMP=vertcat(Series(s).CELLS(Series(s).Ipos).RESPONSE);  %imagesc(TMP); set(gca,'xtick',1:length(Series(s).DOSE),'xticklabel',Series(s).DOSE);
    [~,I]=sort(MeanNoNaN(abs(TMP - MeanNoNaN(TMP,1)).^q, 2))  %sort cells by similarity to the mean
    TMP = TMP(I,:);  %sorted version
    for group=1:3
        P(group).I = I(~isnan(TMP(:,1+group)));  %divvy up by doisng group 1,2,3
        P(group).I = Series(s).Ipos(P(group).I)';
    end
    %HTneg
    TMP=vertcat(Series(s).CELLS(Series(s).Ineg).RESPONSE);
    [~,I]=sort(MeanNoNaN(abs(TMP - MeanNoNaN(TMP,1)).^q, 2))  %sort cells by similarity to the mean
    TMP = TMP(I,:);  %sorted version
    for group=1:3
        G(group).I = I(~isnan(TMP(:,1+group)));  %divvy up by doisng group 1,2,3
        G(group).I = Series(s).Ineg(G(group).I)';
    end
    
    if s==3
        %%create the plot                             Positive exemplar cells          neGative exemplar cells
        PlotCSExemplars_PUBLICATION2(Series(s).CELLS, [P(1).I(1) P(2).I(5) P(3).I(7)], [G(1).I(1) G(2).I(2) G(3).I(3)], bCollapse, 4)
        axis([-10 360 -30 0]);
        %axis off;
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check a specific dose
CHECK_SERIES = length(Series);
CHECK_DOSEGRP = 1;
SORT_DOSE_ =  4;
CHECK_FLD     = 'Ipos';
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
    title([num2str(n) ' of ' num2str(length(CELLCHECK))] );
    colorbar; caxis([-1.1  1.1]);
    disp(CELLCHECK(n))
    pause;
end



