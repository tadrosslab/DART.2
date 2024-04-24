

function [Series, CELLS_ALL, CELLS] = XYW_DoseResponsePlot2(CELLS, OPTIONS, SERIESORDER, SERIESCLR)
%function [Series, CELLS_ALL, CELLS] = XYW_DoseResponsePlot2(CELLS, OPTIONS, SERIESORDER, SERIESCLR)

set(0,'DefaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Default input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    %undo all the proofreading
    CELLS([CELLS.bViewed]==2) = [];     %eliminate hand-drawn cells
    for k=1:length(CELLS)
        CELLS(k).bGood = 1;     %mark all good
    end
end

try
    USEROPT = OPTIONS;
    clear OPTIONS;
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defulat OPTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
OPTIONS.REFTXT = 'RESP(1)-RESP(2)'; OPTIONS.REFLIM = [ -0.2  0.4 ];
OPTIONS.REFTXT = 'log10(delG)'; OPTIONS.REFLIM = [2 5];
OPTIONS.REFTXT = 'RESP(1)';     OPTIONS.REFLIM = [ -0.08  0.16 ];
OPTIONS.REFTXT = 'RESPONSEstd';   OPTIONS.REFLIM = [ 0  0.1 ];
OPTIONS.REFTXT = 'REF0';        OPTIONS.REFLIM = [ -0.2  1.4 ];
OPTIONS.REFTXT = 'REF1';        OPTIONS.REFLIM = [ 0.25  0.65 ];

OPTIONS.DoseVectPlot = [ 10  30  1000  3000  10000  30000];  %additional figures for given doses

OPTIONS.YLIM = [-0.25 1.5];
OPTIONS.bScatterLog = 1;  %for dTomato
OPTIONS.bShowScatterPos = 1;
OPTIONS.bShowScatterNeg = 1;
OPTIONS.bShowLine = 1;
OPTIONS.LegendLocation = 'northwest';

OPTIONS.bUseMedian=0;
OPTIONS.bUseMax = 0;  %use GCaMP max, instead of mean
OPTIONS.MAXGrnClip = inf; %65000;

OPTIONS.NumRepsPerDose = 6;    %THIS IS HARD CODED
OPTIONS.RepVect = 1:OPTIONS.NumRepsPerDose; %use all reps
OPTIONS.MaxNperCS = inf;  %assume every cell is independent, even if on the same CS
OPTIONS.EF_Flag = 0;   %special case for atropine protocol

OPTIONS.RLOW = 2.5;
OPTIONS.RLOW_REF_RANGE = []; %don't use it
OPTIONS.RHIGH = 4.0;
OPTIONS.RHIGH_REF_RANGE = []; %don't use it
OPTIONS.RBIN = 30;

OPTIONS.nREF0 = []; %don't do it
OPTIONS.nREF1 = []; %don't do it
OPTIONS.dwREF0 = 1;  %denominator weigting of REF0; used to address strong correlations between data and REF0/REF1
OPTIONS.ScalingPostREF = 1; %used in concert with dwREF0

OPTIONS.RespMinMax = {};  %Exclude cells based on certain criteria

OPTIONS.bLinearInterpolateFloor = 0; %1= use linear interpoloation; 2= use higher-order interpolation for the floor
OPTIONS.RundownCorrection = [];
OPTIONS.bFixRESP16 = 0;

%overrwite with user-specified options
try
    fld = fieldnames(USEROPT);
    for k=1:length(fld)
        OPTIONS.(fld{k}) = USEROPT.(fld{k});
    end
catch
end

if length(OPTIONS.RLOW)==1
    OPTIONS.RLOW = [-inf OPTIONS.RLOW]; %min/max
end
if length(OPTIONS.RHIGH)==1
    OPTIONS.RHIGH = [OPTIONS.RHIGH inf]; %min/max
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADDITIONAL PROCESSING OF EACH CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(CELLS,'SHAPE') || OPTIONS.EF_Flag  %save time if this is a re-call
    for k=1:length(CELLS)
        %                                          length(grnNorm) typically 480 = 16t x 6rep x 5doses.                                 %                                          typically=
        NUMREP = OPTIONS.NumRepsPerDose*length(CELLS(k).DOSE);    %tyipcally  30 =       6rep x 5doses.  total number of repetitions of the whole assay
        NperDose = length(CELLS(k).grnNorm)/length(CELLS(k).DOSE); %typically  96 = 16t x 6rep            number of grn timepoints for a given dose
        NperRep  = NperDose/OPTIONS.NumRepsPerDose;               %typically       16t                   number of grn timepoints for a single repetition
        
        %Add new field
        try
            if isempty(CELLS(k).red1Pix)
                error('no red1');
            end
            CELLS(k).redMed  = 4*median(double(CELLS(k).red1Pix)); %wider dynamic range (less saturation)            
        catch
            %red1Pix does not exist or is empty
            CELLS(k).redMed  = median(double(CELLS(k).redPix)); %prior to new dual-TRITC
        end
        try
            CELLS(k).blueMed  = median(double(CELLS(k).bluePix)); 
        catch
        end
        
        if OPTIONS.EF_Flag
            %grnMin = min(CELLS(k).grnRaw(OPTIONS.EF_nPRE));  %only consdier nPRE for grnMin (may allow some negative data points)
            grnMin = min(CELLS(k).grnRaw);    %min over entire waveform
            grnMax = max(CELLS(k).grnRaw);    %max over entire waveform
            grnNorm = (CELLS(k).grnRaw - grnMin)/(grnMax - grnMin); %normalize
                                    
            CELLS(k).meanPRE = mean(grnNorm(OPTIONS.EF_nPRE));
            CELLS(k).meanSIG = mean(grnNorm(OPTIONS.EF_nSIG));
            CELLS(k).meanPST = mean(grnNorm(OPTIONS.EF_nPST));
            CELLS(k).stdPRE  =  std(grnNorm(OPTIONS.EF_nPRE));
            CELLS(k).stdSIG  =  std(grnNorm(OPTIONS.EF_nSIG));
            CELLS(k).stdPST  =  std(grnNorm(OPTIONS.EF_nPST));
            
            PRE = CELLS(k).meanPRE;
            SIG = CELLS(k).meanSIG;
            PST = CELLS(k).meanPST;
            
            stdPRE = CELLS(k).stdPRE;
            stdSIG = CELLS(k).stdSIG;
            stdPST = CELLS(k).stdPST;
            
            RED = CELLS(k).redMed;
            SIG = eval(OPTIONS.EF_SIG_Nonlinear); %'1.8 + 0.7*log10(SIG)
            if ~isreal(SIG)
                SIG = NaN;
            end
        else
            %Perform a more robust determination of the floor
            %First determine and subtract min
            grnNorm = CELLS(k).grnRaw;
            if OPTIONS.bLinearInterpolateFloor==2
                %piecewise fitting for each dose.
                %connect the dots through starting value of each repetition
                grnMin = grnNorm*0; %starting
                for q=1:length(CELLS(k).DOSE)
                    grnTmp = grnNorm( (1+(q-1)*NperDose):(q*NperDose) );  %green data for this dose
                    nFit   = 1:NperRep:length(grnTmp);
                    grnFit = grnTmp(nFit);
                    P = polyfit(nFit, grnFit, 2); %2nd order polynomial fit
                    grnMin( (1+(q-1)*NperDose):(q*NperDose) ) = polyval(P, 1:length(grnTmp));  %piecewise
                end
                %nNeg = find(grnMin > grnNorm);
                %if ~isempty(nNeg)
                %    grnMin(nNeg) = grnNorm(nNeg); %lowest
                %    %figure(21); clf; plot(grnNorm,'.'); hold on; plot(grnMin, '.');
                %    %hold on;
                %end
            elseif OPTIONS.bLinearInterpolateFloor==1
                %a linear fit through the minimum of each dose
                X = 1:length(CELLS(k).DOSE);
                Y = X*0; %allocate memory & erase prior use of this variable
                for q=1:length(X)
                    if strcmpi(CELLS(k).DOSE{q},'NaN')
                        Y(q) = NaN;  %missing dose
                    else
                        GrnSort = sort(grnNorm( (1+(q-1)*NperDose):(q*NperDose) ));
                        Y(q) = mean(GrnSort(1:ceil(0.05*length(GrnSort))));    %bottom 5%
                    end
                end
                X = X(~isnan(Y));
                Y = Y(~isnan(Y));
                try
                    P = polyfit(X,Y,1); %linear fit
                    %for q=1:2 %eliminate the two midpoints that are the highest above the linear regression
                    %    [~,I] = max(Y(2:end-1)-polyval(P,X(2:end-1)));
                    %    X(I+1) = [];
                    %    Y(I+1) = [];
                    %    P = polyfit(X,Y,1); %linear fit
                    %end
                    grnMin = polyval(P, (1:length(grnNorm))/NperDose-0.5 );
                end
            else
                %get the global minimum
                GrnSort = sort(grnNorm);
                grnMin = grnNorm*0 + mean(GrnSort(1:ceil(0.05*length(GrnSort))));       %bottom 5%
            end
            
            %subtract min (same code whether it's a constant or a linear interpolation)
            grnNorm = grnNorm - grnMin;
            
            %determine and normalize max
            GrnSort = sort(grnNorm);
            grnMax = mean(GrnSort(end-ceil(0.01*length(GrnSort)):end)); %top 1%
            grnNorm = grnNorm/grnMax;
            grnMin  = grnMin/grnMax;
            
            %DETERMINE Min/Max/Mean of each repeition, as well as Response variability
            SIG = zeros(NperRep, NUMREP); SIG(:)=grnNorm;  %convert to 1 column per repetition
        end
        repMin = min(SIG);  %floor response for each rep
        repMax = max(SIG);  %ceiling response for each rep
        repMean = mean(SIG);  %mean response for each rep
        
        %CALCULATE SHAPE COVARIANCE
        SHAPECOV = zeros(NUMREP, NUMREP);  %measure of shape similarity of each rep to every other rep
        SIGz = SIG - repmat(repMean, NperRep, 1); %normalize
        for r=1:NUMREP
            for q=r:NUMREP
                SHAPECOV(r,q) =  sum(SIGz(:,r).*SIGz(:,q))/(NperRep-1);
                SHAPECOV(q,r) =  SHAPECOV(r,q);  %symmetry
            end
        end
        %normalize to power
        for r=1:NUMREP
            PWR(r) = SHAPECOV(r,r);
        end
        for r=1:NUMREP
            for q=1:NUMREP
                SHAPECOV(r,q) =  SHAPECOV(r,q)/sqrt(PWR(r)*PWR(q));
            end
        end
        bDebug = 0;
        if bDebug
            clf; imagesc(SHAPECOV); axis image; hold on; set(gca,'ydir','normal'); plot(0.5+(1:length(grnNorm))/NperRep, 1+5*grnNorm, 'k', 'linewidth', 2); colorbar; caxis([-1.1  1.1]);
        end
        
        %save this info to the Blob
        CELLS(k).grnMax  = grnMax;
        CELLS(k).grnMin  = grnMin;
        CELLS(k).grnMaxOverMin  = grnMax./grnMin;
        CELLS(k).grnNorm = grnNorm;
        CELLS(k).repMin = repMin;
        CELLS(k).repMax = repMax;
        CELLS(k).repMean = repMean;
        
        CELLS(k).SHAPECOV = SHAPECOV;
        CELLS(k).SHAPE = mean(mean(CELLS(k).SHAPECOV(:,:)));
        
        
        %convert dose to numeric
        TMP = CELLS(k).DOSE;
        if iscell(TMP)
            CELLS(k).DOSE = [];
            for j=1:length(TMP)
                if isnumeric(TMP{j})
                    CELLS(k).DOSE(j) = TMP{j};
                elseif strcmpi(TMP{j},'P')
                    CELLS(k).DOSE(j) = -1;
                elseif strcmpi(TMP{j},'B')
                    CELLS(k).DOSE(j) = 0;
                elseif strcmpi(TMP{j},'W')
                    CELLS(k).DOSE(j) = inf;
                elseif strcmpi(TMP{j},'NaN')
                    CELLS(k).DOSE(j) = -inf; %special code to ignore
                else
                    error('Unrecognized Dose');
                end
            end
        end
    end
end
CELLS_ALL = CELLS; %output parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1)RepVect-->RESPONSE;  (2)RundownCorrection;  (3)REF0/REF1 normalization;  (4)Min/Max criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(CELLS)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%STEP#1 =  RepVect calculation of RESPONSE, RESPONSEstd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RESPONSE      = zeros(1,length(CELLS(k).DOSE));
    RESPONSEcv    = zeros(1,length(CELLS(k).DOSE));
    RESPONSEstd   = zeros(1,length(CELLS(k).DOSE));
    RESPONSEfloor = zeros(1,length(CELLS(k).DOSE));
    ROBUSTfloor   = zeros(1,length(CELLS(k).DOSE));
    r = 0; %total # reps
    for d=1:length(CELLS(k).DOSE)
        %Variability in activity for a single dose (note there are 6 reps per dose):
        %Average activity for a single dose
        if OPTIONS.bUseMax
            TMP = CELLS(k).repMax((d-1)*OPTIONS.NumRepsPerDose + OPTIONS.RepVect);
        else
            TMP = CELLS(k).repMean((d-1)*OPTIONS.NumRepsPerDose + OPTIONS.RepVect);
        end
        RESPONSE(d) = mean(TMP); %activity score is the mean scaled GCaMP signal for this dose                end
        RESPONSEstd(d) = std(TMP); %coefficient of variation; measure of reproducibility on repetitions.
        RESPONSEcv(d) = RESPONSEstd(d)/abs(RESPONSE(d)); %coefficient of variation; measure of reproducibility on repetitions.
        
        RESPONSEfloor(d) = min(CELLS(k).repMin((d-1)*OPTIONS.NumRepsPerDose + OPTIONS.RepVect));
        ROBUSTfloor(d)   = mean(CELLS(k).repMin((d-1)*OPTIONS.NumRepsPerDose + OPTIONS.RepVect));
    end
    FloorDelta = 0;
    if any(ROBUSTfloor>0)
        FloorDelta = FloorDelta + max(ROBUSTfloor);
    end
    if any(ROBUSTfloor<0)
        FloorDelta = FloorDelta - min(ROBUSTfloor);
    end
    
    CELLS(k).RESPONSE = RESPONSE;
    CELLS(k).RESPONSEcv = RESPONSEcv;
    CELLS(k).RESPONSEstd = RESPONSEstd;
    CELLS(k).RESPONSEfloor = RESPONSEfloor;
    CELLS(k).ROBUSTfloor = ROBUSTfloor;
    
    CELLS(k).FloorDelta = FloorDelta;
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%STEP#2 =  RundownCorrection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %NOTE: Rundown correction performed BEFORE RespMinMax
    try
        if ~isempty(OPTIONS.RundownCorrection)
            CELLS(k).RESPONSE = CELLS(k).RESPONSE./OPTIONS.RundownCorrection;
        end
    catch
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%STEP#3 =  R0 & R1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %normalize to the reference response
    if ~isempty(OPTIONS.nREF0)
        REF0          = mean(CELLS(k).RESPONSE(OPTIONS.nREF0));%this code allows for more than one reference
    else
        REF0 = 0;
    end
    if ~isempty(OPTIONS.nREF1)
        REF1          = mean(CELLS(k).RESPONSE(OPTIONS.nREF1));%this code allows for more than one reference
    else
        REF1 = 1;
    end
    %subtract then divide
    CELLS(k).RESPONSE     = ((CELLS(k).RESPONSE - REF0)./(REF1 - OPTIONS.dwREF0*REF0))./OPTIONS.ScalingPostREF;
    CELLS(k).RESPONSEstd  = (      CELLS(k).RESPONSEstd./(REF1 - OPTIONS.dwREF0*REF0))./OPTIONS.ScalingPostREF;
    CELLS(k).RESPONSEcv   =  CELLS(k).RESPONSEstd./abs(CELLS(k).RESPONSE);
    
    CELLS(k).REF0 = REF0;
    CELLS(k).REF1 = REF1;

    try
        if OPTIONS.bFixRESP16
            CELLS(k).RESPONSE(2) = CELLS(k).RESPONSE(2)  - 0.3*(CELLS(k).RESPONSE(1)-CELLS(k).RESPONSE(6));
            CELLS(k).RESPONSE(3) = CELLS(k).RESPONSE(3)  - 0.1*(CELLS(k).RESPONSE(1)-CELLS(k).RESPONSE(6));
            CELLS(k).RESPONSE(4) = CELLS(k).RESPONSE(4)  + 0.1*(CELLS(k).RESPONSE(1)-CELLS(k).RESPONSE(6));
        end
    catch
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%STEP#4 =  RespMinMax exclusion criteria.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if CELLS(k).grnMax > OPTIONS.MAXGrnClip
        %eliminate any cells that saturate GCaMP signal
        CELLS(k).bGood = 0;
    end
    
    for q=1:length(OPTIONS.RespMinMax)
        MinMax = OPTIONS.RespMinMax{q};
        if iscell(MinMax)
            %special code
            field = MinMax{1};
            MinMax = [MinMax{2:end}];
            switch(field)
                case 'REF0/REF1'
                    VAL = CELLS(k).REF0 / CELLS(k).REF1;
                case 'REF1-REF0'
                    VAL = CELLS(k).REF1 - CELLS(k).REF0;
                case 'mean(grnMin)'
                    VAL = mean(CELLS(k).grnMin);
                case 'mean(grnRaw)'
                    VAL = mean(CELLS(k).grnRaw);
                case 'mean(grnNorm)'
                    VAL = mean(CELLS(k).grnRaw);
                case 'grnMinWavyness'
                    VAL = CELLS(k).grnMin;  %Ncell x t
                    VAL = VAL - mean(VAL); %zero-offset
                    VAL = mean(abs(VAL));
                case 'mean(RESPONSEstd)'
                    VAL = CELLS(k).RESPONSEstd;
                    VAL = MeanNoNaN(VAL,2);
                case 'max(RESPONSEstd)'
                    VAL = CELLS(k).RESPONSEstd;
                    VAL = max(VAL,[],2);
                case 'min(RESPONSEstd)'
                    VAL = CELLS(k).RESPONSEstd;
                    VAL = min(VAL,[],2);
                otherwise
                    VAL = CELLS(k).(field);
            end
            
        elseif length(MinMax)==3
            if CELLS(k).DOSE(MinMax(1)) == -Inf
                continue; %skip this condition; there is a bad dose.
            end
            VAL = CELLS(k).RESPONSE(MinMax(1));
            
            
        elseif length(MinMax)==4
            if any(CELLS(k).DOSE(abs(MinMax(1:2))) == -Inf)
                continue; %skip this condition; there is a bad dose.
            end
            if MinMax(2)<0
                %[DOSE1-DOSE2  MIN MAX]
                VAL = CELLS(k).RESPONSE(MinMax(1)) - CELLS(k).RESPONSE(-MinMax(2));
            else
                %[DOSE1/DOSE2  MIN MAX]
                VAL = CELLS(k).RESPONSE(MinMax(1)) / CELLS(k).RESPONSE(MinMax(2));
            end
        else
            error('OPTIONS.RespMinMax must be a cell array containing 3- or 4-element vectors');
        end
        
        if length(VAL)==1
            if VAL < MinMax(end-1)
                CELLS(k).bGood = 0;
                continue;
            end
            if VAL > MinMax(end)
                CELLS(k).bGood = 0;
                continue;
            end
        elseif length(VAL) == length(CELLS(k).RESPONSE)
            Ibad = ( VAL < MinMax(end-1) |  VAL > MinMax(end));
            CELLS(k).RESPONSE(Ibad) = NaN; %block these
        else
            error('VAL min/max criteria has an issue');
        end
    end
end
CELLS = CELLS([CELLS.bGood]>0); %remove all bad cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GROUPING OF CELLS INTO SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find unique series
[USeriesName,~,IUS] = unique({CELLS.SeriesName});
for k=1:length(USeriesName)
    Series(k).Name = USeriesName{k};
    Series(k).CELLS = CELLS(IUS==k); %cells for this SeriesName
    Series(k).DOSE  = unique([Series(k).CELLS.DOSE]);
    Series(k).DOSE(Series(k).DOSE==-Inf) = [];
    %convert each cell to a full DOSE/RESPONSE
    
    for j=1:length(Series(k).CELLS)
        %indices of doses present in this cell, mapped to all the doses over all cells of this series
        [~,B,C] = intersect(Series(k).CELLS(j).DOSE, Series(k).DOSE);
        FldNames    = {'RESPONSEfloor' 'ROBUSTfloor' 'RESPONSEcv' 'RESPONSEstd' 'RESPONSE'};
        for f=1:length(FldNames)
            RESP = Series(k).CELLS(j).(FldNames{f});     %original responses
            Series(k).CELLS(j).(FldNames{f}) = NaN*Series(k).DOSE; %more doses, all initially set to NaN
            Series(k).CELLS(j).(FldNames{f})(C) = RESP(B);   %responses mapped to new doeses
            Series(k).CELLS(j).([FldNames{f} '_']) = RESP;   %original responses
        end
        Series(k).CELLS(j).DOSE = Series(k).DOSE;  %overwrite
        Series(k).CELLS(j).PathCS = [Series(k).CELLS(j).Path '_CS' num2str(Series(k).CELLS(j).CS,'%.2d')];
    end
end
try
    if ~isempty(SERIESORDER)
        Series = Series(SERIESORDER);
    end
catch
end

if nargin==1
    return;  %no need to plot for pre-run
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot
FS = 12; %font size

try
    clrs = SERIESCLR;
catch
    if length(Series)==2
        clrs={[0 0 0]  [1 0 0]};
    elseif length(Series)==3
        clrs={[0.5 0.5 0.5]  [0 0 1]  [1 0 0]};
    elseif length(Series)==4
        clrs={[0 0 0]  [0 0 1]  [0 0.8 1]  [1 0 0]};
    elseif length(Series)==5
        clrs={[0 0 0]  [0 0 1]  [0 0.8 1] [0 0.8 0] [1 0 0]};
    elseif length(Series)==6
        clrs={[0 0 0]  [0 0 1]  [0 0.8 1] [0 0.8 0] [1 0.8 0] [1 0 0]};
    else
        for k=1:length(Series)
            clrs{k} = MapValue2Color(k,[1 length(Series)],jet);
        end
    end
end
for k=1:length(Series)
    dclrs{k} = clrs{k}*0.8; %dark colors
end
%initialize axes
for f=1:1+length(OPTIONS.DoseVectPlot)
    figure(f); clf;
    if f>1
        subplot(2,1,1);
    end
    hold on;
    for k=1:length(Series)
        plot([0 0],[0 1],'o-','color',clrs{k},'linewidth', 2);
        LegTxt{k} = strrep(Series(k).Name,'_',' ');
    end
    if 1 %f==1
        legend(LegTxt,'AutoUpdate','off', 'location', OPTIONS.LegendLocation);
    end
end
%drawnow;


if length(OPTIONS.REFTXT)>6 && strcmpi(OPTIONS.REFTXT(1:6),'log10(') && strcmpi(OPTIONS.REFTXT(end),')')
    bLog10 = 1;
    OPTIONS.REFTXT = OPTIONS.REFTXT(7:end-1);
else
    bLog10 = 0;
end

for k=length(Series):-1:1
    
    REF0  = vertcat(Series(k).CELLS.REF0);
    REF1  = vertcat(Series(k).CELLS.REF1);
    RED  = vertcat(Series(k).CELLS.redMed);
    RESP = vertcat(Series(k).CELLS.RESPONSE);
    RESP_ = vertcat(Series(k).CELLS.RESPONSE_);  %original; no NaN
    GMX  = vertcat(Series(k).CELLS.grnMax);
    
    switch OPTIONS.REFTXT
        case 'RESP(1)'
            REF = RESP_(:,1);
        case 'RESP(2)'
            REF = RESP_(:,2);
        case 'RESP(3)'
            REF = RESP_(:,3);
        case 'RESP(4)'
            REF = RESP_(:,4);
        case 'RESP(5)'
            REF = RESP_(:,5);
        case 'RESP(6)'
            REF = RESP_(:,6);
        case 'RESP(7)'
            REF = RESP_(:,7);
        case 'RESP(1)-RESP(2)'
            REF = RESP_(:,1) - RESP_(:,2);
        case 'RESP(1)-RESP(6)'
            REF = RESP_(:,1) - RESP_(:,6);
        case 'RESP(6)-RESP(1)'
            REF = RESP_(:,6) - RESP_(:,1);
        case 'RESP(1)-RESP(3)'
            REF = RESP_(:,1) - RESP_(:,3);
        case 'RESP(2)-RESP(3)'
            REF = RESP_(:,2) - RESP_(:,3);
        case  'REF0/REF1'
            REF = REF0./REF1;
        case 'REF1-REF0'
            REF = REF1 - REF0;
        case 'mean(grnMin)'
            TMP = vertcat(Series(k).CELLS.grnMin);  %Ncell x t
            REF = mean(TMP,2);
            REF = REF - (RED./GMX)/5;  %spectral crosstalk
        case 'mean(grnRaw)'
            TMP = vertcat(Series(k).CELLS.grnRaw);  %Ncell x t
            REF = mean(TMP,2);
            REF = REF - (RED)/5;  %spectral crosstalk
        case 'mean(grnNorm)'
            TMP = vertcat(Series(k).CELLS.grnNorm);  %Ncell x t
            REF = mean(TMP,2);
            %REF = REF - (RED./GMX)/5;  %spectral crosstalk
        case 'grnMinWavyness'
            TMP = vertcat(Series(k).CELLS.grnMin);  %Ncell x t
            TMP = TMP - repmat(mean(TMP,2),1,size(TMP,2)); %zero-offset
            REF = mean(abs(TMP),2);
        case 'min(RESP)'
            REF = min(RESP,[],2);
        case 'max(RESP)'
            REF = max(RESP,[],2);
        case 'mean(RESPONSEstd)'
            REF = vertcat(Series(k).CELLS.RESPONSEstd);
            REF = MeanNoNaN(REF,2);
        case 'min(RESPONSEstd)'
            REF = vertcat(Series(k).CELLS.RESPONSEstd);
            REF = min(REF,[],2);
        case 'max(RESPONSEstd)'
            REF = vertcat(Series(k).CELLS.RESPONSEstd);
            REF = max(REF,[],2);
        otherwise
            REF = vertcat(Series(k).CELLS.(OPTIONS.REFTXT));
    end
    if bLog10
        REF = log10(REF);
    end
    
    
    [UniquePathCS,~,IndxUPCS] = unique({Series(k).CELLS.PathCS});
    
    figure(1);
    X = Series(k).DOSE;  %uM dose for plotting
    tmp = X(X>0 & X<inf);
    X(X==0) = min(tmp)/10;
    X(X==-1) = min(tmp)/100;
    X(isinf(X)) = max(tmp)*10;
    X = X/1000; %convert to uM
    
    I = RED>=10^(OPTIONS.RLOW(1)) &  RED<=10^(OPTIONS.RLOW(2)); %HTneg
    if ~isempty(OPTIONS.RLOW_REF_RANGE)
        I = I & REF>=OPTIONS.RLOW_REF_RANGE(1) &  REF<=OPTIONS.RLOW_REF_RANGE(2);
    end
    Series(k).Ineg = find(I);
    if ~isempty(Series(k).Ineg)
        Y = RESP(I, :);
        C = IndxUPCS(I);
        N    =  Y*0;
        nCS  =  Y*0;
        for j=1:length(UniquePathCS)
            I = find(C==j); %cells from the same CS
            %ensure that there are no more than MaxNperCS for this CS
            if length(I)<= OPTIONS.MaxNperCS
                nPerCell = 1;  %a cell can count for at most N=1
            else
                nPerCell = OPTIONS.MaxNperCS/length(I);
            end
            N(I,:) = N(I,:) + nPerCell;
            nCS(I,:) = nCS(I,:) + 1/length(I);  %count # of coverslips
        end
        [Ymean, Ysem, ~, Ymid] = WeightedMeanSEMNoNaN(Y, N, 1);
        if ~OPTIONS.bUseMedian
            Ymid=Ymean;
        end
        if OPTIONS.bShowScatterNeg
            Xj = repmat(X,size(Y,1),1);
            for j=1:size(Xj,1)
                Xj(j,:) = Xj(j,:)*10^(0.05*randn);
            end
            plot(Xj,Y,'o','color', dclrs{k}*0.4 + [1 1 1]*0.6);
        end
        if OPTIONS.bShowLine
            sym = 'o-';
            dsym = 'o--';
        else
            sym = 'o';
            dsym = 'o';
        end
        semilogx(X, Ymid, dsym, 'color', dclrs{k}, 'linewidth', 2);
        ErrorPlot(gca, X, 0,0, Ymid, Ysem, Ysem, -0.5, 0, 2, dclrs{k});
        
        %Code added 1/6/24 for NatMeth final
        Series(k).DRneg.X    = X;
        Series(k).DRneg.Ymid = Ymid;
        Series(k).DRneg.Ysem = Ysem;
        Series(k).DRneg.nCS  = sum(nCS,'omitnan');  %total number of cverslips
        Series(k).DRneg.n    = sum(N>0,'omitnan');  %total number of cells
        Series(k).DRneg.N    = sum(N,'omitnan');    %N as defined for this analysis
    end
    
    I = RED>=10^(OPTIONS.RHIGH(1)) &  RED<=10^(OPTIONS.RHIGH(2)); %HTpos
    if ~isempty(OPTIONS.RHIGH_REF_RANGE)
        I = I & REF>=OPTIONS.RHIGH_REF_RANGE(1) &  REF<=OPTIONS.RHIGH_REF_RANGE(2);
    end
    Series(k).Ipos = find(I);
    if ~isempty(Series(k).Ipos)
        Y = RESP(I, :);
        C = IndxUPCS(I);
        N    =  Y*0;
        nCS  =  Y*0;
        for j=1:length(UniquePathCS)
            I = find(C==j); %cells from the same CS
            %N(I,:) = N(I,:) + 1; %1/length(I);
            %ensure that there are no more than MaxNperCS for this CS
            if length(I)<= OPTIONS.MaxNperCS
                nPerCell = 1;  %a cell can count for at most N=1
            else
                nPerCell = OPTIONS.MaxNperCS/length(I);
            end
            N(I,:) = N(I,:) + nPerCell;
            nCS(I,:) = nCS(I,:) + 1/length(I);  %count # of coverslips
        end
        [Ymean, Ysem, ~, Ymid] = WeightedMeanSEMNoNaN(Y, N, 1);
        if ~OPTIONS.bUseMedian
            Ymid=Ymean;
        end
        if OPTIONS.bShowScatterPos
            Xj = repmat(X,size(Y,1),1);
            for j=1:size(Xj,1)
                Xj(j,:) = Xj(j,:)*10^(0.05*randn);
            end
            plot(Xj,Y,'o','color', clrs{k}*0.4 + [1 1 1]*0.6);
        end
        if OPTIONS.bShowLine
            sym = 'o-';
            dsym = 'o--';
        else
            sym = 'o';
            dsym = 'o';
        end
        semilogx(X, Ymid, sym, 'color', clrs{k}, 'linewidth', 2);
        ErrorPlot(gca, X, 0,0, Ymid, Ysem, Ysem, -0.5, 0, 2, clrs{k});
        
        %Code added 1/6/24 for NatMeth final
        Series(k).DRpos.X    = X;
        Series(k).DRpos.Ymid = Ymid;
        Series(k).DRpos.Ysem = Ysem;
        Series(k).DRpos.nCS  = sum(nCS,'omitnan');  %total number of cverslips
        Series(k).DRpos.n    = sum(N>0,'omitnan');  %total number of cells
        Series(k).DRpos.N    = sum(N,'omitnan');    %N as defined for this analysis
    end
    
    set(gca,'xscale','log');
    set(gca,'xlim', [min(X)/2 max(X)*2]);
    set(gca,'ylim', OPTIONS.YLIM);
    grid on;
    for j=1:length(X)
        if Series(k).DOSE(j) <= 0
            XL{j} = 'Baseline';
        elseif isinf(Series(k).DOSE(j))
            XL{j} = 'Wash';
        else
            XL{j} = num2str(X(j));
        end
    end
    set(gca,'xtick',X, 'xticklabel', XL);
    set(gca, 'fontsize', FS);
    xlabel('DART dose (uM)', 'fontweight', 'bold', 'fontsize', FS);
    ylabel('Assay Activity', 'fontweight', 'bold', 'fontsize', FS);
    
    
    for q=1:length(OPTIONS.DoseVectPlot)
        figure(q+1);
        
        subplot(2,1,1);
        if OPTIONS.bScatterLog
            axis([2 5.2 -0.1 1.3]);
            plot([OPTIONS.RLOW(1)  OPTIONS.RLOW(1)], [-10  10],'color',[0 0 0]+0.5, 'linewidth', 2);
            plot([OPTIONS.RLOW(2)  OPTIONS.RLOW(2)], [-10  10],'color',[0 0 0]+0.5, 'linewidth', 2);
            plot([OPTIONS.RHIGH(1) OPTIONS.RHIGH(1)],[-10  10],'color',[0 0 0]+0.5, 'linewidth', 2);
            plot([OPTIONS.RHIGH(2) OPTIONS.RHIGH(2)],[-10  10],'color',[0 0 0]+0.5, 'linewidth', 2);
        else
            axis([0 10^5.2 -0.1 1.3]);
            plot(10.^[OPTIONS.RLOW(1)  OPTIONS.RLOW(1)], [-10 10],'color',[0 0 0]+0.5, 'linewidth', 2);
            plot(10.^[OPTIONS.RLOW(2)  OPTIONS.RLOW(2)], [-10 10],'color',[0 0 0]+0.5, 'linewidth', 2);
            plot(10.^[OPTIONS.RHIGH(1) OPTIONS.RHIGH(1)],[-10 10],'color',[0 0 0]+0.5, 'linewidth', 2);
            plot(10.^[OPTIONS.RHIGH(2) OPTIONS.RHIGH(2)],[-10 10],'color',[0 0 0]+0.5, 'linewidth', 2);
        end
        %get the relevant data to plot
        Y = RESP(:, find(Series(k).DOSE == OPTIONS.DoseVectPlot(q)) ); %column vector; single dose
        if OPTIONS.bScatterLog
            R = log10(RED);
        else
            R = RED;
        end
        C = IndxUPCS; %unique coverslip info
        %eliminate NaNs at this stage
        C = C(~isnan(Y));
        R = R(~isnan(Y));
        Y = Y(~isnan(Y)); %must do this last
        %effective N calculation
        N    =  Y*0 + 1;
        for j=1:length(UniquePathCS)
            if length(OPTIONS.RBIN)==1
                %a coarse approximation; blind to coverslip
                I = find(C==j); %cells from the same CS
                if length(I)> OPTIONS.MaxNperCS
                    N(I) = OPTIONS.MaxNperCS/length(I);
                end
           else
                %specific OPTIONS.RBIN boundaries
                for w=1:length(OPTIONS.RBIN)-1
                    I = find(C==j & R>OPTIONS.RBIN(w) & R<=OPTIONS.RBIN(w+1)); %cells in the same CS and same OPTIONS.RBIN
                    if length(I)> OPTIONS.MaxNperCS
                        N(I) = OPTIONS.MaxNperCS/length(I);
                    end
                end
            end
        end
        
        if OPTIONS.bShowLine
            sym = 'o-';
            dsym = 'o--';
        else
            sym = 'o';
            dsym = 'o';
        end
        PlotBinnedData(OPTIONS.RBIN, R, Y, N, [1 0], clrs{k}, sym, 2, OPTIONS.bUseMedian, 1);
        if OPTIONS.bScatterLog
            xlabel('HaloTag Expression; log10(dTomato)', 'fontweight', 'bold', 'fontsize', FS);
        else
            xlabel('HaloTag Expression; (dTomato)', 'fontweight', 'bold', 'fontsize', FS);
        end
        set(gca,'ylim', OPTIONS.YLIM);
        grid on;
        ylabel(['Assay Activity (' num2str(OPTIONS.DoseVectPlot(q)/1000) ' uM)'], 'fontweight', 'bold', 'fontsize', FS);
        set(gca, 'fontsize', FS);
        
        for w=1:2
            if w==1
                subplot(2,2,3);  %HT-
                I = Series(k).Ineg; %HTneg for this dose only
                title('HT-');
            else
                subplot(2,2,4);  %HT+
                I = Series(k).Ipos; %HTpos
                title('HT+');
            end
            if size(REF,2) == 1
                R=REF(I,:);
            elseif size(REF,2)==size(RESP,2)
                R=REF(I,find(Series(k).DOSE == OPTIONS.DoseVectPlot(q)));
            else
                'debug'
            end
            Y=RESP(I, find(Series(k).DOSE == OPTIONS.DoseVectPlot(q)));
            R = R(~isnan(Y));
            Y = Y(~isnan(Y)); %must do this last
            PlotBinnedData(OPTIONS.REFBINS, R,Y,[], [1 0], clrs{k}, 'o', 1, OPTIONS.bUseMedian, 1);
            set(gca,'ylim', OPTIONS.YLIM);
            set(gca,'xlim', OPTIONS.REFLIM);
            hold on; 
            %drawnow;
            if bLog10
                xlabel(['log10(' OPTIONS.REFTXT ')']);
            else
                xlabel(OPTIONS.REFTXT);
            end
            grid on;
        end
        
    end
    
    %drawnow;
end

for s=1:length(Series)
    disp(Series(s).Name); 
    disp(['HTpos #CS=' num2str(min(Series(s).DRpos.nCS)) ', #Cells=' num2str(min(Series(s).DRpos.n)) ', N=' num2str(min(Series(s).DRpos.N))]);  
    disp(['HTneg #CS=' num2str(min(Series(s).DRneg.nCS)) ', #Cells=' num2str(min(Series(s).DRneg.n)) ', N=' num2str(min(Series(s).DRneg.N))]);  
    disp(newline); 
end

return;


