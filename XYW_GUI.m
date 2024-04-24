function varargout = XYW_GUI(varargin)
%
%
%to erase prior analysis without erasing rotation:
%   for k=1:12, for j=1:2, delete(['CS' num2str(k) '_' num2str(j) '*.mat']); end; end; delete XYW_*.mat;

if nargin == 0
    LaunchGUI('');
elseif nargin==1 && (strcmpi(varargin{1},'XYW') || strcmpi(varargin{1},'GH0')  || strcmpi(varargin{1},'GH1') || strcmpi(varargin{1},'GH2')  || strcmpi(varargin{1},'GH3') || strcmpi(varargin{1},'GH4') || strcmpi(varargin{1},'GH5') || strcmpi(varargin{1},'GH5.2') || strcmpi(varargin{1},'EF1') || strcmpi(varargin{1},'EF2')|| strcmpi(varargin{1},'EF2_WTF'))
    LaunchGUI(varargin{1});
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch ERR
        rethrow(ERR);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LaunchGUI(MODE)

%set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultFigureWindowStyle','normal')

%ERROR CHECK
if ~isempty(whos('global','XYW'))
    Title = 'Warning: is XYW_GUI already running?';
    Question = 'Only one instance of XYW_GUI can run at a time.  If you are sure one is not already running press OK';
    ButtonName=questdlg(Question,Title,'OK','Cancel','Cancel');
    if ~strcmpi(ButtonName, 'OK')
        return;
    end
end
DIRdRXY  = dir('*dRXY.mat');
if isempty(DIRdRXY)
    error(sprintf('Error: must perform alignment before analysis. Here is example code: \n  EF1/EF2:   >> for n=1:8, EF_ALIGN_IX83(n); end;  \n  EF2_WTF:   >> for n=1:8, EF_ALIGN_WTF(n); end;  \n  Otherwise: >> for n=1:12, XYW_ALIGN_IX83(n); end;'));
    return;
end

%INITIALIZE
global XYW
XYW = [];
XYW.PWD = pwd;
XYW.CS = DIRdRXY;
XYW.bEditMode = 1;
XYW.XYPlotMode = 2;
XYW.showROI = 1;

%load GlobalVars
try
    GV = load([XYW.PWD '\XYW_GlobalVars.mat']); % 'NR', 'NG', 'RFPReps', 'GFPReps' and several others ... see below
    Fld = fieldnames(GV);
    for f=1:length(Fld)
        XYW.(Fld{f}) = GV.(Fld{f});
    end
catch
    disp('No GlobalVars.mat');
end
if isempty(MODE) & ~isfield(XYW,'MODE')
    %not specified, not previously saved
    clear global XYW;
    error( [...
        'New analysis folder; Please specify XYW_GUI(''XYW'')  or  another mode (see below)' newline ...
        newline ...
        '% XYW:  very similar to GH3, except the segmentation parameters are less sensitive' newline ...
        newline ...
        '% GH0:  (DISCONTINUED) uses GH protocol (escalating LED intensity); no GABA during baseline (nREF = 1)' newline ...
        '% GH1:  (DISCONTINUED) uses GH protocol (escalating LED intensity); GABA during baseline (nREF = last)' newline ...
        '%*GH2:  *BEST FOR GABAZINE-DART; uses XYW protocol (fixed LED intensity     ); GABA during baseline (nREF = last)' newline ...
        '% GH3:  *BEST FOR DIAZEPAM-DART; uses XYW protocol (fixed LED intensity     ); no GABA during baseline (nREF = 1)' newline ...
        '% GH4:  *** same as GH2 but with 4 reps per dose (instead of 6)' newline ...
        '%*GH5:  uses GH protocol (DE-escalating LED intensity); highest activity at start of assay (nREF = 1)' newline ...
        '%*GH5.2:  uses GH protocol (DE-escalating LED intensity); highest activity at start of assay (nREF = 2)' newline ...
        ] );

elseif ~isempty(MODE) & ~isfield(XYW,'MODE')
    %yes specified, not previously saved
    XYW.MODE = MODE; %define it
elseif isempty(MODE) & isfield(XYW,'MODE')
    %not specified, yes previously saved
    %nothing to do here
elseif ~isempty(MODE) & isfield(XYW,'MODE')
    %yes specified, yes previously saved
    %confirm they match
    if ~strcmpi(MODE, XYW.MODE)
        MODE2 =  XYW.MODE;
        clear global XYW;
        error(['XYW_GUI(''' MODE ''')  called, but saved analysis was with XYW_GUI(''' MODE2 ''')']);
    end
end

%sort this folder according to the first number in each name
for j=1:length(XYW.CS)
    %for sorting purposes
    XYW.CS(j).Num = NumbersInString(XYW.CS(j).name);
    XYW.CS(j).Num = XYW.CS(j).Num(1);

    %determine stage of analysis based on existing .mat files
    RootName = ['CS' num2str(XYW.CS(j).Num) '_'];
    Fname = [RootName '1AutoSegment.mat'];
    if ~isempty(dir([XYW.PWD '\' Fname]))
        XYW.CS(j).name = Fname;
    end
    Fname = [RootName '2RoiProofread.mat'];
    if ~isempty(dir([XYW.PWD '\' Fname]))
        TMP = load([XYW.PWD '\' Fname]);  %load data for this one
        XYW.CS(j).DAT.Blobs = TMP.Blobs;
        if isempty(XYW.CS(j).DAT.Blobs)
            XYW.CS(j).name = ['CS' num2str(XYW.CS(j).Num) '_ no blobs'];
        elseif ~isfield(XYW.CS(j).DAT.Blobs,'bViewed')
            XYW.CS(j).name = ['CS' num2str(XYW.CS(j).Num) '_ no bViewed Flag!'];
        else
            PercentViewedThisCS = 100*sum(([XYW.CS(j).DAT.Blobs.bViewed]>0))/length(XYW.CS(j).DAT.Blobs);
            XYW.CS(j).name = ['CS' num2str(XYW.CS(j).Num) '_' num2str(floor(PercentViewedThisCS)) '% Viewed'];
        end
    else
        %data initialization
        XYW.CS(j).DAT = [];
    end
end
[~,I] = sort([XYW.CS.Num],'ascend');
XYW.CS = XYW.CS(I);



%create the main figure -- default position
% DFW = get(0,'DefaultFigureWindowStyle');
% set(0,'DefaultFigureWindowStyle','normal');
fig = figure('Name', 'XYW_GUI', 'NumberTitle', 'off', 'menubar', 'none');
% set(0,'DefaultFigureWindowStyle',DFW);

set(gcf, 'Pointer', 'watch'); drawnow;

%make handle visibility full so that we can create dynamic uicontrols in this function
%then set back to 'callback' at end of initialization
set(fig, 'HandleVisibility', 'on', 'IntegerHandle', 'on');

%Set the CloseRequestFcn -- to save options on close
set(fig, ...
    'CloseRequestFcn', CallbackString('XYW_GUI', '''MainFigure_CloseRequestFcn'''), ...
    'ResizeFcn', CallbackString('XYW_GUI', '''MainFigure_ResizeFcn'''), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));


%Now create all the main UI objects
XYW.UI = []; %place holder
XYW.bRGB = [true true true]; %toggle to show red/grn/blue
XYW.RGB = [];  %RGB image
XYW.hRGB = []; %handle to image in main axes
XYW.hRGBz = []; %handle to image in zoom axes

XYW.UI.AxIMG = axes('Units', 'pixels'); %image axes
XYW.UI.AxMAG = axes('Units', 'pixels'); %image axes
XYW.UI.AxWAV = axes('Units', 'pixels'); %waveform axes
XYW.UI.AxXVY = axes('Units', 'pixels'); %X vs Y axes

%listbox for each ROI
XYW.UI.ListBoxROI = uicontrol('Style', 'listbox', 'String', '', 'Units', 'pixels', 'max', inf, 'min', 0, ...
    'Callback', CallbackString('XYW_GUI','''ListBoxROI_Callback''', 0), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));
%listbox for each coverslip
XYW.UI.ListBoxCS = uicontrol('Style', 'listbox', 'String', {XYW.CS.name}, 'Units', 'pixels', ...
    'Callback', CallbackString('XYW_GUI','''ListBoxCS_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));

%checkboxes, edits, buttons
XYW.UI.EditModeCheck = uicontrol('Style', 'checkbox', 'String', 'EditMode', 'Value', XYW.bEditMode, 'Units', 'pixels', 'HorizontalAlignment', 'right', ...
    'Callback', CallbackString('XYW_GUI','''EditModeCheck_Callback''',1), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));

XYW.UI.ButtonToggleROI = uicontrol('Style', 'pushbutton', 'String', 'Toggle ROI (enter)', 'Units', 'pixels', ...
    'Callback', CallbackString('XYW_GUI','''ButtonToggleROI_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));

XYW.UI.ButtonNewROI = uicontrol('Style', 'pushbutton', 'String', 'New ROI (quote)', 'Units', 'pixels', ...
    'Callback', CallbackString('XYW_GUI','''ButtonNewROI_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));

XYW.UI.ButtonAutoSegmentAll = uicontrol('Style', 'pushbutton', 'String', 'Auto Segment All', 'Units', 'pixels', ...
    'Callback', CallbackString('XYW_GUI','''ButtonAutoSegmentAll_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));

XYW.UI.ButtonFinalAnalyzeAll = uicontrol('Style', 'pushbutton', 'String', 'Final Analyze All', 'Units', 'pixels', ...
    'Callback', CallbackString('XYW_GUI','''ButtonFinalAnalyzeAll_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));

XYW.UI.ButtonSaveBackup = uicontrol('Style', 'pushbutton', 'String', 'Save', 'Units', 'pixels', ...
    'Callback', CallbackString('XYW_GUI','''XYW_GUI_Save'''), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));

XYW.UI.XYPlotMode = uicontrol('Style', 'popupmenu', 'String', {'Grn/Uniformity' 'Red/Uniformity' 'Red/Grn' 'Blue/Uniformity' 'Blue/Grn' 'Blue/Red'}, 'Value', XYW.XYPlotMode, 'Units', 'pixels', ...
    'Callback', CallbackString('XYW_GUI','''XYPlotMode_Callback'''), ...
    'KeyPressFcn', @(obj, evd)  XYW_GUI('MainFigure_KeyPressFcn', obj, evd));

MainFigure_ResizeFcn;
EditModeCheck_Callback(0);
set(gcf, 'Pointer', 'arrow'); drawnow;


%make handle visibility full so that we can create dynamic uicontrols in this function
%then set back to 'callback' at end of initialization
%set(gcf, 'HandleVisibility', 'callback', 'IntegerHandle', 'off');

uicontrol(XYW.UI.ListBoxCS); %set focus to CS box (only at startup, so that enter key does something)

if isfield(XYW,'nSel')
    %loaded from GlobalVars
    set(XYW.UI.ListBoxCS,'Value',XYW.nSel);
    XYW = rmfield(XYW,'nSel');
    ListBoxCS_Callback;
end
if isfield(XYW,'nROI')
    %loaded from GlobalVars
    set(XYW.UI.ListBoxROI,'Value',XYW.nROI);
    XYW = rmfield(XYW,'nROI');
    ListBoxROI_Callback(0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MainFigure_ResizeFcn
%Called whenever the figure is resized

global XYW

%Layout constants
PosMain = get(gcf,'Position');
FigW = PosMain(3);
FigH = PosMain(4);

AxImgH = FigH;
AxImgW = min(AxImgH, FigW*0.6);  %ideally square, but less wide if need be
set(XYW.UI.AxIMG,'Position',     [1   1    AxImgW AxImgH], 'visible', 'on');

AxWavW = FigW - AxImgW;
AxWavH = round(AxWavW/10); %1:10 proportion
set(XYW.UI.AxWAV,'Position',[AxImgW FigH-AxWavH AxWavW AxWavH], 'visible', 'on');

LBRoiW = min(FigW*0.2, 200);
MagWH = LBRoiW*2;
LBRoiH = FigH - AxWavH - MagWH;
set(XYW.UI.ListBoxROI,'Position',[AxImgW+2 1 LBRoiW LBRoiH]);
set(XYW.UI.AxMAG,'Position',     [AxImgW+2 LBRoiH MagWH MagWH]);

LBcsH = min((FigH - AxWavH)/2, 300);
LBcsW = min(FigW - AxImgW - LBRoiW, 200);
set(XYW.UI.ListBoxCS,            'Position', [AxImgW+LBRoiW+2  1             LBcsW  LBcsH], 'fontsize', 10, 'fontweight', 'bold');
set(XYW.UI.ButtonSaveBackup,     'Position', [AxImgW+LBRoiW+2  LBcsH+3       LBcsW  36   ], 'fontsize', 10);
set(XYW.UI.ButtonFinalAnalyzeAll,'Position', [AxImgW+LBRoiW+2  LBcsH+6+36    LBcsW  36   ], 'fontsize', 10);
set(XYW.UI.ButtonAutoSegmentAll, 'Position', [AxImgW+LBRoiW+2  LBcsH+9+72    LBcsW  36   ], 'fontsize', 10);
set(XYW.UI.ButtonNewROI,         'Position', [AxImgW+LBRoiW+2  LBcsH+12+108  LBcsW  36   ], 'fontsize', 10);
set(XYW.UI.ButtonToggleROI,      'Position', [AxImgW+LBRoiW+2  LBcsH+15+144  LBcsW  36   ], 'fontsize', 10);
set(XYW.UI.EditModeCheck,        'Position', [AxImgW+LBRoiW+10 LBcsH+18+180  100  36   ], 'fontsize', 12, 'fontweight', 'bold');

set(XYW.UI.XYPlotMode,           'Position', [FigW-125       1            125     30   ], 'fontsize', 12);


S = 60; %for axes labels
AxXvyW = FigW - AxImgW - LBRoiW - LBcsW;
AxXvyH = FigH - AxWavH - S;
AxXvyW = max(10, AxXvyW);
AxXvyH = max(10, AxXvyH);
AxXvyW = min(AxXvyH/1.8, AxXvyW);
set(XYW.UI.AxXVY,'Position',[FigW-AxXvyW S AxXvyW AxXvyH], 'visible', 'on');
axis(XYW.UI.AxXVY, 'equal'); 
try
   axis(XYW.UI.AxXVY, XYW.UI.AxXVYminmax); %, = [XMinMax YMinMax];  %for resize
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EditModeCheck_Callback(bFullRefresh)
global XYW

XYW.bEditMode = get(XYW.UI.EditModeCheck, 'value');
if XYW.bEditMode
    set(XYW.UI.ButtonToggleROI,      'Enable', 'on');
    set(XYW.UI.ButtonNewROI,         'Enable', 'on');
    set(XYW.UI.ButtonAutoSegmentAll, 'Enable', 'on');
    set(XYW.UI.ButtonFinalAnalyzeAll,'Enable', 'off');
else
    set(XYW.UI.ButtonToggleROI,      'Enable', 'off');
    set(XYW.UI.ButtonNewROI,         'Enable', 'off');
    set(XYW.UI.ButtonAutoSegmentAll, 'Enable', 'off');
    set(XYW.UI.ButtonFinalAnalyzeAll,'Enable', 'on');
end
if bFullRefresh
    ListBoxCS_Callback;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XYPlotMode_Callback
global XYW

XYW.XYPlotMode = get(XYW.UI.XYPlotMode, 'value');
cla(XYW.UI.AxXVY); %signal that these handles all need to be redrawn
ListBoxCS_Callback;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonAutoSegmentAll_Callback
global XYW

Str = get(XYW.UI.ListBoxCS,'String');
for nSel = 1:length(Str)
    set(XYW.UI.ListBoxCS,'Value',nSel);
    set(XYW.UI.ButtonAutoSegmentAll, 'String', ['**Segmenting: ' num2str(nSel) '/' num2str(length(Str)) ' ...'], 'ForegroundColor', 'r', 'fontweight', 'bold');
    drawnow;
    ListBoxCS_Callback;
end
set(XYW.UI.ButtonAutoSegmentAll, 'String', 'Auto Segment All', 'ForegroundColor', 'k', 'fontweight', 'normal');
XYW_GUI_Save;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonFinalAnalyzeAll_Callback
global XYW

Str = get(XYW.UI.ListBoxCS,'String');
nSelLast = get(XYW.UI.ListBoxCS,'Value'); %last productive CS
for nSel = 1:length(Str)
    set(XYW.UI.ListBoxCS,'Value',nSel);
    set(XYW.UI.ButtonFinalAnalyzeAll, 'String', ['**Analyzing: ' num2str(nSel) '/' num2str(length(Str)) ' ...'], 'ForegroundColor', 'r', 'fontweight', 'bold');
    bAnalysisDone = FullAnalysis;
    if bAnalysisDone
        nSelLast = nSel;
        ListBoxCS_Callback;
    end
    drawnow;
end
set(XYW.UI.ListBoxCS,'Value',nSelLast);
set(XYW.UI.ButtonFinalAnalyzeAll, 'String', 'Final Analyze All', 'ForegroundColor', 'k', 'fontweight', 'normal');
XYW_GUI_Save;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ListBoxCS_Callback

set(gcf, 'Pointer', 'watch'); drawnow;
    
global XYW
    
nSel = get(XYW.UI.ListBoxCS,'Value');

if isfield(XYW,'nSel') && (nSel == XYW.nSel)  && ~isempty(XYW.CS(nSel).DAT)
    nROI = get(XYW.UI.ListBoxROI,'Value');
    BlobsPrevSel = XYW.CS(nSel).DAT.Blobs(nROI); %to keep same selection as before
    bFullRedraw = 0;
else
    BlobsPrevSel = []; %default
    bFullRedraw = 1;
end
XYW.nSel = nSel; %remember for next time

if isempty(XYW.CS(nSel).DAT) || ~isfield(XYW.CS(nSel).DAT, 'Image')   
    bFullRedraw = 1; %need to load data, so trigger full redraw
    
    RootName = ['CS' num2str(XYW.CS(nSel).Num) '_'];
    
    %Load/perform Auto Segmentation
    Fname = [RootName '1AutoSegment.mat'];
    if isempty(dir([XYW.PWD '\' Fname]))
        DoThatThing(RootName, XYW.MODE);
    end    
    WB = waitbar(0,'');
    waitbar(0.5, WB, 'Loading 1AutoSegment ...');        
    try
        ProofreadBlobs = XYW.CS(nSel).DAT.Blobs;
    catch
        ProofreadBlobs = [];
    end
    XYW.CS(nSel).DAT = load([XYW.PWD '\' Fname]);  %load data for this one
    %XYW.CS(nSel).name = Fname;
    %Str{nSel} = Fname;
    %set(XYW.UI.ListBoxCS,'String', Str);

    %Load/perform Auto Segmentation
    Fname = [RootName '2RoiProofread.mat'];
    if ~isempty(dir([XYW.PWD '\' Fname]))
        if ~isempty(ProofreadBlobs)
            XYW.CS(nSel).DAT.Blobs = ProofreadBlobs;
        else
            waitbar(1, WB, 'Loading 2RoiProofread ...');
            TMP = load([XYW.PWD '\' Fname]);  %load data for this one
            XYW.CS(nSel).DAT.Blobs = TMP.Blobs;
        end
        %XYW.CS(nSel).name = Fname;
        %Str{nSel} = Fname;
        %set(XYW.UI.ListBoxCS,'String', Str);
    end
    close(WB);    
end

%pull out some useful variables
%GFPReps = XYW.CS(nSel).DAT.AcqSettings.GFPReps;
%RFPReps = XYW.CS(nSel).DAT.AcqSettings.RFPReps;
%NumRepsPerDose = XYW.CS(nSel).DAT.AcqSettings.NumRepsPerDose;
%nDoseForSegmentation = XYW.CS(nSel).DAT.SegSettings.nDoseForSegmentation;
%nRepForSegmentation = XYW.CS(nSel).DAT.SegSettings.nRepForSegmentation;
%NG = [1:GFPReps/RFPReps] + (nDoseForSegmentation-1)*NumRepsPerDose*(GFPReps/RFPReps) + (nRepForSegmentation-1)*(GFPReps/RFPReps);

%categorize Blobs into good/bad
for k=1:length(XYW.CS(nSel).DAT.Blobs)
    if ~isfield(XYW.CS(nSel).DAT.Blobs(k),'bGood') || isempty(XYW.CS(nSel).DAT.Blobs(k).bGood) 
        %good/bad has not been determined previously
        XYW.CS(nSel).DAT.Blobs(k).bGood = true; %(D(k) >= 0); %default Good/Bad flag
    end
    if ~isfield(XYW.CS(nSel).DAT.Blobs(k),'bViewed') || isempty(XYW.CS(nSel).DAT.Blobs(k).bViewed) 
        %bViewed has not been annotated previously
        XYW.CS(nSel).DAT.Blobs(k).bViewed = false; %keep track of what's been viewed for QC
    end
end

%sort from best to worst
if isempty(XYW.CS(nSel).DAT.Blobs)
    X = [];
    Y = [];
    D = [];
    XMinMax = [ 0  1]; XName = '';
    YMinMax = [ 0  1]; YName = '';
    LineXY = {[0 1], [0 1]};
    XY = vertcat(LineXY{:}); 
else
    Xg = log10([XYW.CS(nSel).DAT.Blobs.delG]./[XYW.CS(nSel).DAT.Blobs.delGstd]);  %green Spatial Uniformity
    try
        Xr = log10([XYW.CS(nSel).DAT.Blobs.red1]./[XYW.CS(nSel).DAT.Blobs.red1Std]);  %red Spatial Uniformity
    catch
        Xr = log10([XYW.CS(nSel).DAT.Blobs.red]./[XYW.CS(nSel).DAT.Blobs.redStd]);  %red Spatial Uniformity
    end
    try
        Xb = log10([XYW.CS(nSel).DAT.Blobs.blue]./[XYW.CS(nSel).DAT.Blobs.blueStd]);  %blue Spatial Uniformity
        X = (Xg + Xr + Xb)/3;
    catch
        X = (Xg + Xr)/2;
    end
    X(isinf(X)) = 1;
    X(X>1.2) = 1.2;
    XMinMax = [-0.2  1.22]; XName = 'Spatial Uniformity';
    if ~isempty(findstr(XYW.MODE,'EF'))
        MinMaxG = [ 2.6  4.5]+0.5;
    else
        MinMaxG = [ 2.6  4.5];
    end
    if XYW.XYPlotMode==1
        Y = log10([XYW.CS(nSel).DAT.Blobs.delG]);                                       %SOMA delG (log)
        YMinMax = MinMaxG; YName = 'delG Magnitude (log)';
        LineXY = {[0.8 2],[0 4]};
    elseif XYW.XYPlotMode==2
        Y = log10([XYW.CS(nSel).DAT.Blobs.red]);                                       %SOMA Red (log)
        YMinMax = [ 1.9  5.1]; YName = 'Red Magnitude (log)';
        LineXY = {[0.8 6], [0 2]};
    elseif XYW.XYPlotMode==3
        X = log10([XYW.CS(nSel).DAT.Blobs.delG]);                                       %SOMA delG (log)
        Y = log10([XYW.CS(nSel).DAT.Blobs.red]);                                       %SOMA Red (log)
        XMinMax = MinMaxG; XName = 'delG Magnitude (log)';
        YMinMax = [ 1.9  5.1]; YName = 'Red Magnitude (log)';
        LineXY = {[2 4], [10 2]};
    elseif XYW.XYPlotMode==4
        Y = log10([XYW.CS(nSel).DAT.Blobs.red]);                                       %SOMA Red (log)
        YMinMax = [ 1.9  5.1]; YName = 'Blue Magnitude (log)';
        LineXY = {[0.8 6], [0 2]};
    elseif XYW.XYPlotMode==5
        X = log10([XYW.CS(nSel).DAT.Blobs.delG]);                                       %SOMA delG (log)
        Y = log10([XYW.CS(nSel).DAT.Blobs.blue]);                                       %SOMA Red (log)
        XMinMax = MinMaxG; XName = 'delG Magnitude (log)';
        YMinMax = [ 1.9  5.1]; YName = 'Blue Magnitude (log)';
        LineXY = {[2 4], [10 2]};
    elseif XYW.XYPlotMode==6
        Y = log10([XYW.CS(nSel).DAT.Blobs.red]);                                       %SOMA Red (log)
        X = log10([XYW.CS(nSel).DAT.Blobs.blue]);                                       %SOMA blue (log)
        YMinMax = [ 1.9  5.1]; YName = 'Red Magnitude (log)';
        XMinMax = [ 1.9  3.6]; XName = 'Blue Magnitude (log)';
        LineXY = {[0 4],[4 7]};
    end
    D = point_to_line([X(:) Y(:)], LineXY{:});  %distance from good/bad line
    XY = vertcat(LineXY{:}); P = polyfit(XY(:,1),XY(:,2),1);
    I = find(Y <= (P(1)*X + P(2))); D(I) = -1*D(I);  %negative
    if XYW.bEditMode==1
        J = find([XYW.CS(nSel).DAT.Blobs.bViewed]==2);  %hand-drawn go last
    else
        J = find([XYW.CS(nSel).DAT.Blobs.bGood]==0);  %bad go last
    end
    if XYW.XYPlotMode==1
        D(J) = D(J) + 10*max(abs(D));
        [~,I] = sort(D,'ascend');
    else
        D(J) = D(J) - 10*max(abs(D));
        [~,I] = sort(D,'descend');
    end
    XYW.CS(nSel).DAT.Blobs = XYW.CS(nSel).DAT.Blobs(I);
    X = X(I);
    Y = Y(I);
    D = D(I);
end

%generate strings for the ROI list box
UpdateROItext;
nSelROI = 1; %find(vertcat(XYW.CS(nSel).DAT.Blobs.bGood));
set(XYW.UI.ListBoxROI, 'value', nSelROI);

%grnDenominator = mean(vertcat(XYW.CS(nSel).DAT.Blobs(I).grnMax)) - mean(vertcat(XYW.CS(nSel).DAT.Blobs(I).grnMin));

if size(XYW.CS(nSel).DAT.SegOutput.SegmentImage, 3)==3
    bNewBLUE = 1;
else
    bNewBLUE = 0;    
end
if bFullRedraw
    %RGB = cat(3, 1*double(XYW.CS(nSel).DAT.Image.RedMean), double(XYW.CS(nSel).DAT.Image.GrnDelta), 0*double(XYW.CS(nSel).DAT.Image.GrnMean));
    if bNewBLUE
        XYW.RGB = XYW.CS(nSel).DAT.SegOutput.SegmentImage;
        for CH=1:3
            if CH==1 %Red
                Min = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(1);
                Max = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(2)/3;  %boost red
                Gam = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(3);
                Sat = 0.8;
            elseif CH==2 %green
                Min = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(1);
                Max = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(2);
                Gam = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(3);
                Sat = 0.8;
            elseif CH==3 %blue
                %CH==3;  
                Min = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{1}(1)/2;
                Max = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{1}(2)/4;  %super-boost blue
                Gam = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{1}(3);
                Sat = 1;
            end
            TMP = (XYW.RGB(:,:,CH)-Min)/(Max-Min);
            TMP(TMP<0)=0;
            TMP(TMP>1)=1;
            TMP=TMP.^Gam;
            XYW.RGB(:,:,CH) = TMP*Sat; %prevent saturation
        end
    else
        XYW.RGB = cat(3, 1*double(XYW.CS(nSel).DAT.SegOutput.SegmentImage(:,:,1)), double(XYW.CS(nSel).DAT.SegOutput.SegmentImage(:,:,2)), 0*double(XYW.CS(nSel).DAT.Image.GrnMean));
        for CH=1:2
            Min = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(1);
            Max = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(2);
            Gam = XYW.CS(nSel).DAT.SegSettings.MinMaxGamSOMA{CH}(3);
            TMP = (XYW.RGB(:,:,CH)-Min)/(Max-Min);
            TMP(TMP<0)=0;
            TMP(TMP>1)=1;
            TMP=TMP.^Gam;
            XYW.RGB(:,:,CH) = TMP*0.8; %prevent saturation
        end
        %XYW.RGB = XYW.RGB/(4000);
        %XYW.RGB(XYW.RGB>1) = 1;
        %XYW.RGB(XYW.RGB<0) = 0;
        %XYW.RGB = XYW.RGB.^0.8;
        %XYW.RGB = XYW.RGB*0.8; %prevent saturation
        %XYW.RGB(repmat(XYW.CS(nSel).DAT.SegOutput.ProcessImBinary==1,1,1,3)) = 0.5;
        XYW.RGB(:,:,3) = 0.5*XYW.CS(nSel).DAT.SegOutput.ProcessImBinary;
    end
    
    TMP = XYW.RGB;
    for k=1:3
        if ~XYW.bRGB(k)
            TMP(:,:,k) = 0;
        end
    end
    cla(XYW.UI.AxIMG);
    cla(XYW.UI.AxMAG);
    cla(XYW.UI.AxXVY);
    XYW.hRGB  = image(XYW.UI.AxIMG, TMP);
    XYW.hRGBz = image(XYW.UI.AxMAG, TMP ); %zoom
    axis(XYW.UI.AxIMG, 'image'); axis(XYW.UI.AxIMG, 'off');
    axis(XYW.UI.AxMAG, 'image'); axis(XYW.UI.AxMAG, 'off');
    set(XYW.hRGB,'ButtonDownFcn', CallbackString('XYW_GUI','''CenterMagToClick''', 1));
    %set(XYW.hRGBz,'ButtonDownFcn', CallbackString('XYW_GUI','''CenterMagToClick''', 2));
end

%set(gca,'position',[0 0 1 1]);
ImH = size(XYW.RGB,1);
ImW = size(XYW.RGB,2);
NumROIRedrawn = 0;
for k=1:length(XYW.CS(nSel).DAT.Blobs)
    try
        if all(isvalid(XYW.CS(nSel).DAT.Blobs(k).h(1:5)))
            bNeedNewROI = 0;
        else
            bNeedNewROI = 1;
        end
    catch
        bNeedNewROI = 1;
    end
    if bNeedNewROI
        NumROIRedrawn = NumROIRedrawn + 1;
        %delete any old handles
        if isfield(XYW.CS(nSel).DAT.Blobs(k),'h')
            for j=1:length(XYW.CS(nSel).DAT.Blobs(k).h)
                try
                    delete(XYW.CS(nSel).DAT.Blobs(k).h(j));
                catch
                end
            end
        end
        
        %need to generate graphics handles
        %crosshair
        xy    = mean(XYW.CS(nSel).DAT.Blobs(k).Boundary,1);
        xyMin = min(XYW.CS(nSel).DAT.Blobs(k).Boundary,[],1);
        xyMax = max(XYW.CS(nSel).DAT.Blobs(k).Boundary,[],1);
        
        %note use .h = to set the type to 'Line' (avoid converting to double)
        XYW.CS(nSel).DAT.Blobs(k).h = line(XYW.UI.AxIMG, [xy(1) xy(1) NaN xy(1) xy(1)   NaN   1 xyMin(1) NaN xyMax(1) ImW],[1 xyMin(2) NaN xyMax(2) ImH   NaN   xy(2) xy(2) NaN xy(2) xy(2)]);
        set(XYW.CS(nSel).DAT.Blobs(k).h(1),'LineStyle',':');
        
        %ROI normal axis
        XYW.CS(nSel).DAT.Blobs(k).h(2) = line(XYW.UI.AxIMG, XYW.CS(nSel).DAT.Blobs(k).Boundary([1:end 1],1),   XYW.CS(nSel).DAT.Blobs(k).Boundary([1:end 1],2));
        
        %ROI zoom axis
        XYW.CS(nSel).DAT.Blobs(k).h(3) = line(XYW.UI.AxMAG, XYW.CS(nSel).DAT.Blobs(k).Boundary([1:end 1],1),   XYW.CS(nSel).DAT.Blobs(k).Boundary([1:end 1],2));
        
        %Code for XvsY:
        Boundary = XYW.CS(nSel).DAT.Blobs(k).Boundary; %[x, y] x (many rows)
        SF = 0.002;
        Boundary(:,1) =  SF*(Boundary(:,1) - mean(Boundary(:,1))) + X(k);
        Boundary(:,2) = -SF*(Boundary(:,2) - mean(Boundary(:,2))) + Y(k);
        
        %crosshair
        xy    = mean(Boundary,1);
        xyMin = min(Boundary,[],1);
        xyMax = max(Boundary,[],1);
        m = P(1); %slope of the LineXY
        b = xy(2) - m*xy(1); %b = y - mx
        if abs(P(1))<1
            XYW.CS(nSel).DAT.Blobs(k).h(4) = line(XYW.UI.AxXVY, ...
                [xy(1)      xy(1)    NaN xy(1)    xy(1)           NaN        [XMinMax(1) xyMin(1) NaN xyMax(1) XMinMax(2)]     ], ...
                [YMinMax(1) xyMin(2) NaN xyMax(2) YMinMax(2)      NaN      m*[XMinMax(1) xyMin(1) NaN xyMax(1) XMinMax(2)]+b   ]);
            %vertical lines                                            %diagonal lines
        else
            XYW.CS(nSel).DAT.Blobs(k).h(4) = line(XYW.UI.AxXVY, ...
                [XMinMax(1) xyMin(1) NaN xyMax(1) XMinMax(2)   NaN     ([YMinMax(1) xyMin(2) NaN xyMax(2) YMinMax(2)]-b)/m  ], ...
                [xy(2)      xy(2)    NaN xy(2)    xy(2)        NaN      [YMinMax(1) xyMin(2) NaN xyMax(2) YMinMax(2)]       ]);
            %horizontal lines                                      %diagonal lines
        end
        set(XYW.CS(nSel).DAT.Blobs(k).h(4),'LineStyle',':');
        
        %ROI
        XYW.CS(nSel).DAT.Blobs(k).h(5) = line(XYW.UI.AxXVY, Boundary([1:end 1],1),  Boundary([1:end 1],2));
    end
    
    %always update ButtonDownFcn (order of Blobs may have changed)
    set(XYW.CS(nSel).DAT.Blobs(k).h([2 3 5]), 'ButtonDownFcn', CallbackString('XYW_GUI','''ROI_Click''',k));
        
    if XYW.showROI==1 || (XYW.showROI==2 && XYW.CS(nSel).DAT.Blobs(k).bGood)
        set(XYW.CS(nSel).DAT.Blobs(k).h([2 3]), 'visible', 'on');
    else
        set(XYW.CS(nSel).DAT.Blobs(k).h([2 3]), 'visible', 'off');
    end
end
disp([num2str(NumROIRedrawn) ' of ' num2str(length(XYW.CS(nSel).DAT.Blobs)) ' ROI Redrawn']);

set(XYW.UI.AxXVY, 'color', [0 0 0]+0.2);
axis(XYW.UI.AxXVY, 'equal'); axis(XYW.UI.AxXVY, 'on');
axis(XYW.UI.AxXVY, [XMinMax YMinMax]);
XYW.UI.AxXVYminmax = [XMinMax YMinMax];  %for resize
xlabel(XYW.UI.AxXVY, XName, 'fontsize', 20);
ylabel(XYW.UI.AxXVY, YName, 'fontsize', 20);

if ~isempty(BlobsPrevSel)
    %try to find matching Blobs
    nROI = [];
    for k=1:length(BlobsPrevSel)
        Pk = BlobsPrevSel(k).PixelIdxList;
        for j=1:length(XYW.CS(nSel).DAT.Blobs)
            Pj = XYW.CS(nSel).DAT.Blobs(j).PixelIdxList;
            if length(Pk)==length(Pj) && all(Pk==Pj)
                nROI = [nROI j]; %found a match for BlobsPrevSel(k)
                break;
            end
        end
    end
    set(XYW.UI.ListBoxROI,'Value', nROI);
    ListBoxROI_Callback(1); %this takes care of arrow/wait
else
    ListBoxROI_Callback(0); %this takes care of arrow/wait
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = point_to_line(pt, v1, v2)
% pt should be nx3
% v1 and v2 are vertices on the line (each 1x3)
% d is a nx1 vector with the orthogonal distances
v1 = repmat(v1,size(pt,1),1);
v2 = repmat(v2,size(pt,1),1);
a = v1 - v2;
b = pt - v2;
if size(a,2)<3
    a(:,3) = 0;
end
if size(b,2)<3
    b(:,3) = 0;
end
d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateROItext
global XYW
nSel = get(XYW.UI.ListBoxCS,'Value');

BlobStr = {};
for k=1:length(XYW.CS(nSel).DAT.Blobs)
    BlobStr{k} = [num2str(k,'%.3d') '--R' num2str(log10(XYW.CS(nSel).DAT.Blobs(k).red),'%.1f')  '--dG' num2str(log10(XYW.CS(nSel).DAT.Blobs(k).delG),'%.1f')];
    if ~XYW.CS(nSel).DAT.Blobs(k).bViewed
        BlobStr{k} = [BlobStr{k}  '==== NEED QC ===='];
    elseif XYW.CS(nSel).DAT.Blobs(k).bViewed==2
        BlobStr{k} = [BlobStr{k}  '---by hand'];
    end
    
    if ~XYW.CS(nSel).DAT.Blobs(k).bGood
        BlobStr{k} = [BlobStr{k}  '---x'];
    end
end
set(XYW.UI.ListBoxROI, 'String', BlobStr);

StrCS = get(XYW.UI.ListBoxCS,'String');
NumTot  = 0;
NumGood  = 0;
NumViewed = 0;
NumGoodViewed = 0;
for k=1:length(XYW.CS)
    try
        NumTot = NumTot + length(XYW.CS(k).DAT.Blobs);
        NumGood = NumGood + sum([XYW.CS(k).DAT.Blobs.bGood]);
        NumViewed = NumViewed + sum(([XYW.CS(k).DAT.Blobs.bViewed]>0));
        NumGoodViewed = NumGoodViewed + sum(([XYW.CS(k).DAT.Blobs.bGood].*([XYW.CS(k).DAT.Blobs.bViewed]>0)));
        
        PercentViewedThisCS = 100*sum(([XYW.CS(k).DAT.Blobs.bViewed]>0))/length(XYW.CS(k).DAT.Blobs);
        StrCS{k} = ['CS' num2str(k) '_' num2str(floor(PercentViewedThisCS)) '% Viewed'];
    catch
        %CS(k).DAT not yet loaded
    end
end
set(XYW.UI.ListBoxCS,'String',StrCS);
%set(gcf, 'Name', ['XYW_GUI; ' num2str(100*NumViewed/NumTot,'%.1f') '% Viewed: ' num2str(NumViewed) '/' num2str(NumTot) ' (' num2str(NumGoodViewed) '/' num2str(NumViewed) ' Good = ' num2str(100*NumGoodViewed/NumViewed,'%.1f') '%)']);
set(gcf, 'Name', ['XYW_GUI; (' num2str(NumGoodViewed) '/' num2str(NumViewed) ' Good = ' num2str(100*NumGoodViewed/NumViewed,'%.1f') '%)']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ROI_Click(nROInew)
global XYW
nSel = get(XYW.UI.ListBoxCS,'Value');

nROIold = get(XYW.UI.ListBoxROI,'Value');
set(XYW.UI.ListBoxROI,'Value', nROInew);
if length(nROInew)==length(nROIold) && all(nROInew==nROIold)
    if ~XYW.bEditMode
        waitfor(errordlg('Must be in Edit Mode to mark ROI good/bad'));
        return;
    end
    
    %toggle good/bad
    for k=nROInew
        if k==nROInew(1)
            %toggle the first one
            XYW.CS(nSel).DAT.Blobs(k).bGood = ~XYW.CS(nSel).DAT.Blobs(k).bGood;
        else
            %make the rest the same as the first
            XYW.CS(nSel).DAT.Blobs(k).bGood = XYW.CS(nSel).DAT.Blobs(nROInew(1)).bGood;
        end
    end
    UpdateROItext;
end

ListBoxROI_Callback(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ListBoxROI_Callback(bResetROIScroll)

global XYW

nSel = get(XYW.UI.ListBoxCS,'Value');
nROI = get(XYW.UI.ListBoxROI,'Value');

if isempty(XYW.CS(nSel).DAT.Blobs)
    uicontrol(XYW.UI.ListBoxCS); %set focus to CS ListBox
    set(gcf, 'Pointer', 'arrow'); drawnow;
    return; %nothing to do
else    
    %pull out some useful variables
    GFPReps = XYW.CS(nSel).DAT.AcqSettings.GFPReps;
    RFPReps = XYW.CS(nSel).DAT.AcqSettings.RFPReps;
    NumRepsPerDose = XYW.CS(nSel).DAT.AcqSettings.NumRepsPerDose;
    NumRuns        = XYW.CS(nSel).DAT.AcqSettings.NumRuns;
    NumDose        = XYW.CS(nSel).DAT.AcqSettings.NumDose;
    
    %nDoseForSegmentation = XYW.CS(nSel).DAT.SegSettings.nDoseForSegmentation;
    %nRepForSegmentation = XYW.CS(nSel).DAT.SegSettings.nRepForSegmentation;
    %NG = [1:GFPReps/RFPReps] + (nDoseForSegmentation-1)*NumRepsPerDose*(GFPReps/RFPReps) + (nRepForSegmentation-1)*(GFPReps/RFPReps);
    NG = XYW.NG; %indices used for getting data
    nTG = XYW.CS(nSel).DAT.AcqSettings.nTG;

    axes(XYW.UI.AxWAV); cla;
    bOffsetRed = 0;
    for k=nROI(end:-1:1)
        %grnNorm = (XYW.CS(nSel).DAT.Blobs(k).grnRaw-XYW.CS(nSel).DAT.Blobs(k).grnMin)./grnDenominator;
        if XYW.bEditMode
            if ~isempty(findstr(XYW.MODE,'EF'))
                %EF1 or EF2 -- non-standard
                Hide = 6:20;  %hard coded
                grnNorm = XYW.CS(nSel).DAT.Blobs(k).grnNorm;
                grnNorm(Hide) = NaN;
            else
                grnNorm = NaN*XYW.CS(nSel).DAT.Blobs(k).grnNorm;
                grnNorm(NG) = XYW.CS(nSel).DAT.Blobs(k).grnNorm(NG);
            end
        else
            grnNorm = XYW.CS(nSel).DAT.Blobs(k).grnNorm;
        end
        if bOffsetRed
            grnNorm = 0.06*grnNorm + log10(XYW.CS(nSel).DAT.Blobs(k).red);
        end
        if ~XYW.CS(nSel).DAT.Blobs(k).bGood
            if length(nROI)>10
                continue;
            else
                lw = 1;
                clr = [0 0 0]+0.8;
            end
        else
            if length(nROI)>1
                lw = 1;
                clr = MapValue2Color(log10(XYW.CS(nSel).DAT.Blobs(k).red),[2 4.5], jet);
            else
                lw = 2;
                clr = 'k'; %single trace
            end
        end
        plot(grnNorm, 'color', clr, 'linewidth', lw); hold on;
    end
    axis(XYW.UI.AxWAV, 'auto');
    ax = axis(XYW.UI.AxWAV);
    ylim = get(gca,'ylim');
    for k=0:GFPReps:nTG
        if mod(k,GFPReps*NumRepsPerDose)==0
            plot([k k],ylim, 'k', 'linewidth', 1);
        else
            plot([k k],ylim, ':k');
        end
    end
    if XYW.bEditMode
        axis(XYW.UI.AxWAV, ax);
    else
        axis(XYW.UI.AxWAV, [1 nTG ax(3:4)]);
    end
    set(XYW.UI.AxWAV, 'XTickLabel', {}, 'YTickLabel', {});
    if length(nROI)==1  && NumDose>1
        NperDose = length(grnNorm)/NumDose;
        
        %a linear fit through the minimum of each dose
        X = 1:NumDose;
        for q=1:length(X)
            GrnSort = sort(grnNorm( (1+(q-1)*NperDose):(q*NperDose) ));
            Y(q) = mean(GrnSort(1:ceil(0.05*length(GrnSort))));    %bottom 5%
        end
        P = polyfit(X,Y,1); %linear fit
        %for q=1:2 %eliminate the two midpoints that are the highest above the linear regression
        %    [~,I] = max(Y(2:end-1)-polyval(P,X(2:end-1)));
        %    X(I+1) = [];
        %    Y(I+1) = [];
        %    P = polyfit(X,Y,1); %linear fit
        %end
        grnMin = polyval(P, (1:length(grnNorm))/NperDose-0.5 );
        plot(grnMin, ':', 'color', clr, 'linewidth', 1);
        
    end
    
    if length(nROI)==1 && XYW.bEditMode
        if ~(XYW.CS(nSel).DAT.Blobs(nROI).bViewed)
            XYW.CS(nSel).DAT.Blobs(nROI).bViewed = 1;
            XYW.UI.topROI = get(XYW.UI.ListBoxROI, 'Listboxtop'); %save scroll position (may get messed up upon changing strings, so save it here)
            UpdateROItext;
        end
    end
    
    for k=1:length(XYW.CS(nSel).DAT.Blobs)
        if XYW.CS(nSel).DAT.Blobs(k).bGood
            set(XYW.CS(nSel).DAT.Blobs(k).h(1),'color',[1 1 1],'linewidth',2); %crosshair
            if XYW.CS(nSel).DAT.Blobs(k).bViewed
                set(XYW.CS(nSel).DAT.Blobs(k).h(2),'color',[1 1 1], 'linestyle', '-');
                set(XYW.CS(nSel).DAT.Blobs(k).h(3),'color',[1 1 1], 'linestyle', '-');
            else
                set(XYW.CS(nSel).DAT.Blobs(k).h(2),'color',[1 1 1], 'linestyle', ':');
                set(XYW.CS(nSel).DAT.Blobs(k).h(3),'color',[1 1 1], 'linestyle', '--');
            end
            set(XYW.CS(nSel).DAT.Blobs(k).h(4),'color',[1 1 1],'linewidth',2); %crosshair
            set(XYW.CS(nSel).DAT.Blobs(k).h(5),'color', MapValue2Color(log10(XYW.CS(nSel).DAT.Blobs(k).red),[2 4.5], jet) );
        else
            set(XYW.CS(nSel).DAT.Blobs(k).h(1),'color',[1 1 1],'linewidth',1); %crosshair
            set(XYW.CS(nSel).DAT.Blobs(k).h(2),'color',[0 0 0]);
            set(XYW.CS(nSel).DAT.Blobs(k).h(3),'color',[0 0 0]);
            set(XYW.CS(nSel).DAT.Blobs(k).h(4),'color',[1 1 1],'linewidth',1); %crosshair
            set(XYW.CS(nSel).DAT.Blobs(k).h(5),'color',[0 0 0]);
        end
        
        if any(k==nROI)
            set(XYW.CS(nSel).DAT.Blobs(k).h(1), 'visible', 'on');
            set(XYW.CS(nSel).DAT.Blobs(k).h(2), 'linewidth', 4);
            set(XYW.CS(nSel).DAT.Blobs(k).h(3), 'linewidth', 4);
            set(XYW.CS(nSel).DAT.Blobs(k).h(4), 'visible', 'on');
            set(XYW.CS(nSel).DAT.Blobs(k).h(5), 'linewidth', 4);
        else
            set(XYW.CS(nSel).DAT.Blobs(k).h(1), 'visible', 'off');
            set(XYW.CS(nSel).DAT.Blobs(k).h(2), 'linewidth', 1);
            set(XYW.CS(nSel).DAT.Blobs(k).h(3), 'linewidth', 1);
            set(XYW.CS(nSel).DAT.Blobs(k).h(4), 'visible', 'off');
            set(XYW.CS(nSel).DAT.Blobs(k).h(5), 'linewidth', 1);
        end
        
        if XYW.showROI==1 || (XYW.showROI==2 && XYW.CS(nSel).DAT.Blobs(k).bGood)
            set(XYW.CS(nSel).DAT.Blobs(k).h([2 3]), 'visible', 'on');
        else
            set(XYW.CS(nSel).DAT.Blobs(k).h([2 3]), 'visible', 'off');
        end
    end
    
    TMP = round(vertcat(XYW.CS(nSel).DAT.Blobs(nROI).Boundary));
    xyMid = round(mean(TMP,1));
    xyMin = min(TMP,[],1);
    xyMax = max(TMP,[],1);
    xyMin = min([xyMin-30; xyMid-100],[],1);
    xyMax = max([xyMax+30; xyMid+100],[],1);
    %xyMin(xyMin<1) = 1;
    %xyMax(xyMax>length(XYW.RGB)) = length(XYW.RGB);
    axis(XYW.UI.AxMAG, 'image'); axis(XYW.UI.AxMAG, 'off');
    try
        axis(XYW.UI.AxMAG, [xyMin(1) xyMax(1) xyMin(2)  xyMax(2)]); %zoom in appropriately
    catch
    end
    

    set(gcf, 'Pointer', 'arrow'); 
    uicontrol(XYW.UI.ListBoxROI); %set focus to ROIListBox
    if bResetROIScroll && ~isempty(XYW.UI.topROI)
        disp([get(XYW.UI.ListBoxROI, 'Listboxtop')  XYW.UI.topROI]); %for debugging
        set(XYW.UI.ListBoxROI, 'Listboxtop', XYW.UI.topROI); %set scroll position
        disp([get(XYW.UI.ListBoxROI, 'Listboxtop')  XYW.UI.topROI]); %for debugging
    else
        XYW.UI.topROI = get(XYW.UI.ListBoxROI, 'Listboxtop'); %save scroll position
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CenterMagToClick(Ax1Or2)
global XYW;

if Ax1Or2==1
    tmp = get(XYW.UI.AxIMG,'CurrentPoint');
else
    tmp = get(XYW.UI.AxMAG,'CurrentPoint');
end
xyMid = round(tmp(1,1:2));
xyMin = xyMid-100;
xyMax = xyMid+100;
axis(XYW.UI.AxMAG, [xyMin(1) xyMax(1) xyMin(2)  xyMax(2)]); %zoom in appropriately

if Ax1Or2==1
    %de-select all ROIs
    nSel = get(XYW.UI.ListBoxCS,'Value');
    nROI = get(XYW.UI.ListBoxROI,'Value'); %prior selection
    for j=nROI
        set(XYW.CS(nSel).DAT.Blobs(j).h([1 4]), 'visible', 'off');
        set(XYW.CS(nSel).DAT.Blobs(j).h([2 3 5]), 'linewidth', 1);
    end
    set(XYW.UI.ListBoxROI,'Value', []); %de-select all ROIs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonToggleROI_Callback
global XYW
if ~isempty(get(XYW.UI.ListBoxROI,'string'))
    XYW.UI.topROI = get(XYW.UI.ListBoxROI, 'Listboxtop'); %save scroll position (may get messed up upon changing strings, so save it here)
    nROI = get(XYW.UI.ListBoxROI,'Value');
    ROI_Click(nROI);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ButtonNewROI_Callback
global XYW

if ~XYW.bEditMode
    waitfor(errordlg('Must be in Edit Mode to draw New ROI'));
    return;
end

nSel = get(XYW.UI.ListBoxCS,'Value');
NumDose = XYW.CS(nSel).DAT.AcqSettings.NumDose;
NumRepsPerDose = XYW.CS(nSel).DAT.AcqSettings.NumRepsPerDose;

h = imellipse(XYW.UI.AxMAG);
if isempty(h)
    return; %user aborted
end
%wait(h);
xy = getVertices(h);
delete(h);

BW = poly2mask(xy(:,1), xy(:,2), size(XYW.RGB,1), size(XYW.RGB,2));

Blob.Centroid = mean(xy);
Blob.Eccentricity = 0;
Blob.FilledArea = sum(BW(:));
Blob.Solidity = 1;
Blob.PixelIdxList = find(BW);
Blob.Boundary = xy;

if Blob.FilledArea < 1
    waitfor(errordlg('Blob had zero area, so it was not added'));
    return;
end

GRN_NG = XYW.CS(nSel).DAT.Image.GRN_NG;
GRN = zeros(size(GRN_NG,1),size(GRN_NG,2),XYW.CS(nSel).DAT.AcqSettings.nTG, class(GRN_NG)); 
GRN(:,:,XYW.NG) = GRN_NG;
Blob = AnalyzeBlobs(Blob, XYW.CS(nSel).DAT.Image, GRN, XYW.CS(nSel).DAT.SegOutput.Ibkg, [], NumRepsPerDose, NumDose);

Blob.bGood = 1;
Blob.bViewed= 2; %special value to indicate it was manually drawn
Blob.h = [];

nROI = get(XYW.UI.ListBoxROI,'Value'); %prior selection
StrROI = get(XYW.UI.ListBoxROI,'String');
if ~isempty(nROI) & ~isempty(StrROI)
    nR = nROI(end); %cursor position; we will insert after this
    nROI(nROI>nR) = nROI(nROI>nR) + 1;  %we will insert a new ROI after nR
else
    nR = length(StrROI); %cursor position
    XYW.UI.topROI = []; %auto-scroll after refreshing everything
end
XYW.CS(nSel).DAT.Blobs = [XYW.CS(nSel).DAT.Blobs(1:nR); Blob; XYW.CS(nSel).DAT.Blobs(nR+1:end)]; %insert adjacent to currently selected ROI (for purposes of listbox scrolling headache)
StrROI = [StrROI(1:nR); {'new'}; StrROI(nR+1:end)];
set(XYW.UI.ListBoxROI,'String', StrROI);
set(XYW.UI.ListBoxROI,'Value', [nR+1 nROI]); %select this one, keeping the prior selection as the cursor position
ListBoxCS_Callback; %resort, redraw 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XYW_GUI_Save
global XYW;

WB = waitbar(0,'Saving ...');
for n=1:length(XYW.CS)
    if ~isempty(XYW.CS(n).DAT)
        if isempty(XYW.CS(n).DAT.Blobs)
            Blobs = [];
        elseif ~isfield(XYW.CS(n).DAT.Blobs, 'h')
            Blobs = XYW.CS(n).DAT.Blobs;  
        else
            Blobs = rmfield(XYW.CS(n).DAT.Blobs,'h'); %don't save handles
        end
        RootName = ['CS' num2str(XYW.CS(n).Num) '_'];
        %rr
        save([XYW.PWD '\' RootName '2RoiProofread.mat'], 'Blobs');
    end
    waitbar(n/length(XYW.CS), WB, ['Saved CS' num2str(n)]);
end
try
    %Save global Vars file
    NR = XYW.NR;
    NG = XYW.NG;
    RFPReps = XYW.RFPReps;
    GFPReps = XYW.GFPReps;
    MODE = XYW.MODE;
    nSel = get(XYW.UI.ListBoxCS,'Value');
    nROI = get(XYW.UI.ListBoxROI,'Value');
    bEditMode = XYW.bEditMode;
    XYPlotMode = XYW.XYPlotMode;
    
    save([XYW.PWD '\XYW_GlobalVars.mat'], 'NR', 'NG', 'RFPReps', 'GFPReps', 'MODE','nSel','nROI', 'bEditMode','XYPlotMode');
catch
    warning('XYW_GlobalVars.mat NOT SAVED!');
end

close(WB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MainFigure_CloseRequestFcn
global XYW;

try
    XYW_GUI_Save;
    clear global XYW;    
    closereq;  %standard matlab close function    
catch
    Title = 'Error closing';
    Question = 'Error closing; Data not saved.  If you are sure you want to exit press OK';
    ButtonName=questdlg(Question,Title,'OK','Cancel','Cancel');
    if ~strcmpi(ButtonName, 'OK')
        return;
    else
        closereq;  %standard matlab close function
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MainFigure_KeyPressFcn(obj, evd)
%processes shortcut keys

global XYW
nSel = get(XYW.UI.ListBoxCS,'Value');

%disp(lower(evd.Key))
switch lower(evd.Key)
    case 'return'
        %enter key pressed
        ButtonToggleROI_Callback;
    case 'quote'
        %new roi
        ButtonNewROI_Callback;
    case 'space'
        %toggle ROI visibility
        try
            XYW.showROI = XYW.showROI + 1;
            if XYW.showROI > 3
                XYW.showROI = 1;
            end
        catch
            XYW.showROI = 1; % show all
            %             2    show only good
            %             3    hide all
        end
        for j=1:length(XYW.CS(nSel).DAT.Blobs)
            if XYW.showROI==1 || (XYW.showROI==2 && XYW.CS(nSel).DAT.Blobs(j).bGood)
                set(XYW.CS(nSel).DAT.Blobs(j).h([2 3]), 'visible', 'on');
            else
                set(XYW.CS(nSel).DAT.Blobs(j).h([2 3]), 'visible', 'off');
            end
        end
    case {'r' 'g' 'b'}
        switch lower(evd.Key)
            case 'r'
                XYW.bRGB(1) = ~XYW.bRGB(1);
            case 'g'
                XYW.bRGB(2) = ~XYW.bRGB(2);
            case 'b'
                XYW.bRGB(3) = ~XYW.bRGB(3);
        end
        
        TMP = XYW.RGB;
        for k=1:3
            if ~XYW.bRGB(k)
                TMP(:,:,k) = 0;
            end
        end
        set(XYW.hRGB, 'CData', TMP);
        set(XYW.hRGBz,'CData', TMP );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = CallbackString(RootFunction, varargin)
%function str = CallbackString(RootFunction, varargin)
%
%generates a string to use as a callback funciton
%the arguments can be function names, strings, or numbers
%NOTE: function names are in single quotes, strings in triple quotes, or in cell-brackets
%Example: CallbackString(MyFunction, 'gcbo', '''hello''', {'bracketed'}, 3) will return (as a string:)
%        MyFunction(gcbo,'hello','bracketed',3)

str = strcat(RootFunction, '(');
for k=1:nargin-1
    if iscell(varargin{k})
        str = strcat(str, '''', char(varargin{k}), '''');    %we must add the quotes around this string
    elseif ischar(varargin{k})
        str = strcat(str, varargin{k});
    else
        %this is numerical, or matrix input -- must convert it into a string.
        %this requires special attention if there are multiple rows
        A = varargin{k};
        if isempty(A)
            str = strcat(str,'[]');
        else
            str = strcat(str, '[');  %opening bracket
            str = strcat(str, num2str(A(1,:)));  %first row
            for nRow=2:size(A,1)
                str = strcat(str, ';', num2str(A(nRow,:)));  %subsequent rows
            end
            str = strcat(str, ']');  %closing bracket
        end
    end
    if k<nargin-1
        str = strcat(str, ',');  %comma
    end
end
str = strcat(str, ')');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DoThatThing(RootName, MODE)
global XYW

disp(['RootName = ' RootName]);
disp(['MODE = ' MODE]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$$ KEY PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bGrnDeltaRiseOnly = 1; %default
if strcmpi(MODE,'XYW')
    %XYW
    NumRepsPerDose = 6;    %this is equal to RFPReps for XYW, but not for GH
    nDoseForSegmentation = 1;
    nRepForSegmentation = 6; %1 to NumRepsPerDose
    disp(['Segmenting from DOSE#' num2str(nDoseForSegmentation) ', REP#' num2str(nRepForSegmentation)]);
    
    %Critical settings for Process detection/elimination
    bBlockOutPROCESS    = sqrt(0.75);  %degree of suppression; 1==full suppression  
    IterLoopPROCESS     = [                      3];  %1 = Red;  2=delGrn;  3=Red+delGrn
    MinMaxGamHolesPROCESS = {[100 10000 0.5  0] [10  10000  0.5  0]}; %for Red then Grn
    GausPrePostPROCESS     = [0.5  1.5]; %(PreFilt - PostFilt); same for both Red/Grn
    bAbsDiffPROCESS     =  0; %binary flag; 0=positive rectify; 1=abs(x)
    FinalGausPROCESS   =  1; %final smoothing; same for Red/Grn
    ThreshPROCESS    = [   0.1         0.01    0.015];  %for red then Grn, then Red+Grn
    MinAreaPROCESS      = 0; %minimum area for process
    DilateRadiusPROCESS = 0; %dilate the ones we keep
    
    %Critical settings for Soma detection/elimination
    IterLoopSOMA     = [2];  %[    3        2 ];  %3 = Red*delGrn;  2=delGrn;  [1 2] = segment with both.
    MinMaxGamSOMA = {[100 10000 0.5]  [10  10000  0.5]}; %for Red then Grn
    ThreshSOMA    = [   0.1                0.08    ];  %for red then Grn, then Red*Grn
    GausPrePostSOMA  = [5  17]; %MUST BE ODD; (PreFilt - PostFilt); same for both Red/Grn
    bFillHolesSOMA    = 1;  %binary
    ErodeRadiusSOMA   = 1;
    DilateRadiusSOMA  = 0;
    
    %Final Min/Max area
    MinMaxArea  = [100  3000];
    MinDelG     = 100;
    EccenMax = 1;
    SolidMin = 0;
    

elseif strcmpi(MODE,'GH0') || strcmpi(MODE,'GH1') || strcmpi(MODE,'GH2') || strcmpi(MODE,'GH3') || strcmpi(MODE,'GH4') || strcmpi(MODE,'GH5') || strcmpi(MODE,'GH5.2')
   if strcmpi(MODE,'GH4')
        NumRepsPerDose = 4;    
    else
        NumRepsPerDose = 6;    %this is equal to RFPReps for XYW, but not for GH
    end
    if strcmpi(MODE,'GH0') || strcmpi(MODE,'GH3') || strcmpi(MODE,'GH5')
        %GH0 or 3 
        nDoseForSegmentation =  1; %first dose; %NO GABA on first dose
    elseif strcmpi(MODE,'GH5.2')
        nDoseForSegmentation =  2; %first dose; %NO GABA on first dose
    else
        %GH1 or 2 or 4
        nDoseForSegmentation = -1; %last dose; %has GABA on first dose
    end
    if strcmpi(MODE,'GH5') || strcmpi(MODE,'GH5.2')
        nRepForSegmentation = 1; %furst rep on this dose
    else
        nRepForSegmentation = NumRepsPerDose; %last rep on this dose
    end
    disp(['Segmenting from DOSE#' num2str(nDoseForSegmentation) ', REP#' num2str(nRepForSegmentation)]);
    
    %Critical settings for Process detection/elimination
    bBlockOutPROCESS    = sqrt(0.75);  %degree of suppression; 1==full suppression  
    IterLoopPROCESS     = [  3];  %1 = Red;  2=delGrn;  
    MinMaxGamHolesPROCESS = {[100 10000 0.5  0] [10  3000  0.5  0]}; %for Red then Grn
    GausPrePostPROCESS     = [0.5  1.5]; %(PreFilt - PostFilt); same for both Red/Grn
    bAbsDiffPROCESS     =  0; %binary flag; 0=positive rectify; 1=abs(x)
    FinalGausPROCESS   =  1; %final smoothing; same for Red/Grn
    ThreshPROCESS    = [   0.1         0.02      .02  ];  %for red then Grn
    MinAreaPROCESS      = 0; %minimum area for process
    DilateRadiusPROCESS = 0; %dilate the ones we keep
    
    %Critical settings for Soma detection/elimination
    IterLoopSOMA     = [2];  %[    3        2 ];  %3 = Red*delGrn;  2=delGrn;  [1 2] = segment with both.
    MinMaxGamSOMA = {[100 10000 0.5]  [10  3000  0.5]}; %for Red then Grn
    ThreshSOMA    = [   0.1     0.1   ];  %for red then Grn, then Red*Grn
    GausPrePostSOMA  = [5  17]; %MUST BE ODD; (PreFilt - PostFilt); same for both Red/Grn
    bFillHolesSOMA    = 1;  %binary
    ErodeRadiusSOMA   = 1;
    DilateRadiusSOMA  = 0;
    
    %Final Min/Max area
    MinMaxArea  = [100  3000];
    MinDelG     = 100;
    EccenMax = 1;
    SolidMin = 0;
    
elseif strcmpi(MODE,'EF1') ||  strcmpi(MODE,'EF2') ||  strcmpi(MODE,'EF2_WTF')
    %EF1
    NumRepsPerDose = 1;    %this is equal to RFPReps for XYW, but not for GH
    nDoseForSegmentation = 1;
    nRepForSegmentation = 1; %1 to NumRepsPerDose
    disp(['Segmenting from DOSE#' num2str(nDoseForSegmentation) ', REP#' num2str(nRepForSegmentation)]);
    
    %Critical settings for Process detection/elimination
    bBlockOutPROCESS    = sqrt(0.75);  %degree of suppression; 1==full suppression  
    IterLoopPROCESS     =  [ 3];  %1 = Red;  2=delGrn;  
    MinMaxGamHolesPROCESS = {[100 10000 0.5  0] [100  30000  0.5  0]}; %for Red then Grn
    GausPrePostPROCESS     = [0.5  1.5]; %(PreFilt - PostFilt); same for both Red/Grn
    bAbsDiffPROCESS     =  0; %binary flag; 0=positive rectify; 1=abs(x)
    FinalGausPROCESS   =  1; %final smoothing; same for Red/Grn
    ThreshPROCESS    = [   0.1         0.01    0.015];  %for red then Grn, then Red+Grn
    MinAreaPROCESS      = 0; %minimum area for process
    DilateRadiusPROCESS = 0; %dilate the ones we keep
    
    %Critical settings for Soma detection/elimination
    IterLoopSOMA     = [2];  %[    3        2 ];  %3 = Red*delGrn;  2=delGrn;  [1 2] = segment with both.
    MinMaxGamSOMA = {[100 10000 0.5]  [100  70000  0.5]}; %for Red then Grn
    ThreshSOMA    = [   0.1             0.1     0.03];  %for red then Grn, then Red*Grn
    GausPrePostSOMA  = [5  17]; %MUST BE ODD; (PreFilt - PostFilt); same for both Red/Grn
    bFillHolesSOMA    = 1;  %binary
    ErodeRadiusSOMA   = 1;
    DilateRadiusSOMA  = 0;
    
    %Final Min/Max area
    MinMaxArea  = [100  3000];
    MinDelG     = 100;
    EccenMax = 1;
    SolidMin = 0;
    
else
    error('Unknown MODE!');    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD RAW IMAGE DATA; only what's needed for segmentation!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WB = waitbar(0,'');

try
    NR = XYW.NR; %save time if we already know which data to load for segmenting
    NG = XYW.NG; %save time if we already know which data to load for segmenting
    RFPReps = XYW.RFPReps;
    GFPReps = XYW.GFPReps;
catch
    NR = []; %need to figure out NR
    NG = []; %need to figure out NG
    RFPReps = [];
    GFPReps = [];
end

if strcmpi(MODE,'EF2_WTF')
    %CODE ADDED Feb 10, 2021.  File mixup error  (WTF)
    [~,~,XL] = xlsread([XYW.PWD '\' RootName 'WTF.xlsx']);  %contains the correct file/frame lookup table to unmix the WTF error
    [VSIstr, ~,VSIn] = unique(XL(:,[1 3]));  %list of vsi files needed
    VSIframe = vertcat(XL{:,[2 4]});  %frame in each file
    
    %load the unique VSI files
    for k=1:length(VSIstr)
        waitbar(k/length(VSIstr), WB, 'Loading Raw Data ...');
        disp(VSIstr{k});
        VSI(k).tmp = bfopenIX83(VSIstr{k}); %load the next data
        fprintf(newline);
    end
    
    %map all the frames
    RED = [];
    GRN = [];
    for k=1:(length(VSIframe)/2)
        j = k+length(VSIframe)/2; %green indices are the second half
        RED(:,:,k) = VSI(VSIn(k)).tmp(:,:,VSIframe(k));
        GRN(:,:,k) = VSI(VSIn(j)).tmp(:,:,VSIframe(j));
    end
    
    for k=1:(length(VSIframe)/2)
    end
    
    clear VSI;  %clean up
    
    
    if isempty(NR)
        XYW.NR = 1:size(RED,3); %save for future
        XYW.NG = 1:size(GRN,3); %save for future
        XYW.RFPReps = size(RED,3); %save for future
        XYW.GFPReps = size(GRN,3); %save for future
    end
    NR = XYW.NR;
    NG = XYW.NG;
    RFPReps = XYW.RFPReps;
    GFPReps = XYW.GFPReps;
    NumDose = 1;
    NumRuns = 1;
    nTR = RFPReps * NumRuns; %total # timepoints
    nTG = GFPReps * NumRuns; %total # timepoints
 
    bNewTRITC = 0; %don't bother with this.
else
    
    D{1} = dir([XYW.PWD '\' RootName '*GFP*.vsi']);   %green data
    D{2} = dir([XYW.PWD '\' RootName '*TRITC*.vsi']); %red data
    D{4} = dir([XYW.PWD '\' RootName '*DAPI*.vsi']); %blue data (if it exists)
    
    if length(D{2})==0 && length(D{4})==2*length(D{1})
        disp('OMG -- I can not believe this -- DAPI is totally not TRITC, except for today');
        SWAPVAR = D{2};
        D{2} = D{4};
        D{4} = SWAPVAR;
    end
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
    if (strcmpi(MODE,'EF1') || strcmpi(MODE,'EF2')) &&  (length(D{2}) == length(D{1}) + 1)
        disp('Aha, you used a 2x mTRITC EF protocol. Nice!');
        D{3} = D{2}(1);     %short exposure time (avoid saturation)
        D{2} = D{2}(2:end); %original exposure time; for alignment
        bNewTRITC = 1; %FOR NOW; this isn't supported very well for EF
    elseif length(D{2}) == 2*length(D{1})
        disp('Aha, you used the new 2 x mTRITC protocol. Nice!');
        %9/2018 -- two mTRITC exposures (5ms=new and 20ms=old)
        D{3} = D{2}(1:2:end);
        D{2} = D{2}(2:2:end); %original exposure time
        bNewTRITC = 1;
    else
        bNewTRITC = 0;
    end
    if length(D{1}) == length(D{4})
        disp('Aha aha!  You too used blue - hi five!');
        bNewBLUE = 1;
    else
        bNewBLUE = 0;        
    end
    if length(D{1}) ~= length(D{2})
        error('mismatch in number of mGFP vs mTRITC files -- gut punch');
    end
    
    NumRuns = length(D{1}); %number of runs=number of files
    

    g=0; r=0;
    if ~isempty(findstr(MODE,'EF'))
        %EF1 or EF2 -- non-standard
        for k=1:length(D{1})
            Gtmp = bfopenIX83(D{1}(k).name); %load the next GFP data
            Rtmp = bfopenIX83(D{2}(k).name); %load the next mTRITC data
            if k==1 && bNewTRITC
                RED1 = bfopenIX83(D{3}(k).name); %load the next mTRITC data
            end
            if k==1 && bNewBLUE
                BLUE = bfopenIX83(D{4}(k).name); %load the next DAPI data
            end
            waitbar(k/length(D{1}), WB, 'Loading Raw Data ...');
            disp(D{1}(k).name);
            
            if ~isempty(GFPReps) && k==1
                %allocate memory, only if we know exactly how big it will be
                GRN = zeros(size(Gtmp,1),size(Gtmp,2),GFPReps, class(Gtmp)); %only allocate memory sufficient for all green data (note: at this stage, we will only populate NG timepoints)
                RED = zeros(size(Rtmp,1),size(Rtmp,2),RFPReps, class(Rtmp)); %only allocate memory sufficient for all red data (note: at this stage, we will only populate NR timepoints)
            end
            GFPReps = size(Gtmp,3); %number of repeitions for a single run
            RFPReps = size(Rtmp,3); %number of repeitions for a single run
            GRN(:,:,g+1:g+GFPReps) = Gtmp;  %copy
            RED(:,:,r+1:r+RFPReps) = Rtmp;  %copy
            r = r + RFPReps;
            g = g + GFPReps;
        end
        if strcmpi(MODE,'EF1')
            %EF1 needs to be re-sorted due to the nested loops; not needed for EF2
            GRN = GRN(:,:,[1:15 31 16:30 32]);  %last two images were actually interleaved
            RED = RED(:,:,[1:15 31 16:30 32]);  %last two images were actually interleaved
        end
        if isempty(NR)
            XYW.NR = 1:size(RED,3); %save for future
            XYW.NG = 1:size(GRN,3); %save for future
            XYW.RFPReps = size(RED,3); %save for future
            XYW.GFPReps = size(GRN,3); %save for future
        end
        NR = XYW.NR;
        NG = XYW.NG;
        RFPReps = XYW.RFPReps;
        GFPReps = XYW.GFPReps;
        NumDose = 1;
        NumRuns = 1;
        nTR = RFPReps * NumRuns; %total # timepoints
        nTG = GFPReps * NumRuns; %total # timepoints
    else
        %XYW or GH -- the standard
        bFirstDataDone = 0;
        for k=1:length(D{1})
            waitbar(k/length(D{1}), WB, 'Loading Raw Data ...');
            
            %   need to figure out NR  ||  need to load this data because it's part of NR
            if (isempty(NR) && k==1) || ~isempty(intersect(r+1:r+RFPReps, NR))
                disp(D{1}(k).name);
                Gtmp = bfopenIX83(D{1}(k).name); %load the next GFP data
                Rtmp = bfopenIX83(D{2}(k).name); %load the next mTRITC data
                if bNewTRITC
                    Rtmp1 = bfopenIX83(D{3}(k).name); %load the next mTRITC data
                end
                if bNewBLUE
                    Btmp = bfopenIX83(D{4}(k).name); %load the next mTRITC data
                end
                bDataLoaded = 1;
            else
                bDataLoaded = 0;
            end
            if bDataLoaded==1 && bFirstDataDone==0
                %we may have never determined NR .. need to figure this out
                
                RFPReps = size(Rtmp,3); %number of repetitions for a single run
                nTR = RFPReps * NumRuns; %total # timepoints
                
                GFPReps = size(Gtmp,3); %number of repetitions for a single run
                nTG = GFPReps * NumRuns; %total # timepoints
                
                %calculate minimum data needed
                NumDose = RFPReps*NumRuns/NumRepsPerDose;
                if nDoseForSegmentation<0
                    %-1 -->  NumDose
                    %-2 -->  NumDose-1 ..
                    nDoseForSegmentation = NumDose +1 + nDoseForSegmentation;
                end
                NR =  1                  + ((nDoseForSegmentation-1)*NumRepsPerDose + (nRepForSegmentation-1));
                NG = [1:GFPReps/RFPReps] + ((nDoseForSegmentation-1)*NumRepsPerDose + (nRepForSegmentation-1))*(GFPReps/RFPReps);
                
                GRN = zeros(size(Gtmp,1),size(Gtmp,2),nTG, class(Gtmp)); %only allocate memory sufficient for all green data (note: at this stage, we will only populate NG timepoints)
                RED = zeros(size(Rtmp,1),size(Rtmp,2),nTR, class(Rtmp)); %only allocate memory sufficient for all red data (note: at this stage, we will only populate NR timepoints)
                if bNewTRITC
                    RED1 = RED;  %allocate memory
                end
                if bNewBLUE
                    BLUE = RED;  %allocate memory
                end
                
                XYW.NR = NR; %save for future
                XYW.NG = NG; %save for future
                XYW.RFPReps = RFPReps; %save for future
                XYW.GFPReps = GFPReps; %save for future
                bFirstDataDone = 1; %done for this CS
            end
            if bDataLoaded
                GRN(:,:,g+1:g+GFPReps) = Gtmp;  %copy
                RED(:,:,r+1:r+RFPReps) = Rtmp;  %copy
                if bNewTRITC
                    RED1(:,:,r+1:r+RFPReps) = Rtmp1;  %copy
                end
                if bNewBLUE
                    BLUE(:,:,r+1:r+RFPReps) = Btmp;  %copy
                end
            end
            r = r + RFPReps;
            g = g + GFPReps;
        end
    end
    clear Gtmp Rtmp Rtmp1 Btmp;  %clean up
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALIGNMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    load([XYW.PWD '\' RootName '_dRXY.mat']);
catch
    error('MUST RUN XYW_ALIGN_IX83 prior to using XYW_GUI');
end

%ALIGN/ROTATE RFP
fprintf('\n*Applying Alignment:');
fprintf('\nR');
RED = RED + 1; %make sure none of the valid red pixels are zero, so that after rotation we know that any zeros were lost to rotation (Iclip)
for k=1:size(RED,3)
    %rotated-out pixels will become zero
    RED(:,:,k) = AlignImageGUI_IX83('ShiftImage', RED(:,:,k), dRXY(k,1), dRXY(k,2), dRXY(k,3));
    if bNewTRITC && k<=size(RED1,3)
        RED1(:,:,k) = AlignImageGUI_IX83('ShiftImage', RED1(:,:,k), dRXY(k,1), dRXY(k,2), dRXY(k,3));
    end
    if bNewBLUE && k<=size(BLUE,3)
        BLUE(:,:,k) = AlignImageGUI_IX83('ShiftImage', BLUE(:,:,k), dRXY(k,1), dRXY(k,2), dRXY(k,3));
    end
    fprintf([num2str(k) '.']);
    if bNewBLUE
        waitbar(k/size(RED,3), WB, 'Aligning RED & BLUE ...');
    else
        waitbar(k/size(RED,3), WB, 'Aligning RED ...');
    end
end
Iclip = find(min(RED,[],3)==0); %excludes pixels lost to rotation/alignment
RED = RED-1; %undo the +1 from before

%ALIGN/ROTATE GFP
fprintf('\nG');
for k=NG %1:size(GRN,3)
    j = floor((k-1)*RFPReps/GFPReps) + 1;
    GRN(:,:,k) = AlignImageGUI_IX83('ShiftImage', GRN(:,:,k), dRXY(j,1), dRXY(j,2), dRXY(j,3));
    if mod(k-1, GFPReps/RFPReps)==0
        fprintf([num2str(j) '.']);
    end
    waitbar(k/max(NG), WB, 'Aligning GRN ...');
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AUTO-FIND the Background --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Image.RED_NR   = RED(:,:,NR);
Image.RedMean  = mean(RED(:,:,NR),3);
if bNewTRITC
    if ~isempty(findstr(MODE,'EF'))
        %EF mode is special; the 5ms exposure TRITC only happens once at
        %the beginning
        Image.RED1_NR   = RED1;  
        Image.Red1Mean  = mean(RED1,3);
    else
        %normal mode -- both tritc exposures happen every iteration
        Image.RED1_NR   = RED1(:,:,NR);
        Image.Red1Mean  = mean(RED1(:,:,NR),3);
    end
end
if bNewBLUE
    %normal mode -- exposures happen every iteration
    Image.BLUE_NR   = BLUE(:,:,NR);
    Image.BlueMean  = mean(BLUE(:,:,NR),3);
end
Image.GRN_NG   = GRN(:,:,NG);  %allows for post-editing of ROIs
Image.GrnMean  = mean(GRN(:,:,NG),3);
if bGrnDeltaRiseOnly
    Image.GrnDelta = max(GRN(:,:,NG),[],3)-GRN(:,:,NG(1));
else
    Image.GrnDelta = max(GRN(:,:,NG),[],3)-min(GRN(:,:,NG),[],3);
end

%code to estimate background
if bNewBLUE
    RPG = double(Image.GrnDelta) + Image.GrnMean + Image.RedMean + Image.BlueMean; %red plus green plus bue over all data
else
    RPG = double(Image.GrnDelta) + Image.GrnMean + Image.RedMean; %red plus green over all data
end
RPG(Iclip) = NaN; %don't use the clipped pixels
RPG(1:round(0.25*size(RPG,1)),  :) = inf; %only use center of image
RPG(round(0.75*size(RPG,1)):end,:) = inf; %only use center of image
RPG(:,1:round(0.25*size(RPG,2)))   = inf; %only use center of image
RPG(:,round(0.75*size(RPG,2)):end) = inf; %only use center of image
[~,Ibkg] = sort(RPG(:));
Ibkg = Ibkg(1:1000); %lowest 1000 pixels


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Segmentation code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        RED                  deltaGRN
if bNewBLUE
    SegmentImage = cat(3, double(Image.RedMean), double(Image.GrnDelta), double(Image.BlueMean)); %two-color channels, Red then deltaGrn
    Iclip2 = [Iclip; Iclip + size(SegmentImage,1)*size(SegmentImage,2); Iclip + 2*size(SegmentImage,1)*size(SegmentImage,2)]; %indices that are clipped for above (covers indices of red and green two-color array above)
    SegmentImage(Iclip2) = 0;
else
    SegmentImage = cat(3, double(Image.RedMean), double(Image.GrnDelta)); %two-color channels, Red then deltaGrn
    Iclip2 = [Iclip; Iclip + size(SegmentImage,1)*size(SegmentImage,2)]; %indices that are clipped for above (covers indices of red and green two-color array above)
    SegmentImage(Iclip2) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIRST FIND PROCESSES TO BLOCK THEM OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProcessImBinary = 0*SegmentImage(:,:,1);
if bBlockOutPROCESS>0
    for ITER=IterLoopPROCESS
        %get the appropriate color channel apply Min/Max Clip Gamma
        if ITER==3
            %hack -- Red+Grn!!
            ThreshImTMP = zeros(size(SegmentImage(:,:,1)));
            for ITR=1:2
                Min = MinMaxGamSOMA{ITR}(1);
                Max = MinMaxGamSOMA{ITR}(2);
                Gam = MinMaxGamSOMA{ITR}(3);
                ThreshImTMPTMP = (SegmentImage(:,:,ITR)-Min)/(Max-Min);
                ThreshImTMPTMP(ThreshImTMPTMP<0)=0;
                ThreshImTMPTMP(ThreshImTMPTMP>1)=1;
                ThreshImTMPTMP=ThreshImTMPTMP.^Gam;
                ThreshImTMP= ThreshImTMP + ThreshImTMPTMP;
            end
            PreFillHolesPROCESS = MinMaxGamHolesPROCESS{ITR}(4);
        else
            %apply Min/Max Clip Gamma
            Min = MinMaxGamHolesPROCESS{ITER}(1);
            Max = MinMaxGamHolesPROCESS{ITER}(2);
            Gam = MinMaxGamHolesPROCESS{ITER}(3);
            ThreshImTMP = (SegmentImage(:,:,ITER)-Min)/(Max-Min);
            ThreshImTMP(ThreshImTMP<0)=0;
            ThreshImTMP(ThreshImTMP>1)=1;
            ThreshImTMP=ThreshImTMP.^Gam;
            PreFillHolesPROCESS = MinMaxGamHolesPROCESS{ITER}(4);
        end
    
        if PreFillHolesPROCESS>0
            ThreshImTMP = im2bw(ThreshImTMP, PreFillHolesPROCESS);
            ThreshImTMP = imfill(ThreshImTMP,'holes');
            ThreshImTMP = double(ThreshImTMP);
        end
        
        %take difference of narrow vs broad gaussian
        ThreshImDIFF = imgaussfilt(ThreshImTMP,GausPrePostPROCESS(1)) - imgaussfilt(ThreshImTMP,GausPrePostPROCESS(2));
        if bAbsDiffPROCESS==1
            ThreshImDIFF = abs(ThreshImDIFF); %accentuate edges
        else
            ThreshImDIFF(ThreshImDIFF<0) = 0; %positive-value rectify
        end
        if FinalGausPROCESS>0
            ThreshImDIFF = imgaussfilt(ThreshImDIFF,FinalGausPROCESS);
        end
        
        %threshold
        Thresh = ThreshPROCESS(ITER);
        ThreshImBinary = im2bw(ThreshImDIFF, Thresh);
        fprintf('Process* ');
        waitbar(find(ITER==IterLoopPROCESS)/length(IterLoopPROCESS), WB, 'Process Thresh ...');
        
        % Eliminate small specks or thin lines
        if MinAreaPROCESS>0
            ThreshImBinary = bwareaopen(ThreshImBinary,MinAreaPROCESS); %removes areas less than MinAreaPROCESS
            fprintf('Despeck* ');
            waitbar(find(ITER==IterLoopPROCESS)/length(IterLoopPROCESS), WB, 'Process Despeck ...');
        end
        ProcessImBinary(ThreshImBinary==1) = 1;
    end
    if DilateRadiusPROCESS>0
        ProcessImBinary = imdilate(ProcessImBinary,strel('disk',DilateRadiusPROCESS));
        fprintf('Dilate* ');
        waitbar(find(ITER==IterLoopPROCESS)/length(IterLoopPROCESS), WB, 'Process Dilate ...');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SECOND FIND SOMAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Blobs = [];
ThreshImBinary(:) = 0;
for ITER=IterLoopSOMA
    ThreshImBinaryPrev = ThreshImBinary;
    
    %get the appropriate color channel apply Min/Max Clip Gamma
    if ITER==3
        %hack -- Red*Grn!!
        ThreshImTMP = ones(size(SegmentImage(:,:,1)));
        for ITR=1:2
            Min = MinMaxGamSOMA{ITR}(1);
            Max = MinMaxGamSOMA{ITR}(2);
            Gam = MinMaxGamSOMA{ITR}(3);
            ThreshImTMPTMP = (SegmentImage(:,:,ITR)-Min)/(Max-Min);
            ThreshImTMPTMP(ThreshImTMPTMP<0)=0;
            ThreshImTMPTMP(ThreshImTMPTMP>1)=1;
            ThreshImTMPTMP=ThreshImTMPTMP.^Gam;
            ThreshImTMP= ThreshImTMP .* ThreshImTMPTMP;
        end
    else
        %apply Min/Max Clip Gamma
        Min = MinMaxGamSOMA{ITER}(1);
        Max = MinMaxGamSOMA{ITER}(2);
        Gam = MinMaxGamSOMA{ITER}(3);
        ThreshImTMP = (SegmentImage(:,:,ITER)-Min)/(Max-Min);
        ThreshImTMP(ThreshImTMP<0)=0;
        ThreshImTMP(ThreshImTMP>1)=1;
        ThreshImTMP=ThreshImTMP.^Gam;
    end
    waitbar(find(ITER==IterLoopSOMA)/length(IterLoopSOMA), WB, 'Soma Gamma ...');
    
    %exlude Processes (keep only somas)
    ThreshImTMP(ProcessImBinary==1) = (1-bBlockOutPROCESS)*ThreshImTMP(ProcessImBinary==1);

    %exlude previously chosen Somas
    ThreshImTMP(ThreshImBinaryPrev==1) = 0;
    
    %take difference of narrow vs broad gaussian
    ThreshImDIFF = imgaussfilt(ThreshImTMP,GausPrePostSOMA(1)) - imgaussfilt(ThreshImTMP,GausPrePostSOMA(2));
    
    %threshold
    Thresh = ThreshSOMA(ITER);
    ThreshImBinary = im2bw(ThreshImDIFF, Thresh);
    fprintf('SomaThresh* ');
    waitbar(find(ITER==IterLoopSOMA)/length(IterLoopSOMA), WB, 'Soma Thresh ...');
    
    %be sure to fully exlude previous processes and somas (note that due to gaussian blurring this needs to be done again.
    ThreshImBinary(ProcessImBinary==1) = (1-bBlockOutPROCESS)*ThreshImBinary(ProcessImBinary==1);
    ThreshImBinary(ThreshImBinaryPrev==1) = 0;
    
    %fill holes -- a necessary step given the PROCESS removal
    if bFillHolesSOMA
        ThreshImBinary = imfill(ThreshImBinary,'holes');
        waitbar(find(ITER==IterLoopSOMA)/length(IterLoopSOMA), WB, 'Soma Holes ...');
    end
    
    % Eliminate small specks or thin lines
    if ErodeRadiusSOMA>0
        ThreshImBinary = imerode(ThreshImBinary,strel('disk',ErodeRadiusSOMA)); %opens image
        fprintf('Erode* ');
        waitbar(find(ITER==IterLoopSOMA)/length(IterLoopSOMA), WB, 'Soma Erode ...');
    end
    if DilateRadiusSOMA>0
        ThreshImBinary = imdilate(ThreshImBinary,strel('disk',DilateRadiusSOMA)); %dilates image
        fprintf('Dilate* ');
        waitbar(find(ITER==IterLoopSOMA)/length(IterLoopSOMA), WB, 'Soma Dilate ...');
    end
    
    % Get all the blob properties
    BlobsTMP  = regionprops(ThreshImBinary, {'PixelIdxList' 'FilledArea' 'Centroid' 'Eccentricity' 'Solidity'});
    fprintf('Blobs* ');
    waitbar(find(ITER==IterLoopSOMA)/length(IterLoopSOMA), WB, 'Soma Blobs ...');
    
    %extract the outer boundary for display purposes
    %NOTE: MUST DO THIS BEFORE AREA!!
    B = bwboundaries(ThreshImBinary);
    for k=1:length(BlobsTMP)
        BlobsTMP(k).Boundary       = B{k}(:,end:-1:1);
    end
    fprintf('Bound* ');
    waitbar(find(ITER==IterLoopSOMA)/length(IterLoopSOMA), WB, 'Soma Bound ...');
    
    if ~isempty(BlobsTMP)
        Blobs = [Blobs; BlobsTMP];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIRD, ANALYZE BLOBS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Blobs)
    
    %eliminate small blobs that are mostly within larger blobs
    Area = [Blobs.FilledArea];
    [~,I] = sort(Area, 'descend');  %from largest to smallest
    Blobs = Blobs(I);
    for k = length(Blobs):-1:2 %from smallest to second largest
        for j=1:k-1  %everything larger than k
            if length(intersect(Blobs(j).PixelIdxList, Blobs(k).PixelIdxList)) >= 0.8*length(Blobs(k).PixelIdxList)
                %most (>=80%) of the pixels in the smaller blob are contained in the larger blob
                fprintf('Redundant Blob*');
                Blobs(k) = []; %eliminate the smaller blob
                waitbar(k/length(Blobs), WB, 'Blobs Redundant ...');
            end
        end
    end
    
    %calculate intensities from all color channels
    Blobs = AnalyzeBlobs(Blobs, Image, GRN, Ibkg, WB, NumRepsPerDose, NumDose);
    fprintf('Mean* ');
    fprintf('\n');

    %eliminate blobs that are too small or too big
    delG = [Blobs.delG];
    Area = [Blobs.FilledArea];
    Eccen = [Blobs.Eccentricity];
    Solid = [Blobs.Solidity];
    Blobs = Blobs(delG>=MinDelG & Area>=MinMaxArea(1) &  Area<=MinMaxArea(2) & Eccen<=EccenMax  & Solid>=SolidMin);
    fprintf('MinMaxArea MinDelG* ');
    waitbar(0, WB, 'Blobs Area/delG ...');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORGANIZE ALL SETTINGS AND DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AcqSettings.NumRuns = NumRuns;
AcqSettings.RFPReps = RFPReps;
AcqSettings.GFPReps = GFPReps;
AcqSettings.NumDose = NumDose;
AcqSettings.nTR = nTR;
AcqSettings.nTG = nTG;
AcqSettings.NumRepsPerDose = NumRepsPerDose;

SegSettings.MODE = MODE;  %XYW vs GH
SegSettings.nDoseForSegmentation = nDoseForSegmentation;
SegSettings.nRepForSegmentation = nRepForSegmentation;
SegSettings.bBlockOutPROCESS = bBlockOutPROCESS;
SegSettings.IterLoopPROCESS = IterLoopPROCESS; %1 = Red;  2=delGrn;  [1 2] = segment with both.
SegSettings.MinMaxGamHolesPROCESS = MinMaxGamHolesPROCESS;
SegSettings.ThreshPROCESS = ThreshPROCESS;
SegSettings.GausPrePostPROCESS = GausPrePostPROCESS;
SegSettings.bAbsDiffPROCESS = bAbsDiffPROCESS;
SegSettings.FinalGausPROCESS = FinalGausPROCESS;
SegSettings.MinAreaPROCESS = MinAreaPROCESS;
SegSettings.DilateRadiusPROCESS = DilateRadiusPROCESS;
SegSettings.IterLoopSOMA = IterLoopSOMA; %1 = Red;  2=delGrn;  [1 2] = segment with both.
SegSettings.MinMaxGamSOMA = MinMaxGamSOMA;
SegSettings.ThreshSOMA = ThreshSOMA;
SegSettings.GausPrePostSOMA = GausPrePostSOMA;
SegSettings.bFillHolesSOMA = bFillHolesSOMA;
SegSettings.ErodeRadiusSOMA = ErodeRadiusSOMA;
SegSettings.DilateRadiusSOMA = DilateRadiusSOMA;
SegSettings.MinMaxArea = MinMaxArea;
SegSettings.MinDelG = MinDelG;
SegSettings.EccenMax = EccenMax;
SegSettings.SolidMin = SolidMin;

SegOutput.Ibkg = Ibkg;
SegOutput.Iclip = Iclip;
SegOutput.SegmentImage = SegmentImage;
SegOutput.ProcessImBinary = ProcessImBinary;

waitbar(1, WB, 'SAVE ...');
save([XYW.PWD '\' RootName '1AutoSegment.mat'], 'Image', 'AcqSettings', 'SegSettings', 'SegOutput', 'Blobs');
fprintf('Save* ');

close(WB);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Blobs = AnalyzeBlobs(Blobs, Image, GRN, Ibkg, WB, NumRepsPerDose, NumDose)
global XYW

SigBkg = zeros(1,size(GRN,3));
for k=1:size(GRN,3)
    if any(k==XYW.NG)
        SigBkg(k) = mean(GRN(Ibkg + (k-1)*size(GRN,1)*size(GRN,2)));
    else
        SigBkg(k) = nan;
    end
    if ~isempty(WB)
        waitbar(k/size(GRN,3), WB, 'Calculating Background ...');
    end
end
fprintf('Background* ');

delGImage = double(Image.GrnDelta);
for k=1:length(Blobs)
    Blobs(k).redPix = Image.RedMean(Blobs(k).PixelIdxList);
    Blobs(k).red = mean(Blobs(k).redPix);
    Blobs(k).redVar = var(Blobs(k).redPix);  %spatial measure of delG variance
    Blobs(k).redStd = std(Blobs(k).redPix);  %spatial measure of delG variance
    
    try
        Blobs(k).red1Pix = Image.Red1Mean(Blobs(k).PixelIdxList);
        Blobs(k).red1 = mean(Blobs(k).red1Pix);
        Blobs(k).red1Var = var(Blobs(k).red1Pix);  %spatial measure of delG variance
        Blobs(k).red1Std = std(Blobs(k).red1Pix);  %spatial measure of delG variance
    catch
    end
    try
        Blobs(k).bluePix = Image.BlueMean(Blobs(k).PixelIdxList);
        Blobs(k).blue = mean(Blobs(k).bluePix);
        Blobs(k).blueVar = var(Blobs(k).bluePix);  %spatial measure of delG variance
        Blobs(k).blueStd = std(Blobs(k).bluePix);  %spatial measure of delG variance
    catch
    end
    
    Blobs(k).delGPix = delGImage(Blobs(k).PixelIdxList);       %delG of each pixel
    Blobs(k).delG    = mean(Blobs(k).delGPix);  %delG averaged over space
    Blobs(k).delGvar = var(Blobs(k).delGPix);  %spatial measure of delG variance
    Blobs(k).delGstd = std(Blobs(k).delGPix);  %spatial measure of delG variance
    
    Blobs(k).grnRawWithBkg = zeros(1,size(GRN,3));
    for j=1:size(GRN,3)
        Blobs(k).grnRawWithBkg(j) = mean(GRN(Blobs(k).PixelIdxList + (j-1)*size(GRN,1)*size(GRN,2)));
    end
    Blobs(k).grnRaw = Blobs(k).grnRawWithBkg - SigBkg;  %perform background subtraction here; adds NaNs!!
    
    %Definition of Min/Max depends on the mode.r
    if strcmpi(XYW.MODE, 'GH0') || strcmpi(XYW.MODE, 'GH1') || strcmpi(XYW.MODE,'GH5') || strcmpi(XYW.MODE,'GH5.2')
        %each repetition has a different GFP power
        %reshape data to be         time(1:16)    x  illum(1:6)  x   Dose
        y=reshape(Blobs(k).grnRaw, [XYW.GFPReps/XYW.RFPReps, NumRepsPerDose, NumDose]);
        
        %take mean for each illum
        Blobs(k).yMax = max(max(y,[],1),[],3); %fcn of illum(1:6)
        Blobs(k).yMin = min(min(y,[],1),[],3); %fcn of illum(1:6)
        
        Blobs(k).yMaxNorm = Blobs(k).yMax/Blobs(k).yMax(end);
        Blobs(k).yMinNorm = Blobs(k).yMin/Blobs(k).yMin(end);
        
        %reepeat this data to be the same size as raw data
        Blobs(k).grnMax = repmat(Blobs(k).yMax, [XYW.GFPReps/XYW.RFPReps 1 NumDose]);
        Blobs(k).grnMin = repmat(Blobs(k).yMin, [XYW.GFPReps/XYW.RFPReps 1 NumDose]);
        Blobs(k).grnMax = Blobs(k).grnMax(:)';
        Blobs(k).grnMin = Blobs(k).grnMin(:)';
    else  %GH2 or GH3 or GH4 or XYW
        if k==1
            try
                Blobs = rmfield(Blobs, 'yMax');
                Blobs = rmfield(Blobs, 'yMin');
                Blobs = rmfield(Blobs, 'yMaxNorm');
                Blobs = rmfield(Blobs, 'yMinNorm');
            catch
            end
        end
        
        Blobs(k).grnMax = max(Blobs(k).grnRaw);  %scalar
        Blobs(k).grnMin = min(Blobs(k).grnRaw);  %scalar
    end
    Blobs(k).grnNorm = (Blobs(k).grnRaw - Blobs(k).grnMin)./(Blobs(k).grnMax - Blobs(k).grnMin);
    
    if ~isempty(WB)
        waitbar(k/length(Blobs), WB, 'Blobs Mean ...');
    end
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bAnlysisDone = FullAnalysis
global XYW
nSel = get(XYW.UI.ListBoxCS,'Value');
bAnlysisDone = 0; %default return

if isempty(XYW.CS(nSel).DAT) || isempty(XYW.CS(nSel).DAT.Blobs)
    disp(['CS#' num2str(nSel) ' Full Analysis Not Done; No Blobs']);
    return;
end

global HACK
bAnyNeeded = 0;
for k = 1:length(XYW.CS(nSel).DAT.Blobs)
    if isempty(HACK) && ~XYW.CS(nSel).DAT.Blobs(k).bViewed
        disp(['CS#' num2str(nSel) ' Full Analysis Not Allowed Until All Viewed']);
        return;
    end
    if ~bAnyNeeded && any(isnan(XYW.CS(nSel).DAT.Blobs(k).grnRaw))
        bAnyNeeded = 1;
    end
end
if ~bAnyNeeded
    %disp(['CS#' num2str(nSel) ' Full Analysis Previously Done']);
    %return;
end

if ~isfield(XYW.CS(nSel).DAT, 'Image')
    %load the auto-segmentation details
    ListBoxCS_Callback;
    drawnow;
end
%XYW.CS(nSel).DAT.AcqSettings.NumRuns;
%XYW.CS(nSel).DAT.AcqSettings.RFPReps;
%XYW.CS(nSel).DAT.AcqSettings.GFPReps;
%XYW.CS(nSel).DAT.AcqSettings.NumDose;
%XYW.CS(nSel).DAT.AcqSettings.nTR;
%XYW.CS(nSel).DAT.AcqSettings.nTG;
%XYW.CS(nSel).DAT.AcqSettings.NumRepsPerDose;

%XYW.CS(nSel).DAT.SegOutput.Ibkg;
%XYW.CS(nSel).DAT.SegOutput.Iclip;
%XYW.CS(nSel).DAT.SegOutput.SegmentImage;
%XYW.CS(nSel).DAT.SegOutput.ProcessImBinary;

try
    XYW.CS(nSel).bFullAnalysisDone;
catch
    XYW.CS(nSel).bFullAnalysisDone = false;
end

if XYW.CS(nSel).bFullAnalysisDone
    disp(['Analysis already done for CS#' num2str(nSel)]);
    return;
end

RootName = ['CS' num2str(XYW.CS(nSel).Num) '_'];

%GET DIRECTORY INFORMATION
D{1} = dir([XYW.PWD '\' RootName '*GFP*.vsi']);   %green data
D{2} = dir([XYW.PWD '\' RootName '*TRITC*.vsi']); %red data
D{4} = dir([XYW.PWD '\' RootName '*DAPI*.vsi']); %blue data
if length(D{2})==0 && length(D{4})==2*length(D{1})
    disp('OMG -- I can not believe this -- DAPI is totally not TRITC, except for today');
    SWAPVAR = D{2};
    D{2} = D{4};
    D{4} = SWAPVAR;
end
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
%Note: EF mode doesn't use this code at all
if length(D{2}) == 2*length(D{1})
    disp('Aha, you used the new 2 x mTRITC protocol. Nice!');
    %9/2018 -- two mTRITC exposures (5ms=new and 20ms=old)
    D{3} = D{2}(1:2:end);
    D{2} = D{2}(2:2:end); %original exposure time
    bNewTRITC = 1;
else
    bNewTRITC = 0;
end
if length(D{1}) == length(D{4})
    disp('Aha aha!  You too used blue - hi five!');
    bNewBLUE = 1;
else
    bNewBLUE = 0;
end
if length(D{1}) ~= length(D{2})
    error('mismatch in number of mGFP vs mTRITC files -- gut punch');
end
%CHECK
if length(D{1}) ~= XYW.CS(nSel).DAT.AcqSettings.NumRuns
    error('mismatch in number of mGFP files vs previously determined NumRuns');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD RAW IMAGE DATA; what's needed for the rest of the waveform analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GRN = zeros( ...
    size(XYW.CS(nSel).DAT.Image.GRN_NG,1), ...
    size(XYW.CS(nSel).DAT.Image.GRN_NG,2), ...
    XYW.CS(nSel).DAT.AcqSettings.nTG, ...
    class(XYW.CS(nSel).DAT.Image.GRN_NG)); %allocate memory sufficient for all green data
%RED = zeros( ...
%    size(XYW.CS(nSel).DAT.Image.RED_NR,1), ...
%    size(XYW.CS(nSel).DAT.Image.RED_NR,2), ...
%    XYW.CS(nSel).DAT.AcqSettings.nTR, ...
%    class(XYW.CS(nSel).DAT.Image.RED_NR)); %allocate memory sufficient for all red data
%if bNewTRITC
%    %short exposure red
%    RED1 = RED;
%end

NR = XYW.NR; %indices of data already loaded 
NG = XYW.NG; %indices of data already loaded 
GFPReps = XYW.CS(nSel).DAT.AcqSettings.GFPReps;
RFPReps = XYW.CS(nSel).DAT.AcqSettings.RFPReps;
NumRepsPerDose = XYW.CS(nSel).DAT.AcqSettings.NumRepsPerDose;
NumRuns        = XYW.CS(nSel).DAT.AcqSettings.NumRuns;
NumDose = RFPReps*NumRuns/NumRepsPerDose;

WB = waitbar(0,'');
g=0; 
%r=0;
for k=1:length(D{1})
    waitbar(k/length(D{1}), WB, 'Loading Raw Data ...');

    if 1 %isempty(intersect(g+1:g+GFPReps, NG))
        %####TEMPORARY SOLUTION--ALWAYS LOAD DATA FROM SCRATCH; NEED TO ACCOUNT FOR NG vs GFPReps discrepancy for future improvements (i.e., we need code that loads a subset of a .vsi file)
        disp(D{1}(k).name);
        GRN(:,:,g+1:g+GFPReps) = bfopenIX83(D{1}(k).name); %load the next GFP data
    else
        disp('NG');
        GRN(:,:,g+1:g+GFPReps) = XYW.CS(nSel).DAT.Image.GRN_NG; %this was previously loaded
    end
    g = g + GFPReps;
    
    %if 1 %isempty(intersect(r+1:r+RFPReps, NR))
    %    %####TEMPORARY SOLUTION--ALWAYS LOAD DATA FROM SCRATCH; NEED TO ACCOUNT FOR NG vs GFPReps discrepancy for future improvements (i.e., we need code that loads a subset of a .vsi file)
    %    disp(D{2}(k).name);
    %    RED(:,:,r+1:r+RFPReps) = bfopenIX83(D{2}(k).name); %load the next RFP data
    %else
    %    disp('NR');
    %    RED(:,:,r+1:r+RFPReps) = XYW.CS(nSel).DAT.Image.RED_NR; %this was previously loaded
    %end
    %if bNewTRITC
    %    disp(D{3}(k).name);
    %    RED1(:,:,r+1:r+RFPReps) = bfopenIX83(D{3}(k).name); %load the next RFP data
    %end
    %r = r + RFPReps;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALIGNMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    load([XYW.PWD '\' RootName '_dRXY.mat']);
catch
    error('MUST RUN XYW_ALIGN_IX83 prior to using XYW_GUI');
end

%ALIGN/ROTATE RFP
fprintf('\n*Applying Alignment:');
fprintf('\nR');

%for k=setdiff(1:size(RED,3), NR) %all but NR
%####TEMPORARY SOLUTION--ALWAYS LOAD DATA FROM SCRATCH; NEED TO ACCOUNT FOR NG vs GFPReps discrepancy for future improvements (i.e., we need code that loads a subset of a .vsi file)
%for k=1:size(RED,3)
%    RED(:,:,k) = AlignImageGUI_IX83('ShiftImage', RED(:,:,k), dRXY(k,1), dRXY(k,2), dRXY(k,3));
%    if bNewTRITC
%        RED1(:,:,k) = AlignImageGUI_IX83('ShiftImage', RED1(:,:,k), dRXY(k,1), dRXY(k,2), dRXY(k,3));
%    end
%    fprintf([num2str(k) '.']);
%    waitbar(k/size(RED,3), WB, 'Aligning RED ...');
%end

%ALIGN/ROTATE GFP
fprintf('\nG');
%for k=setdiff(1:size(GRN,3), NG) %all but NG
%####TEMPORARY SOLUTION--ALWAYS LOAD DATA FROM SCRATCH; NEED TO ACCOUNT FOR NG vs GFPReps discrepancy for future improvements (i.e., we need code that loads a subset of a .vsi file)
for k=1:size(GRN,3)
    j = floor((k-1)*RFPReps/GFPReps) + 1;
    GRN(:,:,k) = AlignImageGUI_IX83('ShiftImage', GRN(:,:,k), dRXY(j,1), dRXY(j,2), dRXY(j,3));
    if mod(k-1, GFPReps/RFPReps)==0
        fprintf([num2str(j) '.']);
    end
    waitbar(k/size(GRN,3), WB, 'Aligning GRN ...');
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AUTO-FIND the Background --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ibkg = XYW.CS(nSel).DAT.SegOutput.Ibkg;
SigBkg = zeros(1,size(GRN,3));
for k=1:size(GRN,3)
    SigBkg(k) = mean(GRN(Ibkg + (k-1)*size(GRN,1)*size(GRN,2)));
    waitbar(k/size(GRN,3), WB, 'Calculating Background ...');
end
fprintf('Background* ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYZE BLOBS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Blobs = XYW.CS(nSel).DAT.Blobs;
%calculate intensities from all color channels
for k = 1 : length(Blobs)           % Loop through all blobs.
    for j=1:size(GRN,3)
        Blobs(k).grnRawWithBkg(j) = mean(GRN(Blobs(k).PixelIdxList + (j-1)*size(GRN,1)*size(GRN,2)));
    end
    Blobs(k).grnRaw = Blobs(k).grnRawWithBkg - SigBkg;  %perform background subtraction here; adds NaNs!!
    
    %Definition of Min/Max depends on the mode.
    if strcmpi(XYW.MODE, 'GH0') || strcmpi(XYW.MODE, 'GH1') || strcmpi(XYW.MODE,'GH5') || strcmpi(XYW.MODE,'GH5.2')
        %each repetition has a different GFP power
        %reshape data to be         time(1:16)    x  illum(1:6)  x   Dose
        y=reshape(Blobs(k).grnRaw, [GFPReps/RFPReps, NumRepsPerDose, NumDose]);
        
        %take mean for each illum
        Blobs(k).yMax = max(max(y,[],1),[],3); %fcn of illum(1:6)
        Blobs(k).yMin = min(min(y,[],1),[],3); %fcn of illum(1:6)
        
        Blobs(k).yMaxNorm = Blobs(k).yMax/Blobs(k).yMax(end);
        Blobs(k).yMinNorm = Blobs(k).yMin/Blobs(k).yMin(end);
        
        %reepeat this data to be the same size as raw data
        Blobs(k).grnMax = repmat(Blobs(k).yMax, [GFPReps/RFPReps 1 NumDose]);
        Blobs(k).grnMin = repmat(Blobs(k).yMin, [GFPReps/RFPReps 1 NumDose]);
        Blobs(k).grnMax = Blobs(k).grnMax(:)';
        Blobs(k).grnMin = Blobs(k).grnMin(:)';
    else  %GH2 or 3
        if k==1
            try
                Blobs = rmfield(Blobs, 'yMax');
                Blobs = rmfield(Blobs, 'yMin');
                Blobs = rmfield(Blobs, 'yMaxNorm');
                Blobs = rmfield(Blobs, 'yMinNorm');
            catch
            end
        end
        
        Blobs(k).grnMax = max(Blobs(k).grnRaw);  %scalar
        Blobs(k).grnMin = min(Blobs(k).grnRaw);  %scalar
    end
    Blobs(k).grnNorm = (Blobs(k).grnRaw - Blobs(k).grnMin)./(Blobs(k).grnMax - Blobs(k).grnMin);
    
    waitbar(k/length(Blobs), WB, 'Blobs Mean ...');
end
fprintf('Mean* ');
fprintf('\n');

%SAVE to global variable:
XYW.CS(nSel).DAT.Blobs = Blobs;

close(WB);
bAnlysisDone = 1;