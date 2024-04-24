
try
    CELLS;
catch
    load CELLS.mat;  %CELLS OPTIONS; gabzine DATA in 2022 paper
end
OPTIONS.bShowScatterPos = 0;
OPTIONS.bShowScatterNeg = 0;
XYW_DoseResponsePlot2(CELLS, OPTIONS, 5+[[1  [ 5   7 ]-1]]);
set(gca,'ydir','reverse', 'ylim', [-0.5 1.5]);
xticks('auto'); xticklabels('auto');

ambient = 10.^(-12:.01:-3); %ambient dose

clr = { [0 0 0]  [0 0 1]  [1 0 0]};
RX50  =   [13 70 180]*1e-6;  %gabazine series; experimental Rx50
%RxHTL50 = [20  50]*1e-9;      %gabazine series; pharmacological HTL50(t=15m)

A=1;  %tradeoff of capture vs trafficking (essentially sfHTP vs Promega's)
Trafficking = A;  
HTL50=A*80e-9;  %capture-only (not pharmacological) HTL50;  note: sfHTP = 30e-9 with Trafficking=1; (actual for this dataset)    promegaHTP = 10e-9 with Trafficking=0.5


ReceptorTotal = 250e-6;    %rEceptor total density, in units of “tethered total concentration”
if 0
   ReceptorTotal = ReceptorTotal/1e12;
end
MaxTetherRxTotal= 650e-6; %RX total density, in units of “tethered total concentration” assuming Trafficking=1

TetherFrac=Trafficking*(1-exp(-log(2)*ambient./HTL50)); %fraction of HTP bound to an Rx

%clf; semilogx(ambient, TetherFrac, '--k'); 
grid on; hold on; 

for k=1:length(RX50)
    %complete binding equation, with limited Rx and Receptor
    TetherRxTotal = MaxTetherRxTotal*TetherFrac;
    a = ReceptorTotal;
    b = -(ReceptorTotal+TetherRxTotal+RX50(k));
    c = TetherRxTotal;
    
    FracReceptorBoundToTetherRx=(-b-sqrt(b.^2 - 4*a*c))/(2*a);  %quadratic solution; fraction of native receptors bound to a tethered Rx
    FracReceptorBoundToAmbientRx = ambient./(ambient+RX50(k));

    h = semilogx([ambient NaN ambient]*1e6, 0.02+[FracReceptorBoundToAmbientRx NaN  FracReceptorBoundToTetherRx ], 'linewidth', 2, 'color', clr{k}); 
    %semilogx([RX50(k)  RxHTL50(k)], [0.5 0.5], 'o', 'color', get(h,'color'), 'linewidth', 2);
end
%ylim([0 1]);
set(gca,'ydir','reverse');
