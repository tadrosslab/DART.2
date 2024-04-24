function PlotCSExemplars_PUBLICATION2(CELLS, Ipos, Ineg, bCollapse, REPSPERDOSE)
%function PlotCSExemplars_PUBLICATION2(CELLS, Ipos, Ineg, bCollapse, REPSPERDOSE)
%
%See code in YM90K_SCRIPT which uses this to create exemplars for the 2023 Nature Methods paper

clf; hold on;
try
    REPSPERDOSE;
catch
    REPSPERDOSE =  6;
end
FRAMEPERREP = 16;
t = [1:FRAMEPERREP];

yOffset=0;
for k=[Ipos(:); Ineg(:)]'
    if ~bCollapse || k==(Ineg(1)) || k==(Ipos(1))
        yOffset = yOffset-1.5;
    end
    
    if(ismember(k,Ineg))
        clr = 'k';
    elseif(ismember(k,Ipos))
        clr = 'r';
    else
        clr = [0.5 0.5 0.5];
    end
    
    
    plot([-10 180]*2,[yOffset yOffset],  'k');
    plot([-10 180]*2,[yOffset yOffset]+1,'k');
    
    
    j=1;
    xOffset = 0;
    
    y = CELLS(k).grnNorm;
    for d=1:length(CELLS(k).DOSE)
        if isnan(CELLS(k).RESPONSE(d))
            %this does is missing
        else
            tmp = t*0;
            for r=1:REPSPERDOSE
                yy = y(j:j+FRAMEPERREP-1);
                plot(t+xOffset, yy+yOffset,'color',[0 0 0]+0.5);
                tmp = tmp + yy + yOffset;
                j = j+FRAMEPERREP;
            end
            plot(t+xOffset, tmp/REPSPERDOSE, 'color', clr, 'linewidth', 3);
        end
        xOffset = xOffset + FRAMEPERREP + 10;
    end
end
axis([-20 180 -12 0]);
axis([-20 220 -12 0]);
axis([-20 320 -30 0]);
axis off;

return;

