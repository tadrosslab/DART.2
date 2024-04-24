
for Q=0:1
    figure(Q+1); clf;
    if Q==0
        %diazepam/flumazenil/DMSO assay controls
        nSettings = 8.3;
        DirVect = { ...
            'W:\JFRC_Data\ASSAYS_IN_NEURONS\flumazenil-DART\2017_05_25a' ...
            'W:\JFRC_Data\ASSAYS_IN_NEURONS\flumazenil-DART\2017_05_25b' ...
            'W:\JFRC_Data\ASSAYS_IN_NEURONS\flumazenil-DART\2017_05_25c' ...
            'W:\JFRC_Data\ASSAYS_IN_NEURONS\flumazenil-DART\2017_05_26' ...
            };
        
        %               diazepam control             regular flumazenil                DMSO or nothing
        PosVectVect = {{[2 7:10]  [1:6] [] [1:3]}   {[3:6 12] [7 8 10 11] [] [5:8]}  {[]  []  [1:12] [9:12]}};
        
        for j=1 %s1:length(PosVectVect)
            PositionVect = PosVectVect{j};
            
            if 0
                %overlay colors
                clr = MapValue2Color(j,[0 4],jet);
                DATA = ChrGCHT_DoseResponse(nSettings, 0, PositionVect, DirVect, 'o', clr, 0.5*clr, [0.5 0.5 0.5]);
            else
                %standard red/black
                DATA = ChrGCHT_DoseResponse(nSettings, 0, PositionVect, DirVect, 'o', 'r', 'k', [0.5 0.5 0.5]);
            end
            
            drawnow;
        end
    elseif Q==1
        %first dosing scheme
        nSettings = 8.1;
        DirVect = { ...
            'W:\JFRC_Data\ASSAYS_IN_NEURONS\flumazenil-DART\2017_05_18_fluDART_pilot2' ...
            };
        
        PositionVect = { ...
            [1:12] ...
            };
        DATA1 = ChrGCHT_DoseResponse(nSettings, 1, PositionVect, DirVect, 'o', 'r', 'k', [0.5 0.5 0.5]);
        
        
        pause;
        %second dosing scheme
        nSettings = 8.2;
        DirVect = { ...
            'W:\JFRC_Data\ASSAYS_IN_NEURONS\flumazenil-DART\2017_05_19_fluDART_pilot3' ...
            };
        
        PositionVect = { ...
            [1:12] ...
            };
        DATA2 = ChrGCHT_DoseResponse(nSettings, 1, PositionVect, DirVect, 'o', 'r', 'k', [0.5 0.5 0.5]);
        
        pause;
        
        OUT = ChrGCHT_DoseResponseMerge([DATA1 DATA2])
        OUT.STATS
        
    end
    subplot(1,6,1:3); 
    set(gca,'xlim',10.^[-8.5  6]);
    pause;
end

%PlotCSExemplars_PUBLICATION(DATA,[],[1]);
