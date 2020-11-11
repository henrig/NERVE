function Emission_Factors_OnRoadAllCond()
global tfold Tyear SSB_Vehicle_dist comps Vehicle_dist Vehicle_weight
global debug_mode
%--------------------------------------------------------------
% COMBINE Municipal VEHICLE DISTRIBUTION & EMISSION FACTORS
%--------------------------------------------------------------
TM = readtable(SSB_Vehicle_dist,'Sheet','MODEL');

LightVehiclesIdx = TM.ClassNum==1|TM.ClassNum==2;
BusesVehiclesIdx = TM.ClassNum==3|TM.ClassNum==4;
HeavyVehiclesIdx = TM.ClassNum==5|TM.ClassNum==6|TM.ClassNum==7;
VD               = Vehicle_dist.Vdist;

% Loop all components to be calculated:
for com = 1:length(comps)
    
    % Load the component file with
    fprintf('<--- %s\n',char(comps(com)))
    TFout = readtable('modelHBEFA.xlsx','Sheet',sprintf('%s_%s_Weight',char(comps(com)),Vehicle_weight),'PreserveVariableNames',1);
    
    for i= 1:height(TM)
        idef(i) = find(ismember(TFout.Properties.VariableNames,TM.Name(i)));
    end
    
    Trout    =  TFout(:,1:7);
    for komm = 1: length(Vehicle_dist.D1_KommNr)
        if debug_mode
            fprintf('Sum Weights: Komm: %04i@\n \tT:%f\n \tL:%f\n \tH:%f\n \tB:%f \n',...
                Vehicle_dist.D1_KommNr(komm),sum(VD(komm,:)),...
                sum(VD(komm,LightVehiclesIdx)),...
                sum(VD(komm,HeavyVehiclesIdx)),...
                sum(VD(komm,BusesVehiclesIdx)))
        end
        if abs(sum(VD(komm,:))-3)>1e-5
            VD(komm,LightVehiclesIdx) = VD(komm,LightVehiclesIdx)/sum(VD(komm,LightVehiclesIdx));
            VD(komm,HeavyVehiclesIdx) = VD(komm,HeavyVehiclesIdx)/sum(VD(komm,HeavyVehiclesIdx));
            VD(komm,BusesVehiclesIdx) = VD(komm,BusesVehiclesIdx)/sum(VD(komm,BusesVehiclesIdx));
            if debug_mode
                fprintf('\tNEW Sum Weights: Komm: %04i@\n \t\tT:%f NV:%i \n \t\tL:%f NV:%i \n \t\tH:%f NV:%i \n \t\tB:%f NV:%i \n',...
                    Vehicle_dist.D1_KommNr(komm),sum(VD(komm,:)),sum(Vehicle_dist.modelNV(komm,:)),...
                    sum(VD(komm,LightVehiclesIdx)),sum(Vehicle_dist.modelNV(komm,LightVehiclesIdx)),...
                    sum(VD(komm,HeavyVehiclesIdx)),sum(Vehicle_dist.modelNV(komm,HeavyVehiclesIdx)),...
                    sum(VD(komm,BusesVehiclesIdx)),sum(Vehicle_dist.modelNV(komm,BusesVehiclesIdx)))
            end
        end
        
        for i= 1:height(TFout)
            EFs = VD(komm,:).*table2array(TFout(i,idef));
            % rows2vars(TFout(1,idef))
            EFLight(i,1) = nansum(EFs(LightVehiclesIdx));
            EFHeavy(i,1) = nansum(EFs(HeavyVehiclesIdx));
            EFBuses(i,1) = nansum(EFs(BusesVehiclesIdx));
        end
        
        Trout.Light = EFLight;
        Trout.Heavy = EFHeavy;
        Trout.Buses = EFBuses;
        
        Trout.Properties.VariableNames(find(ismember(Trout.Properties.VariableNames,'Light'))) = {sprintf('EF_Light_%04i',Vehicle_dist.D1_KommNr(komm))};
        Trout.Properties.VariableNames(find(ismember(Trout.Properties.VariableNames,'Heavy'))) = {sprintf('EF_Heavy_%04i',Vehicle_dist.D1_KommNr(komm))};
        Trout.Properties.VariableNames(find(ismember(Trout.Properties.VariableNames,'Buses'))) = {sprintf('EF_Buses_%04i',Vehicle_dist.D1_KommNr(komm))};
    end
    writetable(Trout,'OnRoadEF_RoadClasses.xlsx','Sheet',sprintf('%s_%i',char(comps(com)),Tyear))
end

end