%--------------------------------------------------------------------------
% This file is part of NERVE
%
% NERVE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3.
%
% NERVE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with NERVE.  If not, see <https://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function Emission_Factors_OnRoadAllCond()
% Function for Road link calculations computes for each municipality an
% emission factor (EF) for Light, Heavy and Buses.
% This needs to be run for each component and each year.
global tfold Tyear SSB_Vehicle_dist comps Vehicle_dist Vehicle_weight
global debug_mode ofiles use_temporary_files input

fprintf('---------------------------------------------------------------\n')
fprintf('in Emission_Factors_OnRoadAllCond *\n')
fprintf('---------------------------------------------------------------\n')

end

TM = readtable(SSB_Vehicle_dist,'Sheet','MODEL');
LightVehiclesIdx = TM.ClassNum==1|TM.ClassNum==2;
BusesVehiclesIdx = TM.ClassNum==3|TM.ClassNum==4;
HeavyVehiclesIdx = TM.ClassNum==5|TM.ClassNum==6|TM.ClassNum==7;
VD               = Vehicle_dist.Vdist;
% Loop all components to be calculated:
for com = 1:length(comps)
    fprintf('<--- %s\n',char(comps(com)))
    ofile   = sprintf('OnRoadEF_RoadClasses_%s.xlsx',char(comps(com)));
    osheet  = sprintf('%s_%i',char(comps(com)),Tyear);
    tmpfile = sprintf('%sEF_On_AllRoadCond_Municipality_%i_%s.mat',tfold,Tyear,char(comps(com)));
    
    if use_temporary_files && exist(tmpfile)
        fprintf('Has all necessary  EF files for :%s \n',char(comps(com)))
    else
        
        oEFfile = sprintf('%sEFA_Table_MODEL_%s.mat',tfold,char(comps(com)));
        load(oEFfile);
        fprintf('Found Temporary HBEFA EF on file, \n------------------------\n')
        if input.options.Subtract_bio_from_CO2
            switch char(comps(com))
                case 'CO2'
                    fprintf('### CO2 component: Factor in Biofuels \n')
                    TFout = Emissions_Factor_Mix_in_biofuels(TFout);
            end
        end
        
        
        % Calculate an EF file that can be used on roads. Light / Heavy / Buses
        % for all road types and all municpalities.
        % It also saves a file that can be used in NERVE statistical output
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
            else
                % fprintf('Komm:%04i;',Vehicle_dist.D1_KommNr(komm))
                if rem(komm,15)==0;fprintf('@%03i/%3i \n',komm,length(Vehicle_dist.D1_KommNr)); end
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
        
        K_EF.distNames  = Trout.Properties.VariableNames;
        K_EF.traffSit   = Trout.Name;
        K_EF.TraffSitID = Trout(:,2:7);
        idl = find(contains(Trout.Properties.VariableNames,'Light'));
        idh = find(contains(Trout.Properties.VariableNames,'Heavy'));
        idb = find(contains(Trout.Properties.VariableNames,'Buses'));
        K_EF.TraffDataL  = table2array(Trout(:,idl));
        K_EF.TraffDataH  = table2array(Trout(:,idl));
        K_EF.TraffDataB  = table2array(Trout(:,idl));
        K_EF.kommNamesL  = Trout.Properties.VariableNames(idl);
        K_EF.kommNamesH  = Trout.Properties.VariableNames(idh);
        K_EF.kommNamesB  = Trout.Properties.VariableNames(idb);
        save(tmpfile,'K_EF')
        fprintf('saved a tmp file \n%s\n',tmpfile)
        
        writetable(Trout,ofile,'Sheet',osheet)
        fprintf('Wrote Sheet: %s \n to file; %s\n',osheet,ofile)
        
        fprintf('%s--- >\n',char(comps(com)))
        
        
        % Use a switch for string as MatLab do not love variable names changes
        % inside loops...
        switch char(comps(com))
            case 'FC'
                OnRoadEF_RoadClasses_FC = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_FC','-append');
            case 'FC_MJ'
                OnRoadEF_RoadClasses_FC_MJ = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_FC_MJ','-append');
            case 'CH4'
                OnRoadEF_RoadClasses_CH4 = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_CH4','-append');
            case 'BC'
                OnRoadEF_RoadClasses_BC = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_BC','-append');
            case 'PM'
                OnRoadEF_RoadClasses_PM = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_PM','-append');
            case 'HC'
                OnRoadEF_RoadClasses_HC = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_HC','-append');
            case 'CO'
                OnRoadEF_RoadClasses_CO = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_CO','-append');
            case 'NOx'
                OnRoadEF_RoadClasses_NOx = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_NOx','-append');
            case 'Be'
                OnRoadEF_RoadClasses_Be = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_Be','-append');
            case 'NMHC'
                OnRoadEF_RoadClasses_NMHC = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_NMHC','-append');
            case 'NO2'
                OnRoadEF_RoadClasses_NO2 = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_NO2','-append');
            case 'NO'
                OnRoadEF_RoadClasses_NO = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_NO','-append');
            case 'PN'
                OnRoadEF_RoadClasses_PN = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_PN','-append');
            case 'CO2'
                OnRoadEF_RoadClasses_CO2 = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_CO2','-append');
            case 'N2O'
                OnRoadEF_RoadClasses_N2O = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_N2O','-append');
            case 'NH3'
                OnRoadEF_RoadClasses_NH3 = Trout;
                save(ofiles.MatlabOutput,'OnRoadEF_RoadClasses_NH3','-append');
        end
    end
    
end

end