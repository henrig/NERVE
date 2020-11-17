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
function Emission_Factor_Model_group_HBEFA()
%--------------------------------------------------------------------------

% 22.10.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global tfold comps Vehicle_source debug_mode
global EFA EF_AVG roads Vehicle_weight ofiles input

if ~ismember(Vehicle_source,{'SSB'})
    return
end
fprintf('---------------------------------------------------------------\n')
fprintf('in Emission_Factor_Model_group_HBEFA    *\n')
fprintf('---------------------------------------------------------------\n')
fprintf('---   Grouping Vehicles according to ---\n')
fprintf('---   Model SSB HBEFA merger specifications')
fprintf('---   \n%s \n',input.files.SSB_Vehicle_dist)

T  = readtable(input.files.SSB_Vehicle_dist,'Sheet','HBEFA778toMODEL');
Tm = readtable(input.files.SSB_Vehicle_dist,'Sheet','MODEL');
fprintf('read Sheet : %s and Sheet %s\n','HBEFA778toMODEL','MODEL')

% ModelNumber HBEFA_Weight    Uniform_Weight
modelN = unique(T.ModelNumber);
fprintf('Found %i Model Vehicles \n',length(modelN))

% Loop all components (comps) to do conversion for
for com =1:length(comps)
    fprintf('<--- \nSpec| %s\n',char(comps(com)))
    iEFfile = sprintf('%sEFA_matrix41_RAW_%s.mat',tfold,char(comps(com)));
    oEFfile = sprintf('%sEFA_Table_MODEL_%s.mat',tfold,char(comps(com)));
    ofile   = 'modelHBEFA.xlsx';
    
    if (~input.options.do_preProcessing_HBEFA && input.options.use_temporary_files && exist(oEFfile))
        fprintf('Using saved %s EFA_Table_MODEL file\n',char(comps(com)))
        load(oEFfile)
    else
        
        load(iEFfile);
        D = size(EF_AVG);
        fprintf('EF_AVG -- Dimensions\n')
        for i=1:length(D)
            fprintf('EF_AVG D%i  : %i\n',i,D(i));
        end

        
        % Tables are slow to work with and requires a lot of memory, but given
        % the state of HBEFA, they offer good control of emission factors and
        % what/where they are missing.
        TFout = table;
        % Loop all model vehicles
        for mod = 1:length(modelN)
            
            % Make the subset of emissions/weights needed to calculate a model
            % vehicle group.
            Tsub  = T(T.ModelNumber==modelN(mod),:);
            % Make vehicles with HBEFA weight
            EFsub = EF_AVG(Tsub.HBEFA_Num,:,:,:,:,:);
            
            % Check how many non-nan EF there are in each Vehicle.
            Nef=[];
            for  i =1:size(EFsub,1)
                Nef(i) = length(find(~isnan(EFsub(i,:,:,:,:,:,:))));
            end
            uef = unique(Nef);
            
            % If there are only one amount (1440) of Emission factors
            if length(uef) == 1
                fprintf('MODEL Vehicle: %3i %-42s Found #%i EF',modelN(mod),char(Tm.Name(modelN(mod))),uef)
                Tout = table;
                k = 1;
                for roa = 1:size(EF_AVG,2)
                    for spd = 1:size(EF_AVG,3)
                        for dec = 1:size(EF_AVG,4)
                            for cog = 1:size(EF_AVG,5)
                                for urb = 1:size(EF_AVG,6)
                                    Trout = table;
                                    if ~isnan(EFsub(1,roa,spd,dec,cog,urb))
                                        Trout.Name        = {sprintf('%s/%s-%i/%i%%/%s',char(roads.RoadEnv(urb)),char(roads.RoadType(roa)),roads.RoadSpeeds(spd),roads.RoadGradient(dec),char(roads.Congestion(cog)))};
                                        %fprintf('%s/%s-%i/%i%%/%s\n',char(roads.RoadEnv(urb)),char(roads.RoadType(roa)),roads.RoadSpeeds(spd),roads.RoadGradient(dec),char(roads.Congestion(cog)))
                                        Trout.Nr        = k;
                                        Trout.RoadNum   = roa;
                                        Trout.SpeedNum  = spd;
                                        Trout.DeclNum   = dec;
                                        Trout.CongNum   = cog;
                                        Trout.EnviNum   = urb;
                                        k = k+1;
                                        Tout = [Tout;Trout];
                                    end
                                end
                            end
                        end
                    end
                end
                
                EFac = nan(height(Tout),1);
                Wght = nan(height(Tout),1);
                for cond = 1:height(Tout)
                    for veh = 1:height(Tsub)
                        if veh == 1
                            switch Vehicle_weight
                                case 'Uniform'
                                    EFac(cond) = EFsub(veh,Tout.RoadNum(cond),Tout.SpeedNum(cond),Tout.DeclNum(cond),Tout.CongNum(cond),Tout.EnviNum(cond))*Tsub.Uniform_Weight(veh);
                                    Wght(cond) = Tsub.Uniform_Weight(veh);
                                case 'NERVE'
                                    EFac(cond) = EFsub(veh,Tout.RoadNum(cond),Tout.SpeedNum(cond),Tout.DeclNum(cond),Tout.CongNum(cond),Tout.EnviNum(cond))*Tsub.NERVE_Weight(veh);
                                    Wght(cond) = Tsub.NERVE_Weight(veh);
                                case 'HBEFA'
                                    EFac(cond) = EFsub(veh,Tout.RoadNum(cond),Tout.SpeedNum(cond),Tout.DeclNum(cond),Tout.CongNum(cond),Tout.EnviNum(cond))*Tsub.HBEFA_Weight(veh);
                                    Wght(cond) = Tsub.HBEFA_Weight(veh);
                            end
                        else
                            switch Vehicle_weight
                                case 'Uniform'
                                    EFac(cond) = EFac(cond) + EFsub(veh,Tout.RoadNum(cond),Tout.SpeedNum(cond),Tout.DeclNum(cond),Tout.CongNum(cond),Tout.EnviNum(cond))*Tsub.Uniform_Weight(veh);
                                    Wght(cond) = Wght(cond)+ Tsub.Uniform_Weight(veh);
                                case 'NERVE'
                                    EFac(cond) = EFac(cond) + EFsub(veh,Tout.RoadNum(cond),Tout.SpeedNum(cond),Tout.DeclNum(cond),Tout.CongNum(cond),Tout.EnviNum(cond))*Tsub.NERVE_Weight(veh);
                                    Wght(cond) = Wght(cond)+ Tsub.NERVE_Weight(veh);
                                case 'HBEFA'
                                    EFac(cond) = EFac(cond) + EFsub(veh,Tout.RoadNum(cond),Tout.SpeedNum(cond),Tout.DeclNum(cond),Tout.CongNum(cond),Tout.EnviNum(cond))*Tsub.HBEFA_Weight(veh);
                                    Wght(cond) = Wght(cond)+ Tsub.HBEFA_Weight(veh);
                            end
                        end
                    end
                end
                if debug_mode
                    sprintf('%s_%sWeight',char(comps(com)),Vehicle_Weight)
                    fprintf(' %7.1f/%7.1f/%7.1f (mean/max/min) \n',mean(EFac./Wght),max(EFac./Wght),min(EFac./Wght))
                    fprintf(' %7.1f/%7.1f/%7.1f (mean/max/min) \n',mean(Wght),max(Wght),min(Wght))
                else
                    fprintf('\n')
                end
                
                Tout.EFac = EFac./Wght;
                Tout.Properties.VariableNames(find(ismember(Tout.Properties.VariableNames,{'EFac'}))) = Tm.Name(modelN(mod));
            else
                % if this is the case, there are NaNs in the emission factor
                % and we will get wrong results.
                for i=1:length(uef)
                    fprintf('\n\n\n\n\n\n\n\n\n\n#### Found #%i EF\n',uef(i))
                end
            end
            
            if mod == 1
                % Add the first part of the table as well
                TFout = Tout;
            else
                % Concatonate only the column of the new vehicle data if the
                % align.
                [a,b,c] = intersect(TFout.Name,Tout.Name);
                addCol  = find(ismember(Tout.Properties.VariableNames,Tm.Name(modelN(mod))));
                
                % Not 100% robust test!!!!?
                if length(b)==length(c)
                    TFout = [TFout,Tout(:,addCol)];
                else
                    fprintf('\n\n\n\n\n\n\n\n\n\n#### Found #%i & %i EF\n',a,b)
                end
            end
            
        end
        
        EFrdCond = table2array(TFout(:,8:end));        
        save(oEFfile,'TFout','roads','EFrdCond')
        fprintf('Saved an excel-file for Emission Factors Model:\n%s\n',oEFfile)
        writetable(TFout,ofile,'Sheet',sprintf('%s_%s_Weight',char(comps(com)),Vehicle_weight))
        fprintf('Saved an excel-file for Emission Factors Model:\n%s\n',ofile)

    end

    switch char(comps(com))
        case 'FC'
            MunicpalHBEFA_RoadsEF_FC = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_FC','roads','-append');
        case 'FC_MJ'
            MunicpalHBEFA_RoadsEF_FC_MJ = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_FC_MJ','roads','-append');
        case 'CH4'
            MunicpalHBEFA_RoadsEF_CH4 = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_CH4','roads','-append');
        case 'BC'
            MunicpalHBEFA_RoadsEF_BC = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_BC','roads','-append');
        case 'PM'
            MunicpalHBEFA_RoadsEF_PM = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_PM','roads','-append');
        case 'HC'
            MunicpalHBEFA_RoadsEF_HC = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_HC','roads','-append');
        case 'CO'
            MunicpalHBEFA_RoadsEF_CO = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_CO','roads','-append');
        case 'NOx'
            MunicpalHBEFA_RoadsEF_NOx = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_NOx','roads','-append');
        case 'Be'
            MunicpalHBEFA_RoadsEF_Be = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_Be','roads','-append');
        case 'NMHC'
            MunicpalHBEFA_RoadsEF_NMHC = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_NMHC','roads','-append');
        case 'NO2'
            MunicpalHBEFA_RoadsEF_NO2 = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_NO2','roads','-append');
        case 'NO'
            MunicpalHBEFA_RoadsEF_NO = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_NO','roads','-append');
        case 'PN'
            MunicpalHBEFA_RoadsEF_PN = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_PN','roads','-append');
        case 'CO2'
            MunicpalHBEFA_RoadsEF_CO2 = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_CO2','roads','-append');
        case 'N2O'
            MunicpalHBEFA_RoadsEF_N2O = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_N2O','roads','-append');
        case 'NH3'
            MunicpalHBEFA_RoadsEF_NH3 = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_NH3','roads','-append');
    end
    
end

end
