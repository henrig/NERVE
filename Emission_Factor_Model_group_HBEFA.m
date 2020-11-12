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
global tfold comps SSB_Vehicle_dist Vehicle_source debug_mode
global EFA EF_AVG roads Vehicle_weight


% global rfid
% global exchfile  debug_mode  SSB_sheet_col HBEFA_temp_Matfile
% global NERVE_EF_file text_div

if ~ismember(Vehicle_source,{'SSB'})
    return
end


fprintf('in Emission_Factor_Model_group_HBEFA *\n')
fprintf('---   Grouping Vehicles according to ---\n---   Model SSB HBEFA merger specifications ---\n%s \n',SSB_Vehicle_dist)

T  = readtable(SSB_Vehicle_dist,'Sheet','HBEFA778toMODEL');
Tm = readtable(SSB_Vehicle_dist,'Sheet','MODEL');
fprintf('read Sheet : %s and Sheet %s\n','HBEFA778toMODEL','MODEL')

% ModelNumber HBEFA_Weight    Uniform_Weight
modelN = unique(T.ModelNumber);
fprintf('Found %i Model Vehicles \n',length(modelN))

% ConVfile = sprintf('%s/EFA_conversion.xlsx',tfold);


% Loop all components (comps) to do conversion for
for com =1:length(comps)
    fprintf('<--- \nSpec| %s\n',char(comps(com)))
    iEFfile = sprintf('%s/EFA_matrix41_RAW_%s.mat',tfold,char(comps(com)));
    oEFfile = sprintf('%s/EFA_Table_MODEL_%s',tfold,char(comps(com)));
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
        
        % check how many non-nan EF there are in each Vehicle.
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
    
    
    switch char(comps(com))
        case 'CO2'
            MunicpalHBEFA_RoadsEF_ = TFout;
            save(ofiles.MatlabOutput,'MunicpalHBEFA_RoadsEF_','roads','-append');
        case 'FC'
        case 'NOx'
        case 'CO2'
        case 'CO2'
        case 'CO2'
    
    end
    
    
    writetable(TFout,'modelHBEFA.xlsx','Sheet',sprintf('%s_%s_Weight',char(comps(com)),Vehicle_weight) )
    fprintf('Saved a temp-file for Emission Factors Model:\n%s\n',oEFfile)

end

% OLD METHOD
% for com =1:length(comps)
%     fprintf('<--- \nSpec| %s\n',char(comps(com)))
%     iEFfile = sprintf('%s/EFA_matrix41_RAW_%s',tfold,char(comps(com)));
%     oEFfile = sprintf('%s/EFA_matrix41_MODEL_%s',tfold,char(comps(com)));
%     load(iEFfile);
%     D = size(EF_AVG);
%     fprintf('EF_AVG -- Dimensions\n')
%     for i=1:length(D)
%         fprintf('EF_AVG D%i  : %i\n',i,D(i));
%     end
%
%     Defines the EFA
%     EFA = nan(length(modelN),size(EF_AVG,2),size(EF_AVG,3),size(EF_AVG,4),size(EF_AVG,5),size(EF_AVG,6));
%     for mod = 1:length(modelN)
%         Tsub  = T(T.ModelNumber==modelN(mod),:);
%
%         Make vehicles with HBEFA weight
%         EFsub = EF_AVG(Tsub.HBEFA_Num,:,:,:,:,:);
%
%         switch Vehicle_weight
%
%             case 'Uniform'
%                 Make vehicles with Uniform weight
%                 Wu = sum(Tsub.Uniform_Weight);
%                 if abs(Wu-1)>1e-5
%                     fprintf('### error Model Vehicle Weight does not add to 1 but %d!\n',Wu)
%                     Tsub
%                 end
%                 for roa = 1:size(EFA,2)
%                     for spd = 1:size(EFA,3)
%                         for dec = 1:size(EFA,4)
%                             for cog = 1:size(EFA,5)
%                                 for urb = 1:size(EFA,6)
%                                     TW = 0;
%                                     has_EF =  0;
%                                     for  i =1:size(EFsub,1)
%                                         if ~isnan(EFsub(i,roa,spd,dec,cog,urb))
%                                             EFA(modelN(mod),roa,spd,dec,cog,urb)= EFA(modelN(mod),roa,spd,dec,cog,urb) + ...
%                                                 EFsub(i,roa,spd,dec,cog,urb)*Tsub.Uniform_Weight(i);
%                                             TW = TW + Tsub.Uniform_Weight(i);
%                                             has_EF =  1;
%                                         end
%
%                                     end
%                                     if has_EF>0;
%                                         fprintf('i=%i;modelN(mod)=%i;roa=%i;spd=%i;dec=%i;cog=%i;urb=%i;  TW = %f\n',i,modelN(mod),roa,spd,dec,cog,urb,TW)
%                                         pause;
%                                     end
%                                     if abs(TW-1)> 1e-3 && has_EF &&  ~isnan(EFA(modelN(mod),roa,spd,dec,cog,urb))
%                                         fprintf('modelN(mod)=%i;roa=%i;spd=%i;dec=%i;cog=%i;urb=%i;  TW = %f\n',modelN(mod),roa,spd,dec,cog,urb,TW)
%                                         EFA(modelN(mod),roa,spd,dec,cog,urb) = EFA(modelN(mod),roa,spd,dec,cog,urb)/TW;
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%
%             case 'HBEFA'
%                 Wh = sum(Tsub.HBEFA_Weight);
%                 if abs(Wh-1)>1e-5
%                     fprintf('### error Model Vehicle Weight does not add to 1 but %d\n',Wh)
%                     Tsub
%                 end
%             case 'NERVE'
%                 Wn = sum(Tsub.NERVE_Weight);
%                 if abs(Wn-1)>1e-5
%                     fprintf('### error Model Vehicle Weight does not add to 1!\n')
%                     Tsub
%                 end
%
%         end % switch
%
%
%         fprintf('#%i EF No NAN EFs found: %i\n',modelN(mod),length(find(~isnan(EFA(modelN(mod),:,:,:,:,:,:)))))
%
%
%
%     end % for mod
%
%     Emission_Factor_Test_HBEFA_to_MODEL_conversion(com,ConVfile)
%
%
%     save(oEFfile,'EFA','T','roads');
%     fprintf('Saved a temp-file for Emission Factors Model:\n%s\n',oEFfile)
%
% end % comps

end

