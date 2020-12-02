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
% The model classes have a static Emission Factor for each driving
% condition. This means that
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
fprintf('Using weight %s \n',Vehicle_weight)
Tout = readtable(input.files.HBEFA_roads,'Sheet','Roads');


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
    else
        fprintf('Calculating New  EFA_Table_MODEL_%s file\n',char(comps(com)))
        
        load(iEFfile);
        D = size(EF_AVG);
        fprintf('EF_AVG -- Dimensions\n')
        for i=1:length(D)
            fprintf(' D%i:%-3i  ',i,D(i));
        end
        fprintf('\n')
        
        % Tables are slow to work with and requires a lot of memory, but given
        % the state of HBEFA, they offer good control of emission factors and
        % what/where they are missing.
        TFout = Tout;
        % Loop all model vehicles
        for mod = 1:length(modelN)
            % Make the subset of emissions/weights needed to calculate a model
            % vehicle group.
            Tsub  = T(T.ModelNumber==modelN(mod),:);
            % Make vehicles with HBEFA weight
            EFsub = EF_AVG(Tsub.HBEFA_Num,:,:,:,:,:);
            
            % Check how many non-nan EF there are in each Vehicle. Missing
            % EF should be as nan and so this filter them out.
            Nef=[];
            for  i =1:size(EFsub,1)
                Nef(i) = length(find(~isnan(EFsub(i,:,:,:,:,:,:))));
            end
            uef = unique(Nef);
            for mr =1:length(uef)
                fprintf('MODEL Vehicle: %3i %-42s Found #%i EF\n',modelN(mod),char(Tm.Name(modelN(mod))),uef(mr))
            end
            
            % If there are only one amount (1440 or No) of Emission factors
            % if length(uef) == 1
            EFac = zeros(height(Tout),1);
            Wght = zeros(height(Tout),1);
            for veh = 1:height(Tsub)
                switch Vehicle_weight
                    case 'Uniform'
                        VW = Tsub.Uniform_Weight(veh);
                        SW = sum(Tsub.Uniform_Weight);
                    case 'NERVE'
                        VW = Tsub.NERVE_Weight(veh);
                        SW = sum(Tsub.NERVE_Weight);
                    case 'HBEFA'
                        VW = Tsub.HBEFA_Weight(veh);
                        SW = sum(Tsub.HBEFA_Weight);
                end
                fprintf('%i of %i %s Weight %3.2f %-30s      of %3.2f \n',veh,height(Tsub),Vehicle_weight,VW,char(Tsub.Name(veh)),SW)

                for cond = 1:height(Tout)
                    EFcond(cond,1) = EFsub(veh,Tout.RoadNum(cond),Tout.SpeedNum(cond),Tout.DeclNum(cond),Tout.CongNum(cond),Tout.EnviNum(cond));
                end
                
                if ~isnan(sum(EFcond))
                    fprintf('Found 1440 emission factors for vehicle\n')
                    EFac = EFac + EFcond*VW;
                    Wght = Wght + VW;
                elseif ~isnan(nansum(EFcond))
                    fprintf('Found some emission factors for vehicle\n')
                elseif isnan(nansum(EFcond))
                    fprintf('Found no emission factors for vehicle\n')
                end
            end
            nanmean(EFac)
            Tout.EFac = EFac;
            pos = find(ismember(Tout.Properties.VariableNames,{'EFac'}));
            Tout.Properties.VariableNames(pos) = Tm.Name(modelN(mod));
            TFout = [TFout,Tout(:,pos)];
            
            if debug_mode
                fprintf(' %s_%s_weight',char(comps(com)),Vehicle_weight)
                fprintf(' %7.1f/%7.1f/%7.1f (mean/max/min)   ',mean(EFac./Wght),max(EFac./Wght),min(EFac./Wght))
                fprintf(' %7.1f/%7.1f/%7.1f (mean/max/min) \n',mean(Wght),max(Wght),min(Wght))
            else
                fprintf('\n')
            end
        end
        
        EFrdCond = table2array(TFout(:,8:end));
        save(oEFfile,'TFout','roads','EFrdCond')
        fprintf('Saved an temp .mat -file for Emission Factors Model:\n%s\n',oEFfile)
        
        % writetable(TFout,ofile,'Sheet',sprintf('%s_%s_Weight',char(comps(com)),Vehicle_weight))
        % fprintf('Saved an excel-file for Emission Factors Model:\n%s\n',ofile)
        
    end
end

end
