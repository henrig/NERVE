function Emission_Factor_Model_group_HBEFA()
%--------------------------------------------------------------------------

% 22.10.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global tfold comps SSB_Vehicle_dist Vehicle_source


% global rfid
% global exchfile  debug_mode  SSB_sheet_col HBEFA_temp_Matfile
% global NERVE_EF_file text_div

if ~ismember(Vehicle_source,{'SSB'})
    return
end

Vehicle_weight = 'Uniform';
    
    
fprintf('\tin Emission_Factor_Model_group_HBEFA\n')
fprintf('---   Grouping Vehicles according to \n---   Model SSB HBEFA merger specifications \n\n')

T = readtable(SSB_Vehicle_dist,'Sheet','HBEFA778toMODEL');

% ModelNumber HBEFA_Weight    Uniform_Weight
modelN = unique(T.ModelNumber);
fprintf('Found %i Model Vehicles \n',length(modelN))

for com =1:length(comps)
    fprintf('Spec| %s\n',char(comps(com)))
    iEFfile = sprintf('%s/EFA_matrix41_RAW_%s',tfold,char(comps(com)));
    oEFfile = sprintf('%s/EFA_matrix41_MODEL_%s',tfold,char(comps(com)));
    load(iEFfile);
    D = size(EF_AVG);
    fprintf('EF_AVG -- Dimensions\n')
    for i=1:length(D)
        fprintf('EF_AVG D%i  : %i\n',i,D(i));
    end
    
    % Defines the EFA
    EFA = zeros(length(modelN),size(EF_AVG,2),size(EF_AVG,3),size(EF_AVG,4),size(EF_AVG,5),size(EF_AVG,6));
    for mod=1:length(modelN)
        Tsub  = T(T.ModelNumber==modelN(mod),:);
        %##################
        % errrrrorrrr  Need to add column "HBEFA_Num" to the sheet!!!!!
        %Tsub.HBEFA_Num = (1:height(Tsub))';
        %##################
        % Make vehicles with HBEFA weight
        EFsub = EF_AVG(Tsub.HBEFA_Num,:,:,:,:,:);
        switch Vehicle_weight
            case 'HBEFA'
                Wh = sum(Tsub.HBEFA_Weight);
                if abs(Wh-1)>1e-5
                    fprintf('### error Model Vehicle Weight does not add to 1 but %d\n',Wh)
                    Tsub
                end
                for  i =1:size(EFsub,1)
                    EFA(modelN(mod),:,:,:,:,:)= EFA(modelN(mod),:,:,:,:,:) + ...
                        EFsub(i,:,:,:,:,:)*Tsub.HBEFA_Weight(i);
                end
                
            case 'Uniform'
                % Make vehicles with Uniform weight
                Wu = sum(Tsub.Uniform_Weight);
                if abs(Wu-1)>1e-5
                    fprintf('### error Model Vehicle Weight does not add to 1 but %d!\n',Wu)
                    %Tsub
                end
                for  i =1:size(EFsub,1)
                    EFA(modelN(mod),:,:,:,:,:)= EFA(modelN(mod),:,:,:,:,:) + ...
                        EFsub(i,:,:,:,:,:)*Tsub.HBEFA_Weight(i);
                end
                
            case 'NERVE'
                Wn = sum(Tsub.NERVE_Weight);
                if abs(Wn-1)>1e-5
                    fprintf('### error Model Vehicle Weight does not add to 1!\n')
                    Tsub
                end
                for  i =1:size(EFsub,1)
                    EFA(modelN(mod),:,:,:,:,:)= EFA(modelN(mod),:,:,:,:,:) + ...
                        EFsub(i,:,:,:,:,:)*Tsub.NERVE_Weight(i);
                end
        end % switch
    end % for mod

    save(oEFfile,'EFA','T');
    fprintf('Saved a temp-file for Emission Factors Model:\n%s\n',oEFfile)
    
end % comps
 
end

