function NERVE_Emissions_Statistical_Oputput()
% Data needed by the MijlDirektoratet sheet
comps  = [{'N2O'}];
comps  = [{'CO2'},{'FC'},{'CH4'},{'N2O'}];
comps  = [{'CO2'},{'FC'},{'CH4'}];
comps  = [{'CO2'},{'FC'}];
%x NV(kom,veh) %#
%x TD(kom,veh) %km
%x ORdComp(kom,veh) % frac
%x FRdComp(kom,veh) % frac
%x L_IN(kom)        % km/yr
%x H_IN(kom)        % km/yr
%x B_IN(kom)        % km/yr

%Check the below values against Roads: 
%s L_FROM(kom)       % km/yr
%s H_FROM(kom)       % km/yr
%s B_FROM(kom)       % km/yr
%s EF_IN(kom,veh)    % g/km
%s EM_IN(kom,veh)    % g/yr
%s EF_FROM(kom,veh)  % g/km
%s EM_FROM(kom,veh)  % g/yr


% Files that must be read once only
% 'Model_Vehicles_Merge_SSB_and_HBEFA_Vehicles.xlsx'
TM = readtable('Input/Model_Vehicles_Merge_SSB_and_HBEFA_Vehicles.xlsx','Sheet','MODEL');
LightVehiclesIdx = TM.ClassNum==1|TM.ClassNum==2;
BusesVehiclesIdx = TM.ClassNum==3|TM.ClassNum==4;
HeavyVehiclesIdx = TM.ClassNum==5|TM.ClassNum==6|TM.ClassNum==7;

for com = 1:length(comps)
    % files that must be read per species
    R_EF_File = sprintf('Temp/EFA_Table_MODEL_%s.mat',char(comps(com)));
    load(R_EF_File,'TFout'); % TFout

    fprintf('%s\n',char(comps(com)))

    for Tyear = 2009:2019
        fprintf('%i\n',Tyear)

        % files that must be read per year
        munFile = sprintf('Temp/Municipal_Traffic_Exchange_%i.mat',Tyear);
        RdDistFile = sprintf('Output/RoadTypeDistanceMunicipal%i.mat',Tyear);
        load(munFile)
        load(RdDistFile) % KDD
        
        % files that must be read for each year and component!
        EF_File = sprintf('Temp/EF_On_AllRoadCond_Municipality_%i_%s.mat',Tyear,char(comps(com)));
        load(EF_File)
                
        RLinks = shaperead(sprintf('Output/Traffic_Emissions_%i',Tyear));

        DaysInYear = datenum([Tyear+1 1 1 0 0 0])-datenum([Tyear 1 1 0 0 0]);
        
        KOMM       = extractfield(RLinks,'KOMMS');
        KOMMe      = extractfield(RLinks,'KOMME');
        IDO        = extractfield(RLinks,'IDO');
        L_ADT      = extractfield(RLinks,sprintf('L_ADT%i',Tyear));
        H_ADT      = extractfield(RLinks,sprintf('H_ADT%i',Tyear));
        B_ADT      = extractfield(RLinks,sprintf('B_ADT%i',Tyear));
        DISTANCE   = extractfield(RLinks,'DISTANCE');
        LW         = DaysInYear*L_ADT.*DISTANCE;
        HW         = DaysInYear*H_ADT.*DISTANCE;
        BW         = DaysInYear*B_ADT.*DISTANCE;
        ND_NERVE_L = sum(LW);
        ND_NERVE_H = sum(HW);
        ND_NERVE_B = sum(BW);
        ND_NERVE   = ND_NERVE_L + ND_NERVE_H + ND_NERVE_B;

        % Ø----------------------------------------------------------------
        % find *L_IN(kom)*
        ukomm = unique(KOMM);
        for k=1:length(ukomm)
            f  = find(KOMM  == ukomm(k));
            f2 = find((KOMMe == ukomm(k))&(KOMM ~= ukomm(k)) );
            L_IN(k) = sum(LW(f).*IDO(f))+sum(LW(f2).*(1-IDO(f2)));
            H_IN(k) = sum(HW(f).*IDO(f))+sum(HW(f2).*(1-IDO(f2)));
            B_IN(k) = sum(BW(f).*IDO(f))+sum(BW(f2).*(1-IDO(f2)));
        end
        
        L_EM      = extractfield(RLinks,sprintf('EM_L_%s',char(comps(com))));
        H_EM      = extractfield(RLinks,sprintf('EM_H_%s',char(comps(com))));
        B_EM      = extractfield(RLinks,sprintf('EM_B_%s',char(comps(com))));
        
        EM = L_EM+H_EM+B_EM; 
        for k=1:length(ukomm)
            f  = find(KOMM  == ukomm(k));
            f2 = find((KOMMe == ukomm(k))&(KOMM ~= ukomm(k)) );
            EM_INr(k) = sum(EM(f).*IDO(f))+sum(EM(f2).*(1-IDO(f2)));
            EM_INr(k) = sum(EM(f).*IDO(f))+sum(EM(f2).*(1-IDO(f2)));
            EM_INr(k) = sum(EM(f).*IDO(f))+sum(EM(f2).*(1-IDO(f2)));
        end
        
        % Vehicle_dist
        %------------------------------------------------------------------
        % find *NV(kom,veh)* from:: SSB data
        % find *TD(kom,veh)* from:: SSB data
        % find *OrdComp(kom,veh)* from:: SSB data
        % find *FrdComp(kom,veh)* from:: SSB data
        NV  = Vehicle_dist.modelNV;
        TD  = Vehicle_dist.modelTD;
        ORdComp = Vehicle_dist.Vdist;
        FRdComp = Vehicle_dist.VdistFROM;        
        
        %------------------------------------------------------------------
        % Estimate *EF_IN(kom,Veh)  * based on:::: EF_IN and exchange
        
        KDD_L = KDD.TraffDataL;
        KDD_H = KDD.TraffDataH;
        KDD_B = KDD.TraffDataB;
        
        % DOES NOT INCLUDE BIO!
        
        RdEFs = table2array(TFout(:,8:end));       
        vehicles = TFout.Properties.VariableNames(8:end);
        for k=1:length(ukomm)
            LtrafficSitW_IN = KDD_L(:,k)/sum(KDD_L(:,k));
            HtrafficSitW_IN = KDD_H(:,k)/sum(KDD_H(:,k));
            BtrafficSitW_IN = KDD_B(:,k)/sum(KDD_B(:,k));
            for veh = 1:length(vehicles)
                if (TM.Model_Class(veh) ==1 || TM.Model_Class(veh) ==2)
                    EF_IN(k,veh) = sum(RdEFs(:,veh).*LtrafficSitW_IN);
                elseif(TM.Model_Class(veh) ==5 ||TM.Model_Class(veh) ==6 ||TM.Model_Class(veh) ==7)
                    EF_IN(k,veh) = sum(RdEFs(:,veh).*HtrafficSitW_IN);
                elseif(TM.Model_Class(veh) ==3 ||TM.Model_Class(veh) ==4)
                    EF_IN(k,veh) = sum(RdEFs(:,veh).*BtrafficSitW_IN);
                else
                    fprintf('fail\n')
                end
            end
        end
        
        %------------------------------------------------------------------
        % Estimate *EM_IN(kom,Veh)  * based on:::: EF_IN and exchange
        
        for k=1:length(ukomm)
            for veh = 1:length(vehicles)
                type = TM.Model_Class(veh);
                if type <= 2
                    EM_IN(k,veh) = L_IN(k)*ORdComp(k,veh)*EF_IN(k,veh);
                elseif (type == 5 || type == 6 || type == 7)
                    EM_IN(k,veh) = H_IN(k)*ORdComp(k,veh)*EF_IN(k,veh);
                elseif (type == 3 || type == 4)
                    EM_IN(k,veh) = B_IN(k)*ORdComp(k,veh)*EF_IN(k,veh);
                end
            end
        end

        
% END IN
%--------------------------------------------------------------------------
% This method assigns the total Emissions based on the exchange.

        for komm = 1:size(TrafficFROM,2)
            I            = find(TrafficFROM(komm,:)>0);
            Tshare       = (TrafficFROM(komm,I)/100);
            L_FROM(komm) = nansum(L_IN(I).*Tshare);
            H_FROM(komm) = nansum(H_IN(I).*Tshare);
            B_FROM(komm) = nansum(B_IN(I).*Tshare);
        end
        
        EF_FROM=zeros(size(EF_IN));
        for komm=1:size(EF_IN,1)
            for Veh=1:size(EF_IN,2)
                 EF_FROM(komm,Veh) = EF_FROM(komm,Veh)+nansum(EF_IN(:,Veh).*(TrafficFROM(:,komm)/100));
             end
        end
 
        %--------------------------------------------------------------------------
        % create *L_FROM(kom)* based on:::: L_IN and exchange
        for k=1:length(ukomm)
            for veh = 1:length(vehicles)
                type = TM.Model_Class(veh);
                if type <= 2
                    EM_FROM(k,veh) = L_FROM(k)*FRdComp(k,veh)*EF_FROM(k,veh);
                elseif (type == 5 || type == 6 || type == 7)
                    EM_FROM(k,veh) = H_FROM(k)*FRdComp(k,veh)*EF_FROM(k,veh);
                elseif (type == 3 || type == 4)
                    EM_FROM(k,veh) = B_FROM(k)*FRdComp(k,veh)*EF_FROM(k,veh);
                end
            end
        end
        

        
        
        
        fileout = sprintf('Output/NERVE_output_%s_%04i.mat',char(comps(com)),Tyear);
        fprintf('Processed Emissions for %s year %i\n',char(comps(com)),Tyear)
        save(fileout,'NV','TD','L_IN','H_IN','B_IN','L_FROM','H_FROM','B_FROM','ORdComp','FRdComp','EF_IN','EF_FROM','EM_IN','EM_FROM')
    end
end


end