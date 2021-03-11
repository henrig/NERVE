function NERVE_Emissions_Statistical_Oputput()
% Data needed by the MijlDirektoratet sheet
% comps  = [{'N2O'}];
comps  = [{'CO2'},{'FC'},{'CH4'},{'N2O'}];
% comps  = [{'CO2'},{'FC'},{'CH4'}];
% comps  = [{'CO2'},{'FC'}];
% comps  = [{'CH4'}];

% x NV(kom,veh)      % #
% x TD(kom,veh)      % km
% x ORdComp(kom,veh) % frac
% x FRdComp(kom,veh) % frac
% x L_IN(kom)        % km/yr
% x H_IN(kom)        % km/yr
% x B_IN(kom)        % km/yr

%Check the below values against Roads:
% s L_FROM(kom)       % km/yr
% s H_FROM(kom)       % km/yr
% s B_FROM(kom)       % km/yr
% s EF_IN(kom,veh)    % g/km
% s EM_IN(kom,veh)    % g/yr
% s EF_FROM(kom,veh)  % g/km
% s EM_FROM(kom,veh)  % g/yr

if ispc
    bp = 'N:\Inby';
else
    bp = '/storage/nilu/Inby';
end
pathn = sprintf('%s/Emission_Group/Emission_Models/HEDGE/',bp);
tfold = sprintf('%sTemp/'  ,pathn);
ofold = sprintf('%sOutputNewCong/',pathn);
ifold = sprintf('%sInput/'  ,pathn);
pfold = sprintf('%sPlots/'  ,pathn);


use_ido = 0;
makeplots = 1;

% Files that must be read once only
% 'Model_Vehicles_Merge_SSB_and_HBEFA_Vehicles.xlsx'
model_class_file = sprintf('%sModel_Vehicles_Merge_SSB_and_HBEFA_Vehicles.xlsx',ifold);
TM               = readtable(model_class_file,'Sheet','MODEL');
LightVehiclesIdx = TM.ClassNum==1|TM.ClassNum==2;
BusesVehiclesIdx = TM.ClassNum==3|TM.ClassNum==4;
HeavyVehiclesIdx = TM.ClassNum==5|TM.ClassNum==6|TM.ClassNum==7;

PCVehiclesIdx  = TM.ClassNum==1;
LCVVehiclesIdx = TM.ClassNum==2;
BusesVehiclesIdx = TM.ClassNum==3|TM.ClassNum==4;
ATVehiclesIdx = TM.ClassNum==6|TM.ClassNum==7;
RTVehiclesIdx = TM.ClassNum==5;

for com = 1:length(comps)
    % files that must be read per species
    
    fprintf('\n\n%s\n',char(comps(com)))
    
    for Tyear = 2009:2019
        % files that must be read for each year and component!
        EF_File = sprintf('%sEF_On_AllRoadCond_Municipality_%i_%s.mat',tfold,Tyear,char(comps(com)));
        load(EF_File)
        fileout = sprintf('%sNERVE_output_%s_%04i.mat',ofold,char(comps(com)),Tyear);
        
        try % this should only work for CO2
            R_EF_File = sprintf('%sEFA_Table_MODEL_%s_Bio%i.mat',tfold,char(comps(com)),Tyear);
            load(R_EF_File,'TFout'); % TFout
            fprintf('Found Bio File! %s',R_EF_File)
        catch % all non CO2 species
            R_EF_File = sprintf('%sEFA_Table_MODEL_%s.mat',tfold,char(comps(com)));
            load(R_EF_File,'TFout'); % TFout
            fprintf('No Bio File! %s',R_EF_File)
        end
        fprintf('\n%i\n',Tyear)
        
        % files that must be read per year
        munFile    = sprintf('%sMunicipal_Traffic_Exchange_%i.mat',tfold,Tyear);
        SSBCPfile  = sprintf('%sSSB_CarPark_%i.mat',tfold,Tyear);
        RdDistFile = sprintf('%sRoadTypeDistanceMunicipal%i.mat',ofold,Tyear);
        load(munFile)    % TrafficIN TrafficFROM kmne
        load(RdDistFile) % KDD
        load(SSBCPfile)  % Vehicle_dist
        
        fprintf('Reading large RoadLink file...')
        RLinks = shaperead(sprintf('%sOutput2/Traffic_Emissions_%i',pathn,Tyear));
        fprintf('Done\n')
        
        DaysInYear = datenum([Tyear+1 1 1 0 0 0])-datenum([Tyear 1 1 0 0 0]);
        
        KOMM       = extractfield(RLinks,'KOMMS');
        KOMMe      = extractfield(RLinks,'KOMME');
        KOMMm      = extractfield(RLinks,'KOMM');
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
            if use_ido
                f  = find(KOMM  == ukomm(k));
                f2 = find((KOMMe == ukomm(k))&(KOMM ~= ukomm(k)));
                L_IN(k) = sum(LW(f).*IDO(f))+sum(LW(f2).*(1-IDO(f2)));
                H_IN(k) = sum(HW(f).*IDO(f))+sum(HW(f2).*(1-IDO(f2)));
                B_IN(k) = sum(BW(f).*IDO(f))+sum(BW(f2).*(1-IDO(f2)));
            else
                f  = find(KOMMm  == ukomm(k));
                L_IN(k) = sum(LW(f));
                H_IN(k) = sum(HW(f));
                B_IN(k) = sum(BW(f));
            end
        end
        fprintf('Trafikkarbeid Road : Oppsummert \n')
        fprintf('Light %12.0f : %12.0f %4.2f%%\n',sum(LW),sum(L_IN),100*sum(LW)/sum(L_IN))
        fprintf('Heavy %12.0f : %12.0f %4.2f%%\n',sum(HW),sum(H_IN),100*sum(HW)/sum(H_IN))
        fprintf('Buses %12.0f : %12.0f %4.2f%%\n',sum(BW),sum(B_IN),100*sum(BW)/sum(B_IN))
        
        L_EM      = extractfield(RLinks,sprintf('EM_L_%s',char(comps(com))));
        H_EM      = extractfield(RLinks,sprintf('EM_H_%s',char(comps(com))));
        B_EM      = extractfield(RLinks,sprintf('EM_B_%s',char(comps(com))));
        
        EM = L_EM+H_EM+B_EM;
        for k=1:length(ukomm)
            if use_ido
                f  = find(KOMM  == ukomm(k));
                f2 = find((KOMMe == ukomm(k))&(KOMM ~= ukomm(k)) );
                EM_INr(k) = sum(EM(f).*IDO(f))+sum(EM(f2).*(1-IDO(f2)));
            else
                f  = find(KOMMm  == ukomm(k));
                EM_INr(k) = sum(EM(f));
            end
        end
        fprintf('Utslipp Road : Oppsummert \n----------------------\n')
        pdata =1e-9*[sum(L_EM),sum(H_EM),sum(B_EM),sum(L_EM+H_EM+B_EM)];
        fprintf('Roads LIGHT:%6.2f,  HEAVY:%6.2f,  BUS:%6.2f, TOT RDs : %8.3f\n',pdata)
        % Vehicle_dist
        %------------------------------------------------------------------
        % find *NV(kom,veh)* from:: SSB data
        % find *TD(kom,veh)* from:: SSB data
        % find *OrdComp(kom,veh)* from:: SSB data
        % find *FrdComp(kom,veh)* from:: SSB data
        NV      = Vehicle_dist.modelNV;
        TD      = Vehicle_dist.modelTD;
        ORdComp = Vehicle_dist.Vdist;
        FRdComp = Vehicle_dist.VdistFROM;
        
        %------------------------------------------------------------------
        % Estimate *EF_IN(kom,Veh)  * based on:::: EF_IN and exchange
        
        KDD_L = KDD.TraffDataL;
        KDD_H = KDD.TraffDataH;
        KDD_B = KDD.TraffDataB;
        
        % DOES NOT INCLUDE BIO!
        RdEFs    = table2array(TFout(:,8:end));
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
            Ltcomp = sum(ORdComp(k,LightVehiclesIdx));
            Htcomp = sum(ORdComp(k,HeavyVehiclesIdx));
            Btcomp = sum(ORdComp(k,BusesVehiclesIdx));
            
            for veh = 1:length(vehicles)
                type = TM.Model_Class(veh);
                if type <= 2
                    EM_IN(k,veh) = L_IN(k)*ORdComp(k,veh)/Ltcomp*EF_IN(k,veh);
                elseif (type == 5 || type == 6 || type == 7)
                    EM_IN(k,veh) = H_IN(k)*ORdComp(k,veh)/Htcomp*EF_IN(k,veh);
                elseif (type == 3 || type == 4)
                    EM_IN(k,veh) = B_IN(k)*ORdComp(k,veh)/Btcomp*EF_IN(k,veh);
                end
            end
        end
        pdata =[1e-9*nansum(nansum(EM_IN(:,LightVehiclesIdx ))),1e-9*nansum(nansum(EM_IN(:,HeavyVehiclesIdx ))),1e-9*nansum(nansum(EM_IN(:,BusesVehiclesIdx ))),1e-9*nansum(nansum(EM_IN)),100*nansum(nansum(EM_IN))/sum(L_EM+H_EM+B_EM)];
        fprintf('NERVE LIGHT:%6.2f,  HEAVY:%6.2f,  BUS:%6.2f, TOT_IN  : %8.3f TOT:diff %4.2f%%\n',pdata)
        
        % END IN
        %------------------------------------------------------------------
        % This method assigns the total Emissions based on the exchange.
        for komm = 1:size(TrafficFROM,2)
            I            = find(TrafficFROM(komm,:)>0);
            Tshare       = (TrafficFROM(komm,I)/100);
            L_FROM(komm) = nansum(L_IN(I).*Tshare);
            H_FROM(komm) = nansum(H_IN(I).*Tshare);
            B_FROM(komm) = nansum(B_IN(I).*Tshare);
        end
        
        EF_FROM = zeros(size(EF_IN));
        for komm=1:size(EF_IN,1)
            for Veh=1:size(EF_IN,2)
                EF_FROM(komm,Veh) = EF_FROM(komm,Veh)+nansum(EF_IN(:,Veh).*(TrafficFROM(:,komm)/100));
            end
        end
        %------------------------------------------------------------------
        % create *L_FROM(kom)* based on :::: L_IN and exchange
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
        pdata =[1e-9*nansum(nansum(EM_FROM)),100*nansum(nansum(EM_FROM))/sum(L_EM+H_EM+B_EM)];
        pdata =[1e-9*nansum(nansum(EM_FROM(:,LightVehiclesIdx ))),1e-9*nansum(nansum(EM_FROM(:,HeavyVehiclesIdx ))),1e-9*nansum(nansum(EM_FROM(:,BusesVehiclesIdx ))),1e-9*nansum(nansum(EM_FROM)),100*nansum(nansum(EM_FROM))/sum(L_EM+H_EM+B_EM)];
        fprintf('NERVE LIGHT:%6.2f,  HEAVY:%6.2f,  BUS:%6.2f, TOT_FROM: %8.3f TOT:diff %4.2f%%\n',pdata)

        
        if  makeplots
            f = figure('visible','off');
            scatter(L_IN,L_FROM,sum(NV,2).^.7,ukomm,'filled')
            set(gca,'XScale','log','YScale','log');
            cbh = colorbar;
            cbh.Label.String = 'Kommunenummer';
            hold on
            plot([min([min(L_IN),min(L_FROM)]),max([max(L_IN),max(L_FROM)])],[min([min(L_IN),min(L_FROM)]),max([max(L_IN),max(L_FROM)])],'k--','LineWidth',3)
            grid on
            ylabel('Distanse kjørt med biler med opphave i kommune x')
            xlabel('Distanse kjørt på veiene i kommune x')
            tn = sprintf('Trafikkarbeid Alle Biler og kommuner %i',Tyear   );
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'IN_FROM.jpg'))
            
            f = figure('visible','off');
            scatter(L_FROM,sum(TD,2),sum(NV,2).^.7,ukomm,'filled')
            set(gca,'XScale','log','YScale','log');
            cbh = colorbar;
            cbh.Label.String = 'Kommunenummer';
            hold on
            plot([min([min(L_IN),min(L_FROM)]),max([max(L_IN),max(L_FROM)])],[min([min(L_IN),min(L_FROM)]),max([max(L_IN),max(L_FROM)])],'k--','LineWidth',3)
            grid on
            xlabel('Distanse kjørt med biler med opphave i kommune x')
            ylabel('Distanse kjørt med biler (SSB)registrert i kommune x')
            tn = sprintf('Trafikkarbeid Alle Biler og kommuner %i',Tyear   );
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'FROM_SSBDD.jpg'))
            
            f = figure('visible','off');
            scatter((L_IN+H_IN+B_IN)',nansum(EM_IN,2)./(L_IN+H_IN+B_IN)',sum(NV,2).^.7,ukomm,'filled')
            set(gca,'XScale','log');
            cbh = colorbar;
            cbh.Label.String = 'Kommunenummer';
            grid on
            xlabel('Distanse kjørt med biler med opphave i kommune x')
            ylabel('Utslippsfaktor (g/km)')
            tn = sprintf('Trafikkarbeid og Utslippsfaktor Alle Biler og kommuner %i %s',Tyear,char(comps(com)));
            hold on
            scatter((L_IN+H_IN+B_IN)',nansum(EM_IN,2)./(L_IN+H_IN+B_IN)',5,'w','filled')
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            
            f = figure('visible','off');
            scatter((L_IN)',nansum(EM_IN(:,LightVehiclesIdx),2)./(L_IN)',sum(NV,2).^.7,ukomm,'filled')
            set(gca,'XScale','log');
            cbh = colorbar;
            cbh.Label.String = 'Kommunenummer';
            grid on
            xlabel('Distanse kjørt med biler med opphave i kommune x')
            ylabel('Utslippsfaktor (g/km)')
            tn = sprintf('Trafikkarbeid og Utslippsfaktor Lette Biler %i %s',Tyear,char(comps(com)));
            hold on
            scatter((L_IN)',nansum(EM_IN(:,LightVehiclesIdx),2)./(L_IN)',5,'w','filled')
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            
            f = figure('visible','off');
            scatter((H_IN)',nansum(EM_IN(:,HeavyVehiclesIdx),2)./(H_IN)',sum(NV,2).^.7,ukomm,'filled')
            set(gca,'XScale','log');
            cbh = colorbar;
            cbh.Label.String = 'Kommunenummer';
            grid on
            xlabel('Distanse kjørt med biler med opphave i kommune x')
            ylabel('Utslippsfaktor (g/km)')
            tn = sprintf('Trafikkarbeid og Utslippsfaktor Tunge Kjøretøy %i %s',Tyear,char(comps(com)));
            hold on
            scatter((H_IN)',nansum(EM_IN(:,HeavyVehiclesIdx),2)./(H_IN)',5,'w','filled')
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            
            f = figure('visible','off');
            scatter((B_IN)',nansum(EM_IN(:,BusesVehiclesIdx),2)./(B_IN)',sum(NV,2).^.7,ukomm,'filled')
            set(gca,'XScale','log');
            cbh = colorbar;
            cbh.Label.String = 'Kommunenummer';
            grid on
            xlabel('Distanse kjørt med biler med opphave i kommune x')
            ylabel('Utslippsfaktor (g/km)')
            tn = sprintf('Trafikkarbeid og Utslippsfaktor Busser %i %s',Tyear,char(comps(com)));
            hold on
            scatter((B_IN)',nansum(EM_IN(:,BusesVehiclesIdx),2)./(B_IN)',5,'w','filled')
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            
            close all

            f = figure('visible','off');
            bh1 = boxplot(EF_FROM(:,PCVehiclesIdx),'positions',(1:sum(PCVehiclesIdx))-.2,'symbol','r.','Color',[.4 .4 .4]);
            hold on
            bh2 = boxplot(EF_IN  (:,PCVehiclesIdx),'positions',(1:sum(PCVehiclesIdx))+.2,'symbol','b.');
            set(gca,'XTick',1:sum(PCVehiclesIdx),'XTickLabel',TM.Name(PCVehiclesIdx),'XTickLabelRotation',45)
            grid on
            tn = sprintf('Utslippsfaktor (EF IN) Personbil %i  %s',Tyear,char(comps(com)) );
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            text('Position',[10.0 0],'Interpreter','latex','String','EM IN','Color','b','FontSize',26);
            text('Position',[20.0 0],'Interpreter','latex','String','EM FROM','Color',[.4 .4 .4],'FontSize',26);
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            
            f = figure('visible','off');
            bh1 = boxplot(EF_FROM(:,LCVVehiclesIdx),'positions',(1:sum(LCVVehiclesIdx))-.2,'symbol','r.','Color',[.4 .4 .4]);
            hold on
            bh2 = boxplot(EF_IN  (:,LCVVehiclesIdx),'positions',(1:sum(LCVVehiclesIdx))+.2,'symbol','b.');
            set(gca,'XTick',1:sum(LCVVehiclesIdx),'XTickLabel',TM.Name(LCVVehiclesIdx),'XTickLabelRotation',45)
            grid on
            tn = sprintf('Utslippsfaktor (EF IN) Varebil %i %s',Tyear,char(comps(com)));
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            text('Position',[10.0 0],'Interpreter','latex','String','EM IN','Color','b','FontSize',26);
            text('Position',[20.0 0],'Interpreter','latex','String','EM FROM','Color',[.4 .4 .4],'FontSize',26);
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            
            f = figure('visible','off');
            bh1 = boxplot(EF_FROM(:,BusesVehiclesIdx),'positions',(1:sum(BusesVehiclesIdx))-.2,'symbol','r.','Color',[.4 .4 .4]);
            hold on
            bh2 = boxplot(EF_IN  (:,BusesVehiclesIdx),'positions',(1:sum(BusesVehiclesIdx))+.2,'symbol','b.');
            set(gca,'XTick',1:sum(BusesVehiclesIdx),'XTickLabel',TM.Name(BusesVehiclesIdx),'XTickLabelRotation',45)
            grid on
            tn = sprintf('Utslippsfaktor (EF IN) Busser %s %i',char(comps(com)),Tyear   );
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            text('Position',[10.0 0],'Interpreter','latex','String','EM IN','Color','b','FontSize',26);
            text('Position',[20.0 0],'Interpreter','latex','String','EM FROM','Color',[.4 .4 .4],'FontSize',26);
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            
            f = figure('visible','off');
            bh1 = boxplot(EF_FROM(:,RTVehiclesIdx),'positions',(1:sum(RTVehiclesIdx))-.2,'symbol','r.','Color',[.4 .4 .4]);
            hold on
            bh2 = boxplot(EF_IN  (:,RTVehiclesIdx),'positions',(1:sum(RTVehiclesIdx))+.2,'symbol','b.');
            set(gca,'XTick',1:sum(RTVehiclesIdx),'XTickLabel',TM.Name(RTVehiclesIdx),'XTickLabelRotation',45)
            tn = sprintf('Utslippsfaktor (EF IN) Lastebiler %s %i',char(comps(com)),Tyear   );
            grid on
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            text('Position',[10.0 0],'Interpreter','latex','String','EM IN','Color','b','FontSize',26);
            text('Position',[20.0 0],'Interpreter','latex','String','EM FROM','Color',[.4 .4 .4],'FontSize',26);
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            
            f = figure('visible','off');
            bh1 = boxplot(EF_FROM(:,ATVehiclesIdx),'positions',(1:sum(ATVehiclesIdx))-.2,'symbol','r.','Color',[.4 .4 .4]);
            hold on
            bh2 = boxplot(EF_IN  (:,ATVehiclesIdx),'positions',(1:sum(ATVehiclesIdx))+.2,'symbol','b.');
            set(gca,'XTick',1:sum(ATVehiclesIdx),'XTickLabel',TM.Name(ATVehiclesIdx),'XTickLabelRotation',45)
            grid on
            tn = sprintf('Utslippsfaktor (EF IN) Trekkbiler %s %i',char(comps(com)),Tyear   );
            title(tn)
            x= gcf;
            x.Position=[1986 1 2495 1361];
            text('Position',[10.0 0],'Interpreter','latex','String','EM IN','Color','b','FontSize',26);
            text('Position',[20.0 0],'Interpreter','latex','String','EM FROM','Color',[.4 .4 .4],'FontSize',26);
            tn = regexprep(tn,' ','_');
            saveas(gcf,strcat(pfold,tn,'.jpg'))
            close all
        end

        
        fprintf('Processed Emissions for %s year %i\n',char(comps(com)),Tyear)
        save(fileout,'NV','TD','L_IN','H_IN','B_IN','L_FROM','H_FROM','B_FROM','ORdComp','FRdComp','EF_IN','EF_FROM','EM_IN','EM_FROM')
        close all
    end
end
addpath('/storage/nilu/Inby/Emission_Group/Emission_Models/HEDGE/MiljodirektoratetData')
NERVE_Miljodirektoratet_Output_generator()
end