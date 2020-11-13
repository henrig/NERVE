
function Emission_Factors_Road_DrivingDistance_IN_Municipalities()
global RLinks tfold Tyear Vehicle_dist debug_mode SSB_Vehicle_dist ofiles


TM = readtable(SSB_Vehicle_dist,'Sheet','MODEL');

LightVehiclesIdx = TM.ClassNum==1|TM.ClassNum==2;
BusesVehiclesIdx = TM.ClassNum==3|TM.ClassNum==4;
HeavyVehiclesIdx = TM.ClassNum==5|TM.ClassNum==6|TM.ClassNum==7;


T = readtable('modelHBEFA.xlsx','Sheet','Roads');
file = sprintf('%s%s',tfold,'roads');
fprintf('### Warning,using roads file not produced by NERVE model\n%s\n',file)
load(file)


% All road links kilometer per type of road Light / Heavy / Bus
nDD_L = zeros(height(T),1);
nDD_H = zeros(height(T),1);
nDD_B = zeros(height(T),1);

%--------------------------------------------------------------
% ROAD LINK Calculations
%--------------------------------------------------------------
% extract fields needed for calculations
fprintf('\nExtracts necessary fields   \n')
fprintf('Traffic and Vehicle Year:  %i \n',Tyear)
fprintf('RoadLinks               :  %i \n',size(RLinks,1))
dayInYear = round(datenum([Tyear+1, 1, 1, 12 0 0 ])-datenum([Tyear, 1, 1, 12, 0 ,0]));
fprintf('Found %3i days in year  :  %i \n',dayInYear,Tyear)


KommunerS = extractfield(RLinks,'KOMMS');
KommunerE = extractfield(RLinks,'KOMME');
Kommuner  = extractfield(RLinks,'KOMM');
uKomm     = unique(KommunerS);

% Congestion weights.  Traffic part delayed (TPD)
% Combine the different congestion levels by volume to determine the
% actual EF:
fprintf('### Warning, using parameterized CONGESTION\n%s\n',file)
TPD_L=[1.0, 0.0, 0.0,  0.0,  0.0;...
    0.6, 0.4, 0.0,  0.0,  0.0;...
    0.6, 0.2, 0.2,  0.0,  0.0;...
    0.5, 0.2, 0.15, 0.15, 0.0];

TPD_H=[1.0, 0.0, 0.0,  0.0,  0.0;...
    0.6, 0.4, 0.0,  0.0,  0.0;...
    0.6, 0.2, 0.2,  0.0,  0.0;...
    0.5, 0.2, 0.15, 0.15, 0.0];

TPD_B=[1.0, 0.0, 0.0,  0.0,  0.0;...
    0.6, 0.4, 0.0,  0.0,  0.0;...
    0.6, 0.2, 0.2,  0.0,  0.0;...
    0.5, 0.2, 0.15, 0.15, 0.0];
%-----------------------------------
TPD  =[1.0, 0.0, 0.0,  0.0,  0.0;...
    0.6, 0.4, 0.0,  0.0,  0.0;...
    0.6, 0.2, 0.2,  0.0,  0.0;...
    0.5, 0.2, 0.15, 0.15, 0.0];

Tout = T;
for komm = 1:length(uKomm)
    % extract all necessary roadLinks:
    f1 = find(KommunerS == uKomm(komm));
    f2 = find(KommunerE == uKomm(komm));
    f3 = find(Kommuner  == uKomm(komm));
    f  = unique([f1,f2,f3]);
    MRLinks = RLinks(f);
    fprintf('\n-%3i Kommune_%04i --- \n',komm,Vehicle_dist.D1_KommNr(komm))
    fprintf('     Has %i (%i) Roads     \n',length(f1),length(f))
    
    L       = extractfield(MRLinks,sprintf('L_ADT%04i',Tyear));  % Traffic Volume (# day-1)
    H       = extractfield(MRLinks,sprintf('H_ADT%04i',Tyear));  % Traffic Volume (# day-1)
    B       = extractfield(MRLinks,sprintf('B_ADT%04i',Tyear));  % Traffic Volume (# day-1)
    HBEFA   = extractfield(MRLinks,'HBEFA_EQIV');                % Type of Road        (Access->MW)
    SPD     = extractfield(MRLinks,'SPEED');                     % Speed on Road       (km hr-1)
    DEC     = extractfield(MRLinks,'SLOPE');                     % Verticality of Road (%)
    ENV     = extractfield(MRLinks,'URBAN');                     % Urbanity of Road    (URB/RUR)
    LEN     = extractfield(MRLinks,'DISTANCE');                  % Length of Road      (in kilometres)
    RU      = extractfield(MRLinks,'RUSH_DELAY');
    
    % Classification of rush hour traffic. 0-4 (4 is congested)
    ru                 = zeros(size(RU));
    ru(RU<=0.1)        = 1;
    ru(RU>0.1 & RU<=4) = 2;
    ru(RU>4 & RU<=15)  = 3;
    ru(RU>15)          = 4;
    
    IDO     = extractfield(MRLinks,'IDO');
    KOMS    = extractfield(MRLinks,'KOMMS');
    
    if debug_mode
        fprintf('     Light Traffic L=%7.1f (1 000 000) Km  TDL=%9.3f   \n',1e-6*sum(L.*IDO.*LEN)*dayInYear,1e-6*sum(Vehicle_dist.modelTD(komm,LightVehiclesIdx)))
        fprintf('     Heavy Traffic H=%7.1f (1 000 000) Km  TDH=%9.3f   \n',1e-6*sum(H.*IDO.*LEN)*dayInYear,1e-6*sum(Vehicle_dist.modelTD(komm,HeavyVehiclesIdx)))
        fprintf('     Buses Traffic B=%7.1f (1 000 000) Km  TDB=%9.3f   \n',1e-6*sum(B.*IDO.*LEN)*dayInYear,1e-6*sum(Vehicle_dist.modelTD(komm,BusesVehiclesIdx)))
    end
    
    % *my is the total distance driven on each ROAD link Type in Tyear
    Lmy = zeros(size(T,1),1);
    Hmy = zeros(size(T,1),1);
    Bmy = zeros(size(T,1),1);
    for r =1:length(L)
        dec  = find(roads.RoadGradient == DEC(r));
        env  = find(roads.RoadEnvID	  == ENV(r)+1);
        SPD(r) = round(SPD(r)/10)*10;
        if SPD(r) < 30;  SPD(r)=30;  end
        if SPD(r) > 110; SPD(r)=110; end
        spd  = find(roads.RoadSpeeds  == SPD(r));
        name = {sprintf('%s/%s-%i/%i%%/%s',char(roads.RoadEnv(env)),char(roads.RoadType(HBEFA(r))),roads.RoadSpeeds(spd),roads.RoadGradient(dec))};
        
        idx = find(contains(T.Name,name));
        if KOMS(r) == uKomm(komm)
            Lmy(idx) = Lmy(idx) + dayInYear *IDO(r)* L(r) *LEN(r)*TPD_L(ru(r),:)';
            Hmy(idx) = Hmy(idx) + dayInYear *IDO(r)* H(r) *LEN(r)*TPD_H(ru(r),:)';
            Bmy(idx) = Bmy(idx) + dayInYear *IDO(r)* B(r) *LEN(r)*TPD_B(ru(r),:)';
        else
            Lmy(idx) = Lmy(idx) + dayInYear *IDO(r)* L(r) *LEN(r)*TPD_L(ru(r),:)';
            Hmy(idx) = Hmy(idx) + dayInYear *IDO(r)* H(r) *LEN(r)*TPD_H(ru(r),:)';
            Bmy(idx) = Bmy(idx) + dayInYear *IDO(r)* B(r) *LEN(r)*TPD_B(ru(r),:)';
        end
    end
    if abs((sum(Lmy)/(sum(L.*IDO.*LEN)*dayInYear))-1)>1e-3
       fprintf('###  more than 1%% shift %f \n',sum(Lmy)/(sum(L.*IDO.*LEN)*dayInYear)) 
    end
    
    Tout.Light_DD = Lmy;
    Tout.Heavy_DD = Hmy;
    Tout.Buses_DD = Bmy;
    Tout.Properties.VariableNames(find(ismember(Tout.Properties.VariableNames,'Light_DD'))) = {sprintf('DD_Light_%04i',Vehicle_dist.D1_KommNr(komm))};
    Tout.Properties.VariableNames(find(ismember(Tout.Properties.VariableNames,'Heavy_DD'))) = {sprintf('DD_Heavy_%04i',Vehicle_dist.D1_KommNr(komm))};
    Tout.Properties.VariableNames(find(ismember(Tout.Properties.VariableNames,'Buses_DD'))) = {sprintf('DD_Buses_%04i',Vehicle_dist.D1_KommNr(komm))};
end

fprintf('\n---- NORGE --- \n')
L     = extractfield(RLinks,sprintf('L_ADT%04i',Tyear));  % Traffic Volume (# day-1)
H     = extractfield(RLinks,sprintf('H_ADT%04i',Tyear));  % Traffic Volume (# day-1)
B     = extractfield(RLinks,sprintf('B_ADT%04i',Tyear));  % Traffic Volume (# day-1)
LEN   = extractfield(RLinks,'DISTANCE');
LW  = 1e-6*sum(L.*LEN)*dayInYear;
HW  = 1e-6*sum(H.*LEN)*dayInYear;
BW  = 1e-6*sum(B.*LEN)*dayInYear;
LTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,LightVehiclesIdx)));
HTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,HeavyVehiclesIdx)));
BTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,BusesVehiclesIdx)));

fprintf('     Light Traffic L=%7.1f (1 000 000) Km  TDL=%7.1f  (%5.1f%%) \n',LW,LTD,100*LW/LTD)
fprintf('     Heavy Traffic H=%7.1f (1 000 000) Km  TDH=%7.1f  (%5.1f%%) \n',HW,HTD,100*HW/HTD)
fprintf('     Buses Traffic B=%7.1f (1 000 000) Km  TDB=%7.1f  (%5.1f%%) \n',BW,BTD,100*BW/BTD)



writetable(Tout,'Municipal_DrivingDistances_per_RoadType.xlsx','Sheet',sprintf('DD_%i',Tyear))

DrivingDistances_per_RoadType = Tout;
save(ofiles.MatlabOutput,'DrivingDistances_per_RoadType','-append')

end

%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%     % EMISS_* are the road link emissions for each link in the
%     % municipality
%     EMISS_L = zeros(size(MRLinks));
%     EMISS_H = zeros(size(MRLinks));
%     EMISS_B = zeros(size(MRLinks));
%
%     for r =1:length(L)
%         dec  = find(roads.RoadGradient == DEC(r));
%         env  = find(roads.RoadEnvID	  == ENV(r)+1);
%         SPD(r) = round(SPD(r)/10)*10;
%         if SPD(r) < 30;  SPD(r)=30;  end
%         if SPD(r) > 110; SPD(r)=110; end
%         spd  = find(roads.RoadSpeeds  == SPD(r));
%         name = {sprintf('%s/%s-%i/%i%%/%s',char(roads.RoadEnv(env)),char(roads.RoadType(HBEFA(r))),roads.RoadSpeeds(spd),roads.RoadGradient(dec))};
%
%         idx = find(contains(Tef.Name,name));
%         if KOMS(r) == uKomm(komm)
%             Lmy(idx) = Lmy(idx) + dayInYear *IDO(r)* L(r) *LEN(r)*TPD_L(ru(r),:)';
%             Hmy(idx) = Hmy(idx) + dayInYear *IDO(r)* H(r) *LEN(r)*TPD_H(ru(r),:)';
%             Bmy(idx) = Bmy(idx) + dayInYear *IDO(r)* B(r) *LEN(r)*TPD_B(ru(r),:)';
%         else
%             Lmy(idx) = Lmy(idx) + dayInYear *IDO(r)* L(r) *LEN(r)*TPD_L(ru(r),:)';
%             Hmy(idx) = Hmy(idx) + dayInYear *IDO(r)* H(r) *LEN(r)*TPD_H(ru(r),:)';
%             Bmy(idx) = Bmy(idx) + dayInYear *IDO(r)* B(r) *LEN(r)*TPD_B(ru(r),:)';
%         end
%         %             EMISS_L(r) = sum(Lmy(idx).*table2array(Tef(idx,Lef)))/sum(Lmy(idx));
%         %             EMISS_H(r) = sum(Hmy(idx).*table2array(Tef(idx,Hef)))/sum(Lmy(idx));
%         %             EMISS_B(r) = sum(Bmy(idx).*table2array(Tef(idx,Bef)))/sum(Lmy(idx));
%         EMISS_L(r) = sum(Lmy(idx).*table2array(Tef(idx,Lef)));
%         EMISS_H(r) = sum(Hmy(idx).*table2array(Tef(idx,Hef)));
%         EMISS_B(r) = sum(Bmy(idx).*table2array(Tef(idx,Bef)));
%     end
%
%
%     tEML = Lmy.*table2array(Tef(:,Lef));
%     tEMH = Hmy.*table2array(Tef(:,Hef));
%     tEMB = Bmy.*table2array(Tef(:,Bef));
%
%
%     %         nDD_L = nDD_L + Lmy
%     %         nDD_H = nDD_H + Hmy
%     %         nDD_B = nDD_B + Bmy
%
%     nEMISS_L(f) = nEMISS_L(f) + EMISS_L;
%     nEMISS_H(f) = nEMISS_H(f) + EMISS_H;
%     nEMISS_B(f) = nEMISS_B(f) + EMISS_B;
%
%
%     if debug_mode
%         fprintf('---- Lette   %11.1f   Kilometer (%3.0f%%) \n'  ,sum(L),100*sum(Lmy)/sum(Lmy+Hmy+Bmy))
%         fprintf('---- Tunge   %11.1f   Kilometer (%3.0f%%) \n'  ,sum(H),100*sum(Hmy)/sum(Lmy+Hmy+Bmy))
%         fprintf('---- Busser  %11.1f   Kilometer (%3.0f%%) \n\n',sum(B),100*sum(Bmy)/sum(Lmy+Hmy+Bmy))
%
%         fprintf('---- Lette   %11.1f   Gram %s (%3.0f%%)\n'  ,sum(tEML),char(comps(com)),100*sum(tEML)/sum(tEML+tEMH+tEMB))
%         fprintf('---- Tunge   %11.1f   Gram %s (%3.0f%%)\n'  ,sum(tEMH),char(comps(com)),100*sum(tEMH)/sum(tEML+tEMH+tEMB))
%         fprintf('---- Busser  %11.1f   Gram %s (%3.0f%%)\n\n',sum(tEMB),char(comps(com)),100*sum(tEMB)/sum(tEML+tEMH+tEMB))
%
%         fprintf('---- Lette   %11.1f   Gram/Kilometer \n'  ,sum(tEML)/sum(Lmy))
%         fprintf('---- Tunge   %11.1f   Gram/Kilometer \n'  ,sum(tEMH)/sum(Hmy))
%         fprintf('---- Busser  %11.1f   Gram/Kilometer \n\n',sum(tEMB)/sum(Bmy))
%     end
%
%
%
%     KEF(komm,1)  = sum(Tef(:,Lef).*Tef.Light_DD)/sum(Tef.Light_DD)
%
%
%     Tef.Properties.VariableNames(find(ismember(Tef.Properties.VariableNames,'Light_DD'))) = {sprintf('DD_Light_%04i',Vehicle_dist.D1_KommNr(komm))};
%     Tef.Properties.VariableNames(find(ismember(Tef.Properties.VariableNames,'Heavy_DD'))) = {sprintf('DD_Heavy_%04i',Vehicle_dist.D1_KommNr(komm))};
%     Tef.Properties.VariableNames(find(ismember(Tef.Properties.VariableNames,'Buses_DD'))) = {sprintf('DD_Buses_%04i',Vehicle_dist.D1_KommNr(komm))};
%
%
%
% end
%
%
% fprintf('\n---- NORGE --- \n')
% fprintf('---- Lette   %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*sum(nEMISS_L),char(comps(com)),100*nansum(nEMISS_L)/nansum(nEMISS_L+nEMISS_H+nEMISS_B))
% fprintf('---- Tunge   %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*nansum(nEMISS_H),char(comps(com)),100*nansum(nEMISS_H)/nansum(nEMISS_L+nEMISS_H+nEMISS_B))
% fprintf('---- Busser  %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*nansum(nEMISS_B),char(comps(com)),100*nansum(nEMISS_B)/nansum(nEMISS_L+nEMISS_H+nEMISS_B))
% fprintf('---- Totalt  %11.1f   (1000)Ton %s \n\n',1e-9*nansum(nEMISS_B+nEMISS_H+nEMISS_L),char(comps(com)))
%
% L     = extractfield(RLinks,sprintf('L_ADT%04i',Tyear));  % Traffic Volume (# day-1)
% H     = extractfield(RLinks,sprintf('H_ADT%04i',Tyear));  % Traffic Volume (# day-1)
% B     = extractfield(RLinks,sprintf('B_ADT%04i',Tyear));  % Traffic Volume (# day-1)
% LEN   = extractfield(RLinks,'DISTANCE');
% LW  = 1e-6*sum(L.*LEN)*dayInYear;
% HW  = 1e-6*sum(H.*LEN)*dayInYear;
% BW  = 1e-6*sum(B.*LEN)*dayInYear;
% LTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,LightVehiclesIdx)));
% HTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,HeavyVehiclesIdx)));
% BTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,BusesVehiclesIdx)));
%
% fprintf('     Light Traffic L=%7.1f (1 000 000) Km  TDL=%7.1f  (%5.1f%%) \n',LW,LTD,100*LW/LTD)
% fprintf('     Heavy Traffic H=%7.1f (1 000 000) Km  TDH=%7.1f  (%5.1f%%) \n',HW,HTD,100*HW/HTD)
% fprintf('     Buses Traffic B=%7.1f (1 000 000) Km  TDB=%7.1f  (%5.1f%%) \n',BW,BTD,100*BW/BTD)
% fprintf('     Light Traffic L=%7.1f g/Km\n',sum(nEMISS_L)/LW)
% fprintf('     Heavy Traffic H=%7.1f g/Km\n',sum(nEMISS_H)/HW)
% fprintf('     Buses Traffic B=%7.1f g/Km\n',sum(nEMISS_B)/BW)
%
%
%
%
%
% TLinks = struct2table(RLinks);
% TLinks.LetEM =  nEMISS_L;
% TLinks.HeaEM =  nEMISS_H;
% TLinks.BusEM =  nEMISS_B;
%
% TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'LetEM'))) = {sprintf('EM_Light_%s',char(comps(com)))};
% TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'HeaEM'))) = {sprintf('EM_Heavy_%s',char(comps(com)))};
% TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'BusEM'))) = {sprintf('EM_Buses_%s',char(comps(com)))};
% RLinks = table2struct(TLinks);
% end
