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
function [Sn] = Emissions_Calculations_SSB()
%--------------------------------------------------------------------------
%
% 22.10.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global tfold Tyear SSB_Vehicle_dist comps RLinks Vehicle_dist Vehicle_weight
global debug_mode

%--------------------------------------------------------------
% MUNICIPALITY EMISSION FACTORS FOR ALL ROADS
%--------------------------------------------------------------



%--------------------------------------------------------------
% ROAD LINK
%--------------------------------------------------------------
% extract fields needed for calculations
fprintf('\nExtracts necessary fields   \n')
fprintf('Traffic and Vehicle Year:  %i \n',Tyear)
fprintf('RoadLinks               :  %i \n',size(RLinks,1))

file = sprintf('%s%s',tfold,'roads');
fprintf('Warning,using roads file not produced by NERVE model\n%s\n',file)
load(file)

TM = readtable(SSB_Vehicle_dist,'Sheet','MODEL');
LightVehiclesIdx = TM.ClassNum==1|TM.ClassNum==2;
BusesVehiclesIdx = TM.ClassNum==3|TM.ClassNum==4;
HeavyVehiclesIdx = TM.ClassNum==5|TM.ClassNum==6|TM.ClassNum==7;




% HBEFA = extractfield(RLinks,'HBEFA_EQIV');                % Type of Road        (Access->MW)
% SPD   = extractfield(RLinks,'SPEED');                     % Speed on Road       (km hr-1)
% DEC   = extractfield(RLinks,'SLOPE');                     % Verticality of Road (%)
% ENV   = extractfield(RLinks,'URBAN');                     % Urbanity of Road    (URB/RUR)
% LEN   = extractfield(RLinks,'DISTANCE')*1000;             % Length of Road      (in metres)
% RU    = extractfield(RLinks,'RUSH_DELAY');

dayInYear = round(datenum([Tyear+1, 1, 1, 12 0 0 ])-datenum([Tyear, 1, 1, 12, 0 ,0]));
fprintf('Found %i days in year: %i\n',dayInYear,Tyear)

% Make a cleaned version for writing to file
Tout = struct2table(RLinks);


% roads
KommunerS = extractfield(RLinks,'KOMMS');
KommunerE = extractfield(RLinks,'KOMME');
uKomm    = unique(KommunerS);

% Congestion weights.  Traffic part delayed (TPD)
% Combine the different congestion levels by volume to determine the
% actual EF:
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



for com = 1:length(comps)
    fprintf('<--- %s ---\n',char(comps(com)))

    TEF = readtable('OnRoadEF_RoadClasses.xlsx','Sheet',sprintf('%s_%i',char(comps(com)),Tyear),'PreserveVariableNames',1);    
    % load(oEFfile,'TFout','roads');
    
    ttef = TEF(:,1:7);
    
    % All road links emissions Light / Heavy / Bus
    nEMISS_L = zeros(size(RLinks));
    nEMISS_H = zeros(size(RLinks));
    nEMISS_B = zeros(size(RLinks));

    
    % All road links kilometer per type of road Light / Heavy / Bus
    nDD_L = zeros(height(ttef));
    nDD_H = zeros(height(ttef));
    nDD_B = zeros(height(ttef));    
    
    for komm = 1:length(uKomm)
        % extract all necessary roadLinks:
        f1 = find(KommunerS == uKomm(komm));
        f2 = find(KommunerE == uKomm(komm));
        f  = unique([f1,f2]);
        
        MRLinks = RLinks(f);        
        fprintf('\n-%3i Kommune_%04i --- \n',komm,Vehicle_dist.D1_KommNr(komm))
        fprintf('     Has %i Roads     \n',length(f))
        
        % extract all roadEFs:
        ref = find(contains(TEF.Properties.VariableNames,{sprintf('%i',uKomm(komm))}));
        Tef = [ttef,TEF(:,ref)];
        
        Lef = find(ismember(Tef.Properties.VariableNames,sprintf('EF_Light_%04i',Vehicle_dist.D1_KommNr(komm))));
        Hef = find(ismember(Tef.Properties.VariableNames,sprintf('EF_Heavy_%04i',Vehicle_dist.D1_KommNr(komm))));
        Bef = find(ismember(Tef.Properties.VariableNames,sprintf('EF_Buses_%04i',Vehicle_dist.D1_KommNr(komm))));

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
            fprintf('     Light Traffic L=%7.1f (1 000 000) Km  TDL=%7.1f   \n',1e-6*sum(L.*IDO.*LEN)*dayInYear,1e-6*sum(Vehicle_dist.modelTD(komm,LightVehiclesIdx)))
            fprintf('     Heavy Traffic H=%7.1f (1 000 000) Km  TDH=%7.1f   \n',1e-6*sum(H.*IDO.*LEN)*dayInYear,1e-6*sum(Vehicle_dist.modelTD(komm,HeavyVehiclesIdx)))
            fprintf('     Buses Traffic B=%7.1f (1 000 000) Km  TDB=%7.1f   \n',1e-6*sum(B.*IDO.*LEN)*dayInYear,1e-6*sum(Vehicle_dist.modelTD(komm,BusesVehiclesIdx)))
        end
        
        % *my is the total distance driven on each link in Tyear
        Lmy = zeros(size(Tef,1),1);
        Hmy = zeros(size(Tef,1),1);
        Bmy = zeros(size(Tef,1),1);
        
        % EMISS_* are the road link emissions for each link in the
        % municipality
        EMISS_L = zeros(size(MRLinks));
        EMISS_H = zeros(size(MRLinks));
        EMISS_B = zeros(size(MRLinks));
        
        for r =1:length(L)
            dec  = find(roads.RoadGradient == DEC(r));
            env  = find(roads.RoadEnvID	  == ENV(r)+1);
            SPD(r) = round(SPD(r)/10)*10;
            if SPD(r) < 30;  SPD(r)=30;  end
            if SPD(r) > 110; SPD(r)=110; end
            spd  = find(roads.RoadSpeeds  == SPD(r));
            name = {sprintf('%s/%s-%i/%i%%/%s',char(roads.RoadEnv(env)),char(roads.RoadType(HBEFA(r))),roads.RoadSpeeds(spd),roads.RoadGradient(dec))};
            
            idx = find(contains(Tef.Name,name));
            if KOMS(r) == uKomm(komm)
                Lmy(idx) = Lmy(idx) + dayInYear *IDO(r)* L(r) *LEN(r)*TPD_L(RU(r)+1,:)';
                Hmy(idx) = Hmy(idx) + dayInYear *IDO(r)* H(r) *LEN(r)*TPD_H(RU(r)+1,:)';
                Bmy(idx) = Bmy(idx) + dayInYear *IDO(r)* B(r) *LEN(r)*TPD_B(RU(r)+1,:)';
            else
                Lmy(idx) = Lmy(idx) + dayInYear *IDO(r)* L(r) *LEN(r)*TPD_L(RU(r)+1,:)';
                Hmy(idx) = Hmy(idx) + dayInYear *IDO(r)* H(r) *LEN(r)*TPD_H(RU(r)+1,:)';
                Bmy(idx) = Bmy(idx) + dayInYear *IDO(r)* B(r) *LEN(r)*TPD_B(RU(r)+1,:)';
            end
            EMISS_L(r) = sum(Lmy(idx).*table2array(Tef(idx,Lef)));
            EMISS_H(r) = sum(Hmy(idx).*table2array(Tef(idx,Hef)));
            EMISS_B(r) = sum(Bmy(idx).*table2array(Tef(idx,Bef)));
        end
        tEML = Lmy.* table2array(Tef(:,Lef));
        tEMH = Hmy.* table2array(Tef(:,Hef));
        tEMB = Bmy.* table2array(Tef(:,Bef));
        
        nEMISS_L(f) = nEMISS_L(f) + EMISS_L;
        nEMISS_H(f) = nEMISS_H(f) + EMISS_H;
        nEMISS_B(f) = nEMISS_B(f) + EMISS_B;
        
 
        if debug_mode
            fprintf('---- Lette   %11.1f   Kilometer (%3.0f%%) \n'  ,sum(L),100*sum(Lmy)/sum(Lmy+Hmy+Bmy))
            fprintf('---- Tunge   %11.1f   Kilometer (%3.0f%%) \n'  ,sum(H),100*sum(Hmy)/sum(Lmy+Hmy+Bmy))
            fprintf('---- Busser  %11.1f   Kilometer (%3.0f%%) \n\n',sum(B),100*sum(Bmy)/sum(Lmy+Hmy+Bmy))
            
            fprintf('---- Lette   %11.1f   Gram %s (%3.0f%%)\n'  ,sum(tEML),char(comps(com)),100*sum(tEML)/sum(tEML+tEMH+tEMB))
            fprintf('---- Tunge   %11.1f   Gram %s (%3.0f%%)\n'  ,sum(tEMH),char(comps(com)),100*sum(tEMH)/sum(tEML+tEMH+tEMB))
            fprintf('---- Busser  %11.1f   Gram %s (%3.0f%%)\n\n',sum(tEMB),char(comps(com)),100*sum(tEMB)/sum(tEML+tEMH+tEMB))
            
            fprintf('---- Lette   %11.1f   Gram/Kilometer \n'  ,sum(tEML)/sum(Lmy))
            fprintf('---- Tunge   %11.1f   Gram/Kilometer \n'  ,sum(tEMH)/sum(Hmy))
            fprintf('---- Busser  %11.1f   Gram/Kilometer \n\n',sum(tEMB)/sum(Bmy))
        end
        
        Tef.Light_DD = Lmy;
        Tef.Heavy_DD = Hmy;
        Tef.Buses_DD = Bmy;
        
        Tef.Properties.VariableNames(find(ismember(Tef.Properties.VariableNames,'Light_DD'))) = {sprintf('DD_Light_%04i',Vehicle_dist.D1_KommNr(komm))};
        Tef.Properties.VariableNames(find(ismember(Tef.Properties.VariableNames,'Heavy_DD'))) = {sprintf('DD_Heavy_%04i',Vehicle_dist.D1_KommNr(komm))};
        Tef.Properties.VariableNames(find(ismember(Tef.Properties.VariableNames,'Buses_DD'))) = {sprintf('DD_Buses_%04i',Vehicle_dist.D1_KommNr(komm))};
    end
    fprintf('\n---- NORGE --- \n')
    fprintf('---- Lette   %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*sum(nEMISS_L),char(comps(com)),100*nansum(nEMISS_L)/nansum(nEMISS_L+nEMISS_H+nEMISS_B))
    fprintf('---- Tunge   %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*nansum(nEMISS_H),char(comps(com)),100*nansum(nEMISS_H)/nansum(nEMISS_L+nEMISS_H+nEMISS_B))
    fprintf('---- Busser  %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*nansum(nEMISS_B),char(comps(com)),100*nansum(nEMISS_B)/nansum(nEMISS_L+nEMISS_H+nEMISS_B))
    fprintf('---- Totalt  %11.1f   (1000)Ton %s \n\n',1e-9*nansum(nEMISS_B+nEMISS_H+nEMISS_L),char(comps(com)))
    
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

    TLinks = struct2table(RLinks);
    TLinks.LetEM =  nEMISS_L;
    TLinks.HeaEM =  nEMISS_H;
    TLinks.BusEM =  nEMISS_B;
    
    TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'LetEM'))) = {sprintf('EM_Light_%s',char(comps(com)))};
    TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'HeaEM'))) = {sprintf('EM_Heavy_%s',char(comps(com)))};
    TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'BusEM'))) = {sprintf('EM_Buses_%s',char(comps(com)))};
    RLinks = table2struct(TLinks);
end
Sn = RLinks;
end

% Adj      = SPD/70;
% Km       = LEN*1e-3;
% wearL    = Adj.*MP_wear_Light.*(365*L).*Km;
% wearH    = Adj.*MP_wear_Heavy.*(365*(H+B)).*Km;
% Tout.EM_TW10  = MP_airbornefraction*(wearL + wearH)';
%
% % fileout = sprintf('%sOutput_%i.mat',tfold,Tyear);
% % %save(fileout,'Class_Weight','TPD','aEM_L','aEM_H','aEM_B','aEM_2W','aEM_Tot','nEM_L','nEM_H','nEM_B','nEM_2W','nEM_Tot')
% % save(fileout,'Class_Weight','TPD','aEM_L','aEM_H','aEM_B','aEM_2W','aEM_Tot')
%
%
% % Finally write the output to a shape file.
% fprintf('Re-Structuring Roads...\n')
% Sn = table2struct(Tout);
%end
