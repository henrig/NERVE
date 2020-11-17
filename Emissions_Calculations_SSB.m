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
global debug_mode ofiles input
fprintf('---------------------------------------------------------------\n')
fprintf('in Emissions_Calculations_SSB    *\n')
fprintf('---------------------------------------------------------------\n')

%--------------------------------------------------------------
% ROAD LINK EMISSIONS Calculations
%--------------------------------------------------------------
% extract fields needed for calculations
fprintf('\nExtracts necessary fields   \n')
fprintf('Traffic and Vehicle Year:  %i \n',Tyear)
fprintf('RoadLinks               :  %i \n',size(RLinks,1))

file = sprintf('%s%s',tfold,'roads');
fprintf('### Warning, using roads file not produced by NERVE model\n%s\n',file)
load(file)

TM = readtable(SSB_Vehicle_dist,'Sheet','MODEL');
LightVehiclesIdx = TM.ClassNum==1|TM.ClassNum==2;
BusesVehiclesIdx = TM.ClassNum==3|TM.ClassNum==4;
HeavyVehiclesIdx = TM.ClassNum==5|TM.ClassNum==6|TM.ClassNum==7;
dayInYear = round(datenum([Tyear+1, 1, 1, 12 0 0 ])-datenum([Tyear, 1, 1, 12, 0 ,0]));
fprintf('Found %i days in year: %i\n',dayInYear,Tyear)


% roads
KommunerS = extractfield(RLinks,'KOMMS');
KommunerE = extractfield(RLinks,'KOMME');
uKomm     = unique(KommunerS);

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

% EXTRACT 
L       = extractfield(RLinks,sprintf('L_ADT%04i',Tyear));  % Traffic Volume (# day-1)
H       = extractfield(RLinks,sprintf('H_ADT%04i',Tyear));  % Traffic Volume (# day-1)
B       = extractfield(RLinks,sprintf('B_ADT%04i',Tyear));  % Traffic Volume (# day-1)
LEN     = extractfield(RLinks,'DISTANCE');                  % Length of Road      (in kilometres)

HBEFA   = extractfield(RLinks,'HBEFA_EQIV');                % Type of Road        (Access->MW)
SPD     = extractfield(RLinks,'SPEED');                     % Speed on Road       (km hr-1)
DEC     = extractfield(RLinks,'SLOPE');                     % Verticality of Road (%)
ENV     = extractfield(RLinks,'URBAN');                     % Urbanity of Road    (URB/RUR)
RU      = extractfield(RLinks,'RUSH_DELAY');

% Calculated Properties
LW  = 1e-6*sum(L.*LEN)*dayInYear;
HW  = 1e-6*sum(H.*LEN)*dayInYear;
BW  = 1e-6*sum(B.*LEN)*dayInYear;
LTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,LightVehiclesIdx)));
HTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,HeavyVehiclesIdx)));
BTD = 1e-6*sum(sum(Vehicle_dist.modelTD(:,BusesVehiclesIdx)));

% Classification of rush hour traffic. 0-4 (4 is congested)
ru                 = zeros(size(RU));
ru(RU<=0.1)        = 1;
ru(RU>0.1 & RU<=4) = 2;
ru(RU>4 & RU<=15)  = 3;
ru(RU>15)          = 4;

IDO     = extractfield(RLinks,'IDO');
KOMS    = extractfield(RLinks,'KOMMS');

for com = 1:length(comps)
    % All road links emissions Light / Heavy / Bus
    EMISS_L = zeros(size(RLinks));
    EMISS_H = zeros(size(RLinks));
    EMISS_B = zeros(size(RLinks));
    Link_emission_factor= zeros(size(RLinks));
    fprintf('<--- %s ---\n',char(comps(com)))
    fprintf('Loading large file\n...')
    %%%%%%%%%%%%%%%%%%%%%%%
    % try
    %    TEF = load(ofiles.MatlabOutput,sprintf('OnRoadEF_RoadClasses_%i',char(coms(com))))
    % catch
    %    fprintf('reading excel file\n')
    TEF = readtable(sprintf('OnRoadEF_RoadClasses_%s.xlsx',char(comps(com))),'Sheet',sprintf('%s_%i',char(comps(com)),Tyear),'PreserveVariableNames',1);
    tmpfile = sprintf('%sEF_On_AllRoadCond_Municipality_%i_%s.mat',tfold,Tyear,char(comps(com)));
    load(tmpfile)
    
    % end
    %%%%%%%%%%%%%%%%%%%%%%
    fprintf('Loaded.\n')
    tef = TEF.Name;
    for r =1:length(L)
        % Find the municipality to use Emission factor from:
        % Find the correct emission factor column in TEF:
        komm = extractfield(RLinks(r),'KOMMS');
        Lef  = find(ismember(TEF.Properties.VariableNames,sprintf('EF_Light_%04i',komm)));
        Hef  = find(ismember(TEF.Properties.VariableNames,sprintf('EF_Heavy_%04i',komm)));
        Bef  = find(ismember(TEF.Properties.VariableNames,sprintf('EF_Buses_%04i',komm)));
        
        % Compose the name type of the road:
        dec  = find(roads.RoadGradient == DEC(r));
        env  = find(roads.RoadEnvID	   == ENV(r)+1);
        SPD(r) = round(SPD(r)/10)*10;
        if SPD(r) < 30;  SPD(r)=30;  end
        if SPD(r) > 110; SPD(r)=110; end
        spd  = find(roads.RoadSpeeds  == SPD(r));
        name = {sprintf('%s/%s-%i/%i%%/%s',char(roads.RoadEnv(env)),char(roads.RoadType(HBEFA(r))),roads.RoadSpeeds(spd),roads.RoadGradient(dec))};
        idx = find(contains(tef,name));
        
        % Calculate the annual driving distance on the roadlink in
        % Kilometers:
        DL = L(r) *dayInYear* LEN(r) *IDO(r);
        HL = H(r) *dayInYear* LEN(r) *IDO(r);
        BL = B(r) *dayInYear* LEN(r) *IDO(r);
        
        % Divide into given congestion levels
        Lmy(idx) =   DL*TPD_L(ru(r),:)';
        Hmy(idx) =   HL*TPD_H(ru(r),:)';
        Bmy(idx) =   BL*TPD_B(ru(r),:)';
        
        % Calculate the emissions:
        EMISS_L(r) = sum(Lmy(idx).*table2array(TEF(idx,Lef))');
        EMISS_H(r) = sum(Hmy(idx).*table2array(TEF(idx,Hef))');
        EMISS_B(r) = sum(Bmy(idx).*table2array(TEF(idx,Bef))');
        
        % Calculate roadLink Emission Factor: (Not used as of now)
        Link_Light_emission_factor(r) = EMISS_L(r)/sum(Lmy(idx));
        Link_Heavy_emission_factor(r) = EMISS_H(r)/sum(Hmy(idx));
        Link_Buses_emission_factor(r) = EMISS_B(r)/sum(Bmy(idx));
        Link_emission_factor(r)       = (EMISS_L(r)+EMISS_H(r)+EMISS_B(r))/(sum(Lmy(idx))+sum(Hmy(idx))+sum(Bmy(idx)));
        if rem(r,10000)==0; fprintf('Roads Calculated %i/%i\n  L:%f\n  H:%f\n  B:%f\nEF: %f\n',r,length(RLinks),sum(EMISS_L)*1e-9,sum(EMISS_H)*1e-9,sum(EMISS_B)*1e-9,Link_emission_factor(r)); end
    end
    TLinks = struct2table(RLinks);
    TLinks.LetEM =  EMISS_L;
    TLinks.HeaEM =  EMISS_H;
    TLinks.BusEM =  EMISS_B;
    
    TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'LetEM'))) = {sprintf('EM_L_%s',char(comps(com)))};
    TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'HeaEM'))) = {sprintf('EM_H_%s',char(comps(com)))};
    TLinks.Properties.VariableNames(find(ismember(TLinks.Properties.VariableNames,'BusEM'))) = {sprintf('EM_B_%s',char(comps(com)))};
    RLinks = table2struct(TLinks);
    
    % Print some statistical output
    fprintf('\n---- NORGE --- %i\n',Tyear)
    fprintf('---- Lette   %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*sum(EMISS_L),char(comps(com)),100*nansum(EMISS_L)/nansum(EMISS_L + EMISS_H + EMISS_B))
    fprintf('---- Tunge   %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*nansum(EMISS_H),char(comps(com)),100*nansum(EMISS_H)/nansum(EMISS_L + EMISS_H + EMISS_B))
    fprintf('---- Busser  %11.1f   (1000)Ton %s (%3.0f%%)\n'  ,1e-9*nansum(EMISS_B),char(comps(com)),100*nansum(EMISS_B)/nansum(EMISS_L + EMISS_H + EMISS_B))
    fprintf('---- Totalt  %11.1f   (1000)Ton %s \n\n',1e-9*nansum(EMISS_B+EMISS_H+EMISS_L),char(comps(com)))
    
    
    fprintf('     Light Traffic L=%7.1f g/Km\n',1e-6*sum(EMISS_L)/LW)
    fprintf('     Heavy Traffic H=%7.1f g/Km\n',1e-6*sum(EMISS_H)/HW)
    fprintf('     Buses Traffic B=%7.1f g/Km\n',1e-6*sum(EMISS_B)/BW)
    fprintf('     Light Traffic L=%7.1f (1 000 000) Km  TDL=%7.1f  (%5.1f%%) \n',LW,LTD,100*LW/LTD)
    fprintf('     Heavy Traffic H=%7.1f (1 000 000) Km  TDH=%7.1f  (%5.1f%%) \n',HW,HTD,100*HW/HTD)
    fprintf('     Buses Traffic B=%7.1f (1 000 000) Km  TDB=%7.1f  (%5.1f%%) \n',BW,BTD,100*BW/BTD)

    fprintf('--- %s --->\n',char(comps(com)))
end

save(ofiles.MatlabOutput,'RLinks','-append')
fprintf('Added RLinks to output\n')
Sn     = RLinks;
end
