function [Sn] = Emissions_Calculations_SSB(Calc_Links,Vehicle_dist)
%--------------------------------------------------------------------------
%
% 22.10.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global tfold Tyear SSB_Vehicle_dist comps

% 1) Where does variable roads come from?
% 1b) Does it need to be updated in Emission Factor stuff.

% Need Vehicle Distribution for municipalities. 


%--------------------------------------------------------------
% ROAD LINK 
%--------------------------------------------------------------
% extract fields needed for calculations
fprintf('\nExtracts necessary fields   \n')
fprintf('Traffic and Vehicle Year:  %i \n',Tyear)
fprintf('RoadLinks               :  %i \n',size(Calc_Links,1))

file = sprintf('%s%s',tfold,'roads');
fprintf('Warning,using roads file not produced by NERVE model\n%s\n',file)
load(file)

L     = extractfield(Calc_Links,sprintf('L_ADT%04i',Tyear));  % Traffic Volume (# day-1)
H     = extractfield(Calc_Links,sprintf('H_ADT%04i',Tyear));  % Traffic Volume (# day-1)
B     = extractfield(Calc_Links,sprintf('B_ADT%04i',Tyear));  % Traffic Volume (# day-1)
HBEFA = extractfield(Calc_Links,'HBEFA_EQIV');                % Type of Road        (Access->MW)
SPD   = extractfield(Calc_Links,'SPEED');                     % Speed on Road       (km hr-1)
DEC   = extractfield(Calc_Links,'SLOPE');                     % Verticality of Road (%)
ENV   = extractfield(Calc_Links,'URBAN');                     % Urbanity of Road    (URB/RUR)
LEN   = extractfield(Calc_Links,'DISTANCE')*1000;             % Length of Road      (in metres)
RU    = extractfield(Calc_Links,'RUSH_DELAY');

% Classification of rush hour traffic. 0-4 (4 is congested)
ru                 = zeros(size(RU));
ru(RU<=0.1)        = 1;
ru(RU>0.1 & RU<=4) = 2;
ru(RU>4 & RU<=15)  = 3;
ru(RU>15)          = 4;


dayInYear = round(datenum([Tyear+1, 1, 1, 12 0 0 ])-datenum([Tyear, 1, 1, 12, 0 ,0])); 
fprintf('Using %i days in year\n',dayInYear)
% Make a cleaned version for writing to file
Tout = struct2table(Calc_Links);

MP_airbornefraction = 0.1;
MP_wear_Light       = 0.1;
MP_wear_Heavy       = 0.5;

% Congestion weights.  Traffic part delayed (TPD)
% Combine the different congestion levels by volume to determine the
% actual EF:
TPD=[1.0,0.0,0.0,0.0,0.0;...
    0.6,0.4,0.0,0.0,0.0;...
    0.6,0.2,0.2,0.0,0.0;...
    0.6,0.15,0.1,0.1,0.05];



%%%%% Internal weight for traffic type going from Class -> Traffic type
%    Class_Weight = [1, 0.8, 0.2, 0.5, 0.5, 0.5, 0.5];

%--------------------------------------------------------------
% VEHICLE DISTRIBUTION 
%--------------------------------------------------------------
TM = readtable(SSB_Vehicle_dist,'Sheet','MODEL');

LightVehiclesIdx =  TM.ClassNum==1|TM.ClassNum==2;
HeavyVehiclesIdx =  TM.ClassNum==3|TM.ClassNum==4;
BusesVehiclesIdx =  TM.ClassNum==5|TM.ClassNum==6|TM.ClassNum==7;



VD               = Vehicle_dist.Vdist;

%--------------------------------------------------------------
% EMISSION FACTORS 
%--------------------------------------------------------------

aEM_L   = zeros(length(comps),length(Calc_Links));
aEM_B   = zeros(size(aEM_L));
aEM_H   = zeros(size(aEM_L));
aEM_2W  = zeros(size(aEM_L));
aEM_Tot = zeros(size(aEM_L));

nEM_L   = zeros(size(aEM_L));
nEM_B   = zeros(size(aEM_L));
nEM_H   = zeros(size(aEM_L));
nEM_2W  = zeros(size(aEM_L));
nEM_Tot = zeros(size(aEM_L));


% set which compound (com)
for com = 1:length(comps)
    Zroads = 0;
    % load the emission factor
    EFfile = sprintf('%sEFA_matrix41_MODEL_%s.mat',tfold,char(comps(com)));
    load(EFfile)
    fprintf('Loaded:\n%s\n',EFfile)
    fprintf('******* %s **********\n',char(comps(com)))
    fprintf('Loaded: %s Vehicles     : %i \n',char(comps(com)), size(EFA,1))
    
    % loop through Roads individually
    for f = 1:length(Calc_Links)
        komm = extractfield(Calc_Links(f),'KOMM');
        Kidx = find(Vehicle_dist.D1_KommNr == komm);
        Km       = LEN(f)*1e-3;
        % Set the empty streets speedlimit to 40
        % correct whatever have slipped through the original
        % classification of the roads file.
        if rem(SPD(f),10) ~= 0; SPD(f) = round(SPD(f)/10)*10; end
        if SPD(f)< 30; SPD(f)=40; end
        if SPD(f)>130; SPD(f)=40; end
        
        % find the speed index on the road
        SpIdx(f) = find(roads.RoadSpeeds==SPD(f));
        
        % Find the slope index of the road
        SlIdx(f) = find(roads.RoadGradient==DEC(f));
        
        % find the environment
        if ENV(f) == 1
            EnIdx(f) = 2;
        elseif ENV(f) == 0
            EnIdx(f) = 1;
        else
            EnIdx(f) = 2;
            warning(sprintf('No environment set ENV = %i',ENV(f)))
        end
        
        % Fix to update road definitions to what is included in HBEFAv4.1
        if HBEFA(f) == 10 && SPD(f) == 40
            HBEFA(f) = 3;
        end
        
        % calculate the emissions per vehicle (empv) driving the link.
        %   Length(Km) * EF(g/Km) = gram
        for cg = 1:5
            RLEM(cg,:) = VD(Kidx,:)'.*EFA(:,HBEFA(f),SpIdx(f),SlIdx(f),cg,EnIdx(f))*Km;
        end
        
        % Test the average emissions for the different classes.
        Lempv = sum(RLEM(:,LightVehiclesIdx),2);
        Hempv = sum(RLEM(:,HeavyVehiclesIdx),2);
        Bempv = sum(RLEM(:,BusesVehiclesIdx),2); 

        
        % use the rushour integer weight to weight the different emissions.
        LemAve =sum(Lempv.*TPD(ru(f),:)');
        HemAve =sum(Hempv.*TPD(ru(f),:)');
        BemAve =sum(Bempv.*TPD(ru(f),:)');

        EM_L(f) = LemAve *L(f)*dayInYear;
        EM_H(f) = HemAve *L(f)*dayInYear;       
        EM_B(f) = BemAve *L(f)*dayInYear;
        EM_2W(f) = 0;        
        EM_Tot(f)  = EM_L(f)+EM_H(f)+EM_B(f)+EM_2W(f);
                
        if isnan(EM_Tot(f))
            if isnan(EM_L(f))
             fprintf('road %i No emissions  \n',f)                
            end
            Zroads = Zroads+1;
%             fprintf('road %i No emissions  ',f)
%             fprintf('HBEFA %i; ',HBEFA(f))
%             fprintf('SpIdx %i; ',SpIdx(f))
%             fprintf('SlIdx %i; ',SlIdx(f))
%             fprintf('EnIdx %i\n',EnIdx(f))
%             fprintf(' %i\n',EnIdx(f))
%             fprintf('%04i %s/%s/%i/%i%%\n',f,char(roads.RoadEnv(EnIdx(f))),...
%                 char(roads.RoadType(HBEFA(f))),...
%                 roads.RoadSpeeds(SpIdx(f)),...
%                 roads.RoadGradient(SlIdx(f)))
        end
        %if rem(f,100000)==0; fprintf('Road = %i of %i\n',f, size(Calc_Links,1)); end
    end % road link loop.
    
    aEM_L(com,:)   = EM_L;
    aEM_B(com,:)   = EM_B;
    aEM_H(com,:)   = EM_H;
    aEM_2W(com,:)  = EM_2W;
    aEM_Tot(com,:) = EM_Tot;
    
    % Assign the output table an rename the variable added to contain
    % component information.
    Tout.EM_Tot = EM_Tot';
    Tout.Properties.VariableNames(end) = {sprintf('EM_%s',char(comps(com)))};    
    fprintf('_________________________________________\n')
    fprintf('Light emissions of %s  : %8.2f Ton\n',char(comps(com)),nansum(EM_L)*1e-6)
    fprintf('Buses emissions of %s  : %8.2f Ton\n',char(comps(com)),nansum(EM_B)*1e-6)
    fprintf('Heavy emissions of %s  : %8.2f Ton\n',char(comps(com)),nansum(EM_H)*1e-6)
    fprintf('_________________________________________\n')
    fprintf('Total emissions of %s  : %8.2f Ton\n',char(comps(com)),nansum(EM_Tot)*1e-6)
    fprintf('==========================================\n')
    fprintf('ZeroRoads %i\n',Zroads)
    
%     EM_Class(com,1)= sum(EM_L)*1e-6;
%     EM_Class(com,2)= sum(EM_H)*1e-6;
%     EM_Class(com,3)= sum(EM_B)*1e-6;
end % end loop emissions species.
% in the end, add PM10 emissions of microplastics from tyre wear.
fprintf('Calculating Tyre-wear...\n')
Adj      = SPD/70;
Km       = LEN*1e-3;
wearL    = Adj.*MP_wear_Light.*(365*L).*Km;
wearH    = Adj.*MP_wear_Heavy.*(365*(H+B)).*Km;
Tout.EM_TW10  = MP_airbornefraction*(wearL + wearH)';

% fileout = sprintf('%sOutput_%i.mat',tfold,Tyear);
% %save(fileout,'Class_Weight','TPD','aEM_L','aEM_H','aEM_B','aEM_2W','aEM_Tot','nEM_L','nEM_H','nEM_B','nEM_2W','nEM_Tot')
% save(fileout,'Class_Weight','TPD','aEM_L','aEM_H','aEM_B','aEM_2W','aEM_Tot')


% Finally write the output to a shape file.
fprintf('Re-Structuring Roads...\n')
Sn = table2struct(Tout);
end
