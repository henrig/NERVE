function [Sn] = Emissions_Calculations_HBEFA(Calc_Links)
global Tyear comps tfold Vehicle_source

% extract fields needed for calculations
fprintf('\nExtracts necessary fields   \n')
fprintf('Traffic and Vehicle Year:  %i \n',Tyear)
fprintf('RoadLinks               :  %i \n',size(Calc_Links,1))

L     = extractfield(Calc_Links,sprintf('L_ADT%04i',Tyear));  % Traffic Volume (# day-1)
H     = extractfield(Calc_Links,sprintf('H_ADT%04i',Tyear));  % Traffic Volume (# day-1)
B     = extractfield(Calc_Links,sprintf('B_ADT%04i',Tyear));  % Traffic Volume (# day-1)
HBEFA = extractfield(Calc_Links,'HBEFA_EQIV'); % Type of Road        (Access->MW)
SPD   = extractfield(Calc_Links,'SPEED');      % Speed on Road       (km hr-1)
DEC   = extractfield(Calc_Links,'SLOPE');      % Verticality of Road (%)
ENV   = extractfield(Calc_Links,'URBAN');      % Urbanity of Road    (URB/RUR)
LEN   = extractfield(Calc_Links,'DISTANCE')*1000; % Length of Road      (in metres)
RU    = extractfield(Calc_Links,'RUSH_DELAY');

% Classification of rush hour traffic. 0-4 (4 is congested)
ru=zeros(size(RU));
ru(RU<=0.1)        = 1;
ru(RU>0.1 & RU<=4) = 2;
ru(RU>4 & RU<=15)  = 3;
ru(RU>15)          = 4;

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



    % Internal weight for traffic type going from Class -> Traffic type
    Class_Weight = [1, 0.8, 0.2, 0.5, 0.5, 0.5, 0.5];
end


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
    
    if ismember(Vehicle_source,{'HBEFA'}
        % load the emission factor
        EFfile = sprintf('%sHBEFA41_%i_PROCESSED_%s',tfold,Tyear,char(comps(com)));
        load(EFfile)
        fprintf('Loaded:\n%s\n',EFfile)
        fprintf('******* %s **********\n',char(comps(com)))
        fprintf('Loaded: %s Vehicles     : %i \n',char(comps(com)), size(EF_AVG,1))
        fprintf('Loaded: %s Classes      : %i \n',char(comps(com)), size(EFCLASS,1))
    end
    
    % loop through Roads individually
    for f = 1:length(Calc_Links)
        
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
        
        % Use the Indexes to find (for each Vehicle) The EF to apply.
        ROAD_EF_CLASS(1:size(EFCLASS,1),1:size(EFCLASS,5)) = EFCLASS(:,HBEFA(f),SpIdx(f),SlIdx(f),:,EnIdx(f));
        
        % need to do some replacement stuff for the urban roads missing EF for
        % TT/AT and HDV
        if isnan(ROAD_EF_CLASS(end-1,1))
            if EnIdx(f) == 2
                ROAD_EF_CLASS(6:7,1:size(EFCLASS,5)) = EFCLASS(6:7,HBEFA(f),SpIdx(f),SlIdx(f),:,1);
            end
        end
        

        CEF = zeros(size(EFCLASS,1),1);
        for cl = 1:size(ROAD_EF_CLASS,1)
            CEF(cl) = sum(TPD(ru(f),:).*ROAD_EF_CLASS(cl,:),2);
        end
        
        nCEF = zeros(size(EFCLASS,1),1);
        for cl = 1:size(ROAD_EF_CLASS,1)
            nCEF(cl) = sum([1,0,0,0,0].*ROAD_EF_CLASS(cl,:),2);
        end
        
        
        EM_2W(f)     = 0   *365*LEN(f)*1e-3*CEF(1) *Class_Weight(1);
        EM_L(f)      = L(f)*365*LEN(f)*1e-3*(CEF(2)*Class_Weight(2) + CEF(3)*Class_Weight(3));
        EM_B(f)      = B(f)*365*LEN(f)*1e-3*(CEF(4)*Class_Weight(4) + CEF(5)*Class_Weight(5));
        EM_H(f)      = H(f)*365*LEN(f)*1e-3*(CEF(6)*Class_Weight(6) + CEF(7)*Class_Weight(7));
        EM_Tot(f)    = EM_2W(f) + EM_L(f) + EM_B(f) + EM_H(f);
        
%         % Repeat calculations without congestion
%         n1EM_2W(f)     = 0   *365*LEN(f)*1e-3*CnEF(1) *Class_Weight(1);
%         n1EM_L(f)      = L(f)*365*LEN(f)*1e-3*(nCEF(2)*Class_Weight(2) + nCEF(3)*Class_Weight(3));
%         n1EM_B(f)      = B(f)*365*LEN(f)*1e-3*(nCEF(4)*Class_Weight(4) + nCEF(5)*Class_Weight(5));
%         n1EM_H(f)      = H(f)*365*LEN(f)*1e-3*(nCEF(6)*Class_Weight(6) + nCEF(7)*Class_Weight(7));
%         n1EM_Tot(f)    = EM_2W(f) + EM_L(f) + EM_B(f) + EM_H(f);
        
        if isnan(EM_Tot(f))
            fprintf('road %i No emissions  ',f)
            fprintf('HBEFA %i; ',HBEFA(f))
            fprintf('SpIdx %i; ',SpIdx(f))
            fprintf('SlIdx %i; ',SlIdx(f))
            fprintf('EnIdx %i\n',EnIdx(f))
            fprintf(' %i\n',EnIdx(f))
            fprintf('%04i %s/%s/%i/%i%%\n',f,char(roads.RoadEnv(EnIdx(f))),...
                char(roads.RoadType(HBEFA(f))),...
                roads.RoadSpeeds(SpIdx(f)),...
                roads.RoadGradient(SlIdx(f)))
        end
        %if rem(f,100000)==0; fprintf('Road = %i of %i\n',f, size(Calc_Links,1)); end
    end % road link loop.
    aEM_L(com,:)   = EM_L;
    aEM_B(com,:)   = EM_B;
    aEM_H(com,:)   = EM_H;
    aEM_2W(com,:)  = EM_2W;
    aEM_Tot(com,:) = EM_Tot;

%     nEM_L(com,:)   = n1EM_L;
%     nEM_B(com,:)   = n1EM_B;
%     nEM_H(com,:)   = n1EM_H;
%     nEM_2W(com,:)  = n1EM_2W;
%     nEM_Tot(com,:) = n1EM_Tot;
    
    % Assign the output table an rename the variable added to contain
    % component information.
    Tout.EM_Tot = EM_Tot';
    Tout.Properties.VariableNames(end) = {sprintf('EM_%s',char(comps(com)))};    
    fprintf('_________________________________________\n')
    fprintf('Light emissions of %s  : %8.2f Ton\n',char(comps(com)),sum(EM_L)*1e-6)
    fprintf('Buses emissions of %s  : %8.2f Ton\n',char(comps(com)),sum(EM_B)*1e-6)
    fprintf('Heavy emissions of %s  : %8.2f Ton\n',char(comps(com)),sum(EM_H)*1e-6)
    fprintf('_________________________________________\n')
    fprintf('Total emissions of %s  : %8.2f Ton\n',char(comps(com)),sum(EM_Tot)*1e-6)
    fprintf('==========================================\n')
    
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

fileout = sprintf('%sOutput_%i.mat',tfold,Tyear);
%save(fileout,'Class_Weight','TPD','aEM_L','aEM_H','aEM_B','aEM_2W','aEM_Tot','nEM_L','nEM_H','nEM_B','nEM_2W','nEM_Tot')
save(fileout,'Class_Weight','TPD','aEM_L','aEM_H','aEM_B','aEM_2W','aEM_Tot')


% Finally write the output to a shape file.
fprintf('Re-Structuring Roads...\n')
Sn = table2struct(Tout);
end
