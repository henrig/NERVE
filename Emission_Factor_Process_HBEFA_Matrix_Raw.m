function Emission_Factor_Process_HBEFA_Matrix_Raw()
%--------------------------------------------------------------------------

% 22.10.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------

% Script for reading and processing HBEFA .csv file of subsegments.
%
%
% Basically taken from NERVE, but includes some updates for the format of
% HBEFA4.1.
%
% Makes a 6-D matix of the emissions factor, convenient for speedy
% processing.
%--------------------------------------------------------------------------
%
% % Path to emissions
% HBEFA_path = '/storage/nilu/Inby/Emission_Group/Emission_Factors/Traffic/HBEFA_Raw_data/HOT_HBEFA_41/';
%
% % List the avalable species SLA have downloaded Emission factors for
% comps      = [{'NOx'},{'CO2'},{'BC'},{'Be'},{'CO'},{'HC'},{'NMHC'},{'NO2'},{'PM'},{'PN'},{'NH3'},{'N2O'},{'CH4'},{'FC'},{'FC_MJ'}];

global do_preProcessing_HBEFA debug_mode
global HBEFA_path tfold comps

fprintf('* call PreProcess_HBEFA_Matrix_Raw   *\n')
if ~do_preProcessing_HBEFA
    fprintf('We Already Have Necessary HBEFA Emission Factors \n')
    fprintf('... Continuing Without new HBEFA calculations \n ... \n')
    return
else
    fprintf('Remaking HBEFA COMBINED files for Use in HEDGE \n')    
    fprintf('Combine the "HGV" EF file with the "others" files\n')    
end



% cell fields:
SfldList = [{'Case'},{'VehCat'},{'Component'},{'TrafficSit'},{'Gradient'},{'Subsegment'},{'Technology'},{'SizeClasse'},{'EmConcept'},{'x_OfSubsegment'}];
% Numerical fields:
NfldList = [{'Year'},{'TrafficScenario'},{'RoadCat'},{'IDSubsegment'},{'KM'},{'x_OfSubsegment'},{'V'},{'V_0_'},{'V_100_'},...
    {'EFA'},{'EFA_0_'},{'EFA_100_'},{'V_weighted'},{'V_weighted_0_'},{'V_weighted_100_'},{'EFA_weighted'},{'EFA_weighted_0_'},...
    {'EFA_weighted_100_'},{'EFA_WTT'},{'EFA_WTT_0_'},{'EFA_WTT_100_'},{'EFA_WTW'},{'EFA_WTW_0_'},{'EFA_WTW_100_'},{'AmbientCondPattern'}];


redo = 1;
for com = 1:length(comps)
    
    % assumed a naming convention
    ifile1  = sprintf('%sEFA_HOT_Subsegm_%s.csv',HBEFA_path,char(comps(com)));
    ifile2  = sprintf('%sHGV/EFA_HOT_Subsegm_%s_hgv.csv',HBEFA_path,char(comps(com)));
    
    % avoid doing process repeatedly (save a combined file)
    ifile3  = sprintf('%sCOMBINED_EFA_HOT_Subsegm_%s.csv',HBEFA_path,char(comps(com)));
    
    if ~exist(ifile3) || redo ==1
        fprintf('ifile1: %s\n',ifile1)
        Tn1     = readtable(ifile1);
        
        fprintf('ifile2: %s\n',ifile2)
        Tn2     = readtable(ifile2);
        
        % Use Pre-defined lists to find out which fields should be numeric.
        for i = 1:width(Tn1)
            idn = find(ismember(NfldList,Tn1.Properties.VariableNames(i)));
            ids = find(ismember(SfldList,Tn1.Properties.VariableNames(i)));
            if ~isempty(idn)
                ShoulBeNumeric(i) =1;
            else
                ShoulBeNumeric(i) =0;
            end
        end
        ShoulBeNumeric = logical(ShoulBeNumeric);

        
        % Thers a problem in some of the files where columns are read as
        % different types of variablse. Therefore access each table and check
        % if the right columns are numeric.
        numericVars1 = varfun(@isnumeric,Tn1,'output','uniform');
        numericVars2 = varfun(@isnumeric,Tn2,'output','uniform');
        
        % table 1 first
        NT1 =table;
        for i = 1:width(Tn1)
            fprintf('%02i :: %s\n',i,char(Tn1.Properties.VariableNames(i)))
            if ShoulBeNumeric(i)
                if  ~numericVars1(i)
                    NT1(:,i) =  array2table(str2double(table2array(Tn1(:,i))));
                else
                    NT1(:,i) = Tn1(:,i);
                end
            else
                NT1(:,i) = Tn1(:,i);
            end
        end
        NT1.Properties.VariableNames = Tn1.Properties.VariableNames;
        
        % table 2 second
        NT2 =table;
        for i = 1:width(Tn2)
            fprintf('%02i :: %s\n',i,char(Tn2.Properties.VariableNames(i)))
            if ShoulBeNumeric(i)
                if  ~numericVars2(i)
                    NT2(:,i) =  array2table(str2double(table2array(Tn2(:,i))));
                else
                    NT2(:,i) = Tn2(:,i);
                end
            else
                NT2(:,i) = Tn2(:,i);
            end
        end
        NT2.Properties.VariableNames = Tn2.Properties.VariableNames;
        
        % combine the two tables
        Tn = [NT1;NT2];
        % Test that
        Old = [Tn1(1,:);Tn2(1,:);Tn1(end,:);Tn2(end,:)]
        New = [NT1(1,:);NT2(1,:);NT1(end,:);NT2(end,:)]
        
        writetable(Tn,ifile3)
        clear Tn1 Tn2 NT1 NT2
    else % read already combined file 
        fprintf('Skipping RAW files \n reading already combined file:\n')
        fprintf('ifile3: %s\n',ifile3)
        Tn     = readtable(ifile3);
    end
    %%%%%% FINISHED PAIRING HGV with Others
    %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find the unique subsegments
    expected_veh = 778; % for HBEFA4.1
    sub    = unique(Tn.Subsegment);
    fprintf('Number of subsegments in %s %i\n',char(comps(com)), length(sub))
    if length(sub) ~= expected_veh
        warning('May have flawed input data!')
        fprintf('Normally EXPECTING  %i Vehicles in HBEFA4.1 \n', expected_veh)
    end
    
    % Gradients are cell arrays 0%, 2% 4% etc; make them numeric.
    fprintf('Making gradients numeric ...\n')
    S         = table2struct(Tn);
    pct       = extractfield(S,'Gradient')';
    b         = regexp(pct,'\d+(\.)?(\d+)?','match');
    Gradient  = str2double([b{:}]);
    G         = unique(Gradient);
    
    % roads object structure keeps track of the 6-D matrix dimensions.
    roads.RoadGradient = G;
    roads.Dimensions =[{'Subsegment'},{'RoadType'},{'RoadSpeed'},{'Road'},{''},{''}]
    for i=1:length(G)
        index = find(Gradient==G(i));
        roads.RoadGradientID(i) = i;
        fprintf('%i%% Gradient:: %i \n',G(i),length(index))
    end
    fprintf('done...\n')
    clear T S b pct
    
    % Create a table that contains also numeric values. The coulmn title
    % of the each column is converted to the same format as in HBEFA V3.3.
    T=table();
    T.SPECIES       = Tn.Component;
    T.VEHCAT        = Tn.VehCat;
    T.ROAD          = Tn.RoadCat;
    T.ROAD_COND     = Tn.TrafficSit;
    T.CAR           = Tn.Subsegment;
    T.FUEL          = Tn.Technology;
    T.ENGINE        = Tn.SizeClasse;
    T.EUROCAT       = Tn.EmConcept;
    
    % Numerical columns: with silightly simplified names.
    T.GRADIENT      = Gradient';
    
    if  ~isnumeric(Tn.IDSubsegment)
        T.CASEID        = str2double(Tn.IDSubsegment);
    else
        T.CASEID        = Tn.IDSubsegment;
    end
    
    if  ~isnumeric(Tn.IDSubsegment)
        T.ODOMETER      = str2double(Tn.KM);
    else
        T.ODOMETER      = Tn.KM;
    end
    
    if  ~isnumeric(Tn.x_OfSubsegment)
        T.WEIGHT        = str2double(Tn.x_OfSubsegment);
    else
        T.WEIGHT        =  Tn.x_OfSubsegment;
    end
    
    fprintf('Making speeds numeric ...\n')
    if  ~isnumeric(Tn.V)
        T.SPD_AVG       = str2double(Tn.V);
        T.SPD_0         = str2double(Tn.V_0_);
        T.SPD_100       = str2double(Tn.V_100_);
    else
        T.SPD_AVG       = Tn.V;
        T.SPD_0         = Tn.V_0_;
        T.SPD_100       = Tn.V_100_;
    end
    
    fprintf('Making EF numeric ...\n')
    if  ~isnumeric(Tn.EFA)
        T.EFA_AVG       = str2double(Tn.EFA);
        T.EFA_0         = str2double(Tn.EFA_0_);
        T.EFA_100       = str2double(Tn.EFA_100_);
    else
        T.EFA_AVG       = Tn.EFA;
        T.EFA_0         = Tn.EFA_0_;
        T.EFA_100       = Tn.EFA_100_;
    end
    
    if  ~isnumeric(Tn.V_weighted)
        T.WSPD_AVG      = str2double(Tn.V_weighted);
        T.WSPD_0        = str2double(Tn.V_weighted_0_);
        T.WSPD_100      = str2double(Tn.V_weighted_100_);
    else
        T.WSPD_AVG      = Tn.V_weighted;
        T.WSPD_0        = Tn.V_weighted_0_;
        T.WSPD_100      = Tn.V_weighted_100_;
    end
    
    fprintf('Making weighted EF numeric ...\n')
    if ~isnumeric(Tn.EFA_weighted)
        T.WEFA_AVG      = str2double(Tn.EFA_weighted);
        T.WEFA_0        = str2double(Tn.EFA_weighted_0_);
        T.WEFA_100      = str2double(Tn.EFA_weighted_100_);
    else
        T.WEFA_AVG      = Tn.EFA_weighted;
        T.WEFA_0        = Tn.EFA_weighted_0_;
        T.WEFA_100      = Tn.EFA_weighted_100_;
    end
    
    fprintf('Splitting road situations ...\n')
    DC         = (T.ROAD_COND);
    C          = split(DC,'/');
    T.ROAD_ENV = C(:,1);
    
    % Urban or Rural
    fprintf('Setting Environment to ...\n')
    RT = unique(T.ROAD_ENV);
    roads.RoadEnv = RT;
    for i=1:length(RT)
        index               = find(ismember(T.ROAD_ENV, RT(i)));
        roads.RoadEnvID(i)  = i;
        RTID(index,1)       = i;
        fprintf('%i::%s %i \n',i,char(RT(i)),length(index))
    end
    T.ROAD_ENVID = RTID;
    
    % From local access road to Motorway
    fprintf('Setting RoadType to ...\n')
    T.ROAD_TYPE    = join(C(:,1:3),'/');
    T.ROAD_TYPE    = C(:,2);
    RT             = unique(T.ROAD_TYPE);
    roads.RoadType = RT;
    for i=1:length(RT)
        index               = find(ismember(T.ROAD_TYPE, RT(i)));
        roads.RoadTypeID(i) = i;
        RTID(index,1)       = i;
        fprintf('%i::%14s %i \n',i,char(RT(i)),length(index))
    end
    T.ROAD_TYPEID  = RTID;
    
    % Speedlimit of the road, contains some non numeric entries, clen them
    % up.
    fprintf('Finding Road Speedlimits ...\n')
    T.ROAD_SPEED   = str2double(C(:,3));
    
    f = find(isnan(T.ROAD_SPEED));
    if size(f,1)>0
        n               = split(C(f,3),'>');
        T.ROAD_SPEED(f) = str2double(n(:,2));
    end
    
    
    Sp=unique(T.ROAD_SPEED);
    roads.RoadSpeeds = Sp;
    for i=1:length(Sp)
        index                  = find(T.ROAD_SPEED==Sp(i));
        roads.RoadSpeedsID(i)  = i;
        fprintf('%i::%3i %i \n',i,Sp(i),length(index))
    end
    
    % added a 5th traffic situation
    fprintf('Finding Congestion Levels ...\n')
    % Congestion level at the time of driving on the road
    Congestion       = unique(C(:,4));
    roads.Congestion = Congestion;
    for i=1:length(Congestion)
        index = find(ismember(C(:,4), Congestion(i)));
        if i==1
            CONG(index,1)=1;
        elseif i==2
            CONG(index,1)=3;
        elseif i==3
            CONG(index,1)=2;
        elseif i==4
            CONG(index,1)=4;
        elseif i==5
            CONG(index,1)=5;
        end
        roads.CongestionID(i) = i;
    end
    T.ROAD_CONG=CONG;
    clear CONG
    Co=unique(T.ROAD_CONG);
    
    
    roads.Vehicles          = unique(T.CAR);
    roads.VehiclesID        = 1:length(roads.Vehicles);
    roads.EmissionComponent = comps(com);
    
    
    veh = roads.VehiclesID;
    roa = roads.RoadTypeID;
    spd = roads.RoadSpeeds;
    gra = roads.RoadGradient;
    cog = roads.CongestionID;
    env = roads.RoadEnvID;
    
    TCAR  = table;
    
    for v = 1:length(veh)
        Tt = T(ismember(T.CAR ,roads.Vehicles(v)),:);
        TCAR(1,:) = Tt(1,:);
        count(1:height(Tt)) = 0;
        for r = 1:length(roa)
            iroa = Tt.ROAD_TYPEID == roa(r);
            for s = 1:length(spd)
                ispd = Tt.ROAD_SPEED  == spd(s);
                for g = 1:length(gra)
                    igra = Tt.GRADIENT    == gra(g);
                    for c = 1:length(cog)
                        icog = Tt.ROAD_CONG   == cog(c);
                        for e = 1:length(env)
                            ienv = Tt.ROAD_ENVID  == env(e);
                            I=iroa&ispd&igra&icog&ienv;
                            
                            if sum(I)==1
                                if ~isnumeric(Tt.EFA_AVG(I)); stop; end
                                
                                EF_AVG(v,r,s,g,c,e) = Tt.EFA_AVG(I);
                                EF_000(v,r,s,g,c,e) = Tt.EFA_0(I);
                                EF_100(v,r,s,g,c,e) = Tt.EFA_100(I);
                                
                                if isnan(Tt.EFA_AVG(I))
                                    
                                end
                                count(I) = 1;
                                
                            elseif sum(I)>1
                                if spd(s) ==130 && sum(I)==2 % special case
                                    EF_AVG(v,r,s,g,c,e) = mean(Tt.EFA_AVG(I));
                                    EF_000(v,r,s,g,c,e) = mean(Tt.EFA_0(I));
                                    EF_100(v,r,s,g,c,e) = mean(Tt.EFA_100(I));
                                    count(I) = 1;
                                else
                                    warning('oops, multiply defined EF')
                                    Tt(I,:)
                                end
                            else
                                EF_AVG(v,r,s,g,c,e) = NaN;
                                EF_000(v,r,s,g,c,e) = NaN;
                                EF_100(v,r,s,g,c,e) = NaN;
                                
                            end
                        end
                    end
                end
            end
        end
        Driving(v)   = sum(sum(sum(sum(sum(~isnan(EF_AVG(v,:,:,:,:,:)))))));
        DriveOrig(v) = sum(~isnan(Tt.EFA_AVG));
        
        if debug_mode
        fprintf('%03i %36s %i / %i :::: ',v, char(roads.Vehicles(v)),sum(count),height(Tt))
        fprintf('Numeric EF %i / %i\n',Driving(v),sum(~isnan(Tt.EFA_AVG)))
        end
    end
    
    % save(sprintf('EFA_matrix_RAW_%s',char(comps(com))),'roads','EF_AVG','EF_000','EF_100')
    oEFfile = sprintf('%s/EFA_matrix41_RAW_%s',tfold,char(comps(com)));
    save(oEFfile,'roads','EF_AVG','EF_000','EF_100');
    fprintf('Saved temporary file\n%s\n',oEFfile)
    
    clear Gradient Tn T RTID roads EF_*
end

end
