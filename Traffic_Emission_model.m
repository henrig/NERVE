clear all; close all; clc
%--------------------------------------------------------------------------
% GLOBAL PARAMETERS
% options
global do_preProcessing_HBEFA debug_mode use_temporary_files ofiles tfiles
global Vehicle_source remove_NoTrafficRoads Vehicle_weight Subtract_bio_from_CO2
% folders
global ifold tfold ofold HBEFA_path HBEFA_roads Kartverket
% files
global traffile tettfile Komm_shape Light_traff_years Heavy_traff_years
global Buses_traff_years HBEFA_vehicles SSB_vehicle tfiles traff_exchange
global traff_exchange_sh SSB_Vehicle_dist input Bio_mix_file
% variables
global Ryear Tyear Myear comps Vehicle_dist RLinks

%--------------------------------------------------------------------------
% PATHS
%--------------------------------------------------------------------------
% Toggle Linux/Windows path
if ispc
    bp = 'N:\Inby\Emission_Group\Emission_Models\HEDGE\';
else
    bp =  '/storage/nilu/Inby/Emission_Group/Emission_Models/HEDGE/';
end
% Subfolders within Model
ifold       = sprintf('%s%s',bp,'Input/');
tfold       = sprintf('%s%s',bp,'Temp/');
ofold       = sprintf('%s%s',bp,'OutputNewCong/');
addpath(strcat(pwd,'/Source'))

% Set the year to which traffic volume should be scaled:
for Tyear = 2019:2019
    % Tyear = 2018;
    
    % set the year from which traffic volume is collected:
    Ryear = 2017;
    % path outside Model "sphere"
    HBEFA_path = sprintf('%s%s','/storage/nilu','/Inby/Emission_Group/Emission_Factors/Traffic/HBEFA_Raw_data/HOT_HBEFA_41/');
    %--------------------------------------------------------------------------
    % FILES
    %--------------------------------------------------------------------------
    % set the year from which Municipalities are collected:
    Myear = 2020;
    % All input file names within "Model sphere"
    traffile          = sprintf('%s%s',ifold,'Roads/NORWAY_testFilledRoads');
    % traffile          = sprintf('%s%s',ifold,'Roads/Norway02Dec2019');
    % traffile          = sprintf('%s%s',ifold,'Roads/Norway25Feb2020');
    
    tettfile          = sprintf('%s%s',ifold,'Administrative_Borders/Tettsted2016');
    Komm_shape        = sprintf('%s%s%4i',ifold,'Administrative_Borders/Kommuner',Myear);
    Light_traff_years = sprintf('%s%s',ifold,'Annual_Scaling_Light_2009_2019.csv');
    Heavy_traff_years = sprintf('%s%s',ifold,'Annual_Scaling_Heavy_2009_2019.csv');
    Buses_traff_years = sprintf('%s%s',ifold,'Annual_Scaling_Heavy_2009_2019.csv');
    traff_exchange    = sprintf('%s%s',ifold,'RTM_Trafikkutveksling_2020.xlsx');
    traff_exchange_sh = '2020_Utveksling';
    HBEFA_roads       = sprintf('%s%s',ifold,'Convert_RTM_roads_to_Match_HBEFA.xlsx');
    SSB_vehicle       = sprintf('%s%s',ifold,'SSB_Kjorelengder/Kjørelengder_2005_2019_prikket.csv');
    HBEFA_vehicles    = sprintf('%s%s',ifold,'HBEFA_4_National_Vehicle_Distribution.xlsx');
    SSB_Vehicle_dist  = sprintf('%s%s',ifold,'Model_Vehicles_Merge_SSB_and_HBEFA_Vehicles2.xlsx');
    Kartverket        = sprintf('%s%s',ifold,'Kartverket_fylker_og_kommuner_2020.xlsx');
    Bio_mix_file      = sprintf('%s%s',ifold,'Andel_biodrivstoff_2009_2019.xlsx');
    %Bio_mix_sheet     =
    %--------------------------------------------------------------------------
    % OPTIONS
    %--------------------------------------------------------------------------
    % % Set Model options:
    % List of compounds calculations should be done for:
    % %
    % comps = [{'FC_MJ'},{'BC'},{'PM'},{'HC'},{'CO'},{'NOx'},{'Be'},{'NMHC'},{'NO2'},{'PN'},{'N2O'},{'NH3'},{'CO2'},{'FC'},{'CH4'}];
    comps = [{'FC'},{'CH4'},{'CO2'},{'NO2'},{'BC'},{'FC_MJ'},{'PM'},{'HC'},{'CO'},{'NOx'},{'Be'},{'NMHC'}];
    % comps = [{'N2O'},{'NH3'},{'PN'}];
    % comps = [{'CO2'},{'FC'},{'N2O'},{'CH4'}];
    % comps = [{'CO2'},{'FC'}];
    % comps = [{'NOx'},{'NO2'},{'PM'}];
    % Vehicle_weight = 'Uniform';
    % Vehicle_weight = 'HBEFA';
    Vehicle_weight = 'NERVE';
    Subtract_bio_from_CO2  = 1;
    Vehicle_source         = {'SSB'};
    use_temporary_files    = 1;
    remove_NoTrafficRoads  = 1;
    do_preProcessing_HBEFA = 0;
    debug_mode             = 0; % debug_mode = Toggle what is written to screen (1= Alot 0 = Less)
    % Heavy_load_factor      = 0.6;
    % Buses_load_factor      = 0.6;
    % Set the filenames to be used:
    [input,tfiles,ofiles]   = Set_FileName_Struct();
    %--------------------------------------------------------------------------
    PreProcess_Emission_Factors_HBEFA();
    
    %--------------------------------------------------------------------------
    [Vehicle_dist] = PreProcess_Vehicle_Distribution();
    
    RLinks  = PreProcess_Roads_Input_Traffic();
    
    [EM_Links] = Emissions_Calculations();
    
    Save_shape(EM_Links,ofiles.RLShape,traffile)
    
end
% %--------------------------------------------------------------------------
% % Find the Distribution of Vehicles from SSB  :
% % Get or make the Roadlinks:
%
% %--------------------------------------------------------------------------
% Bio as before, December2019 roads
% ---- NORGE --- 2019 CO2 NAN
% ---- NORGE --- 2018 CO2
% ---- NORGE --- 2017 CO2
% ---- NORGE --- 2016 CO2
% ---- NORGE --- 2015 CO2
% ---- NORGE --- 2014 CO2
% ---- NORGE --- 2013 CO2
% ---- NORGE --- 2012 CO2
% ---- NORGE --- 2011 CO2
% ---- NORGE --- 2010 CO2
% ---- NORGE --- 2009 CO2


% Bio as before, December2019 roads
% ---- NORGE --- 2019 CO2 6366.4
% ---- NORGE --- 2018 CO2 6895.1
% ---- NORGE --- 2017 CO2 6654.8
% ---- NORGE --- 2016 CO2 7165.1
% ---- NORGE --- 2015 CO2 7585.0
% ---- NORGE --- 2014 CO2 7644.6
% ---- NORGE --- 2013 CO2 7641.1
% ---- NORGE --- 2012 CO2
% ---- NORGE --- 2011 CO2
% ---- NORGE --- 2010 CO2
% ---- NORGE --- 2009 CO2 7858.1

% NEW ROADS NEW Bio
% ---- NORGE --- 2019 CO2
% ---- NORGE --- 2018 CO2
% ---- NORGE --- 2017 CO2
% ---- NORGE --- 2016 CO2
% ---- NORGE --- 2015 CO2
% ---- NORGE --- 2014 CO2
% ---- NORGE --- 2013 CO2 7704.8
% ---- NORGE --- 2012 CO2 7667.6
% ---- NORGE --- 2011 CO2 7751.8
% ---- NORGE --- 2010 CO2 7826.5
% ---- NORGE --- 2009 CO2 7917.6

% NEW ROADS NEW Bio + Light
% ---- NORGE --- 2019 CO2 6574.0
% ---- NORGE --- 2018 CO2 7023.7
% ---- NORGE --- 2017 CO2 6796.0
% ---- NORGE --- 2016 CO2 7297.1
% ---- NORGE --- 2015 CO2 7676.1
% ---- NORGE --- 2014 CO2 7732.3
% ---- NORGE --- 2013 CO2 7714.5
% ---- NORGE --- 2012 CO2 7674.6
% ---- NORGE --- 2011 CO2 7755.9
% ---- NORGE --- 2010 CO2 7827.9
% ---- NORGE --- 2009 CO2 7917.6

% NEW ROADS NEW Bio + Light + Congestion
% ---- NORGE --- 2019 CO2 6718.6
% ---- NORGE --- 2018 CO2 7183.6
% ---- NORGE --- 2017 CO2 6964.9
% ---- NORGE --- 2016 CO2 7480.2
% ---- NORGE --- 2015 CO2 7871.4
% ---- NORGE --- 2014 CO2 7932.2
% ---- NORGE --- 2013 CO2 7918.3
% ---- NORGE --- 2012 CO2 7881.8
% ---- NORGE --- 2011 CO2 7971.3
% ---- NORGE --- 2010 CO2 8049.5
% ---- NORGE --- 2009 CO2 8142.6

%
%
%
% ---- NORGE --- 2019 NOx
% ---- Lette          17.9   (1000)Ton NOx ( 61%)
% ---- Tunge           9.9   (1000)Ton NOx ( 34%)
% ---- Busser          1.5   (1000)Ton NOx (  5%)
% ---- Totalt         29.3   (1000)Ton NOx        27.121 gG 2018ceip
%
%
% ---- NORGE --- 2019 NO2
% ---- Lette           6.0   (1000)Ton NO2 ( 82%)
% ---- Tunge           1.1   (1000)Ton NO2 ( 16%)
% ---- Busser          0.2   (1000)Ton NO2 (  3%)
% ---- Totalt          7.4   (1000)Ton NO2
%
%
% ---- NORGE --- 2019 PM
% ---- Lette           0.4   (1000)Ton PM ( 61%)
% ---- Tunge           0.2   (1000)Ton PM ( 34%)
% ---- Busser          0.0   (1000)Ton PM (  5%)
% ---- Totalt          0.6   (1000)Ton PM        1.114 gG 2018ceip
%
















