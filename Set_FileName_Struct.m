function [input, tfiles, ofiles] = Set_FileName_Struct()
global do_preProcessing_HBEFA debug_mode use_temporary_files ofiles tfiles
global Vehicle_source remove_NoTrafficRoads Vehicle_weight Bio_mix_file
% folders
global ifold tfold ofold HBEFA_path HBEFA_roads Kartverket
% files
global traffile tettfile Komm_shape Light_traff_years Heavy_traff_years
global Buses_traff_years HBEFA_vehicles SSB_vehicle tfiles traff_exchange
global traff_exchange_sh SSB_Vehicle_dist Subtract_bio_from_CO2
% variables
global Ryear Tyear Myear comps Vehicle_dist RLinks
fprintf('---------------------------------------------------------------\n')
fprintf('in Set_FileName_Struct\n')
fprintf('---------------------------------------------------------------\n')

input.Tyear      = Tyear;
input.Myear      = Myear;
input.Components = comps; 
input.VehcileSource = Vehicle_source;

input.options.do_preProcessing_HBEFA = do_preProcessing_HBEFA;
input.options.debug_mode             = debug_mode;
input.options.use_temporary_files    = use_temporary_files; 
input.options.Vehicle_source         = Vehicle_source; 
input.options.remove_NoTrafficRoads  = remove_NoTrafficRoads; 
input.options.Vehicle_weight         = Vehicle_weight; 
input.options.Subtract_bio_from_CO2  = Subtract_bio_from_CO2;

input.ifold                   = ifold;
input.tfold                   = tfold;
input.tfolde                  = tfold;
input.HBEFA_path              = HBEFA_path;

input.files.HBEFA_roads       = HBEFA_roads;
input.files.Kartverket        = Kartverket;
input.files.Trafficfile       = traffile;
input.files.Tettstedfile      = tettfile;
input.files.Kommune_shape     = Komm_shape;
input.files.Light_traff_scale = Light_traff_years;
input.files.Heavy_traff_scale = Heavy_traff_years;
input.files.Buses_traff_scale = Buses_traff_years;
input.files.HBEFA_vehicles    = HBEFA_vehicles;
input.files.SSB_vehicle       = SSB_vehicle; 
input.files.traff_exchange    = traff_exchange; 
input.files.traff_exchange_sh = traff_exchange_sh; 
input.files.SSB_Vehicle_dist  = SSB_Vehicle_dist; 
input.files.Bio_mix_file      = Bio_mix_file;

tfiles.EF           = sprintf('%sHEDGE_EF_%4i',tfold,Tyear);
tfiles.CarParks     = sprintf('%sSSB_Vehicle_data',tfold);
tfiles.CarPark      = sprintf('%sSSB_CarPark_%4i',tfold,Tyear);
tfiles.Exchange     = sprintf('%sMunicipal_Traffic_Exchange_%04i.mat',tfold,Tyear);
tfiles.DD_Municipal = sprintf('%sRoadTypeDistanceMunicipal%4i',ofold,Tyear);

if input.options.remove_NoTrafficRoads
    tfiles.CleanRoads = sprintf('%s_Cleaned0Removed',input.files.Trafficfile);
else
    tfiles.CleanRoads = sprintf('%s_Cleaned',input.files.Trafficfile);
end    

a = strfind(traffile,'/');
if ~isempty(a)
    tfiles.RL      = sprintf('%sHEDGE_RL_%4i_%s_temp',tfold,Tyear,traffile(a(end)+1:end));
else
    tfiles.RL      = sprintf('%sHEDGE_RL_%4i_%s_temp',tfold,Tyear,traffile);
end

ofiles.RLShape              = sprintf('%sTraffic_Emissions_%4i',ofold,Tyear);
ofiles.SSB_Stat             = sprintf('%sNational_SSB_Vehicle_Number_Distribution.xlsx',ofold); 
ofiles.MatlabOutput         = sprintf('%sOutput_%4i_Simulation_%s.mat',ofold,Tyear,datestr(now,'yyyymmmdd_HH'));
ofiles.OnRoadEF_RoadClasses = sprintf('%s%s',ofold,'OnRoadEF_RoadClasses.xlsx');

fprintf('initiated new output file: \n %s\n',ofiles.MatlabOutput)
save(ofiles.MatlabOutput,'input','tfiles','ofiles')


end
