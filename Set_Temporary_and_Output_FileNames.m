function [tfiles ofiles] = Set_Temporary_and_Output_FileNames()

global Tyear ofold tfold traffile


tfiles.EF      = sprintf('%sHEDGE_EF_%4i',tfold,Tyear);

tfiles.CarPark = sprintf('%sSSB_CarPark_%4i',ofold,Tyear);


a = strfind(traffile,'/');
if ~isempty(a)
    tfiles.RL      = sprintf('%sHEDGE_RL_%4i_%s_temp',tfold,Tyear,traffile(a(end)+1:end));
else
    tfiles.RL      = sprintf('%sHEDGE_RL_%4i_%s_temp',tfold,Tyear,traffile);
end

ofiles.RLShape  = sprintf('%sTraffic_Emissions_%4i',ofold,Tyear);
ofiles.SSB_Stat = sprintf('%sNational_SSB_Vehicle_Number_Distribution.xlsx',ofold); 

ofiles.MatlabOutput = sprintf('%sOutput.mat',ofold);

end