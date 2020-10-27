function [Sn] = Emissions_Calculations(Calc_Links,Vehicle_dist)
%--------------------------------------------------------------------------
%
% 22.10.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------

% Calculate_Emissions calculates emissions 
global Vehicle_source

if ismember(Vehicle_source,{'HBEFA'})
    [Sn] = Emissions_Calculations_HBEFA(Calc_Links);
    return
end

if ismember(Vehicle_source,{'SSB'})
    [Sn] = Emissions_Calculations_SSB(Calc_Links,Vehicle_dist);
    return
end


end
