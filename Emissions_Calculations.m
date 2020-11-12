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
function [Sn] = Emissions_Calculations()
%--------------------------------------------------------------------------
%
% 22.10.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global RLinks 
% Calculate_Emissions calculates emissions 
global Vehicle_source

if ismember(Vehicle_source,{'HBEFA'})
    [Sn] = Emissions_Calculations_HBEFA(RLinks);
    return
end

if ismember(Vehicle_source,{'SSB'})
    Emission_Factors_OnRoadAllCond()
    Emission_Factors_Road_DrivingDistance_IN_Municipalities()
    [Sn] = Emissions_Calculations_SSB();
    return
end


end
