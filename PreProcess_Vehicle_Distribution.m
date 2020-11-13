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
function [Vehicle_dist] = PreProcess_Vehicle_Distribution()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global Vehicle_source 

fprintf('\t in PreProcess_Vehicle_Distribution\n')
Vehicle_dist=[];

if ismember(Vehicle_source,{'SSB'})
    fprintf('Using Municipal file from SSB for:\nVehcile_distribution\n\n')
    [Vehicle_dist] = Vehicle_Distribution_per_Municipality_SSB();
end

if ismember(Vehicle_source,{'HBEFA'})
    fprintf('Using National file from HBEFA for:\nVehcile_distribution\n\n')
    [Vehicle_dist] = Vehicle_Distribution_National_HBEFA();
    %%%% Vehicle,x2020,CLASS
    HBEFA_Make_Class_EF(Vehicle_dist)
end


if ismember(Vehicle_source,{'Scenario'})
    fprintf('Using Costum-made file for:\nVehcile_distribution\n\n')
    [Vehicle_dist] = Vehcile_Distribution_National_Scenario();
end
save(ofiles.MatlabOutput,'Vehicle_dist','-append')
end

