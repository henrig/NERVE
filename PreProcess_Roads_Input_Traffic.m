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
function [RLinks] = PreProcess_Roads_Input_Traffic()
%--------------------------------------------------------------------------
% Module to read in an original Road Network and traffic file.
% Processes data to calculate the required fields and splits it to
% municipal level data.
% Function is not generic and require input of specified format (.shp) and
% also require existing fieldnames specified.
%
% -INPUT  : this function reads road and traffic network from file (ESRI shape)
% -OUTPUT : this function adds the road link shapes per municipality to the
%           temporary folder (tpath).
%
% Road sources are generally considered to be ESRI shapefiles by the model.
% They have an assumed set of fields ('Geometry','BoundingBox','X','Y') and
% traffic volumes given as ADT (Annual daily traffic) for Light, Heavy and
% Buses. In addition the Model wants some advanced road properties such as the
% slope, type of road, speed, congestion etc that it uses for Advanced EF
% calculations on each road link.
%
% Functions also read ancillary files to define fields.
%   read_projection()
%   Roads_Add_Urban()
%   Roads_Add_HBEFA_Parameters()
%   Roads_Add_Advanced_Properties() 
%   Roads_Add_Width()
%   Roads_Scale_Traffic_to_Year()
%   Roads_Congestion_Parameters()
%   Roads_Fix_Roadfields
%
% 20.09.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global use_temporary_files tfiles Tyear
global traffile Listfields Inputfields
fprintf('---------------------------------------------------------------\n')
fprintf('in PreProcess_Input_Traffic *\n')
fprintf('---------------------------------------------------------------\n')

% List of input fields used
Inputfields = [{'Geometry'},{'X'},{'Y'},{'KOLL_ADT'},{'LETTE_BILE' },...
    {'GODS_ADT'},{'SUM_ADT'},{'FM_TIME'},{'KO_MORGEN'},{'KO_ETTERM'},...
    {'FM_SPEED'},{'SHAPE_LENG'},{'VK'},{'STIGNING_P'},{'SPEED'},{'HP_ID'},{'LANES'},{'DEKKEBREDD'}];

% List of output fields a finished file has:
Listfields = [{'Geometry'},{'X'},{'Y'},{'BoundingBox'},{'DISTANCE'},{'KOMMS'},{'KOMME'},{'KOMM'},{'SPEED'},{'SLOPE'},...
    {'URBAN'},{'RUSH_DELAY'},{'HBEFA_EQIV'},{'WIDTH'},{'CAPACITY'},{'N_LANES'},{'IDO'}];

if use_temporary_files
    try
        % First do tests if there is a temporary file saved:
        fprintf('Trying to read temp file\n%s ...\n',tfiles.RL)
        RLinks = shaperead(tfiles.RL);
        fprintf('Found Pre-Processed file:\n %s \n',tfiles.RL)
        fields = fieldnames(RLinks);
        fprintf('Road Links found %i\n',length(RLinks))
        fprintf('Fields:\n')
        for i=1:length(fields)
            fprintf('%s\n',char(fields(i)))
        end
        return
    catch
        fprintf('### Temporary file Scaled to Year %i NOT found\n',Tyear)
        try
            RLinks = shaperead(tfiles.CleanRoads);
            fprintf('Temporary CLEANED file found \n')
            RLinks = Roads_Scale_Traffic_to_Year(RLinks);
            RLinks = Roads_Clean_Annual(RLinks);
            Save_shape(RLinks,tfiles.RL)

            return
        end
        prj    = read_projection(traffile);
        fprintf('Reading large file: \n\t %s  ...',traffile)
        RLinks = shaperead(traffile);
        fprintf('done\n')
        fields = fieldnames(RLinks);
        fprintf('Fields:\n')
        for i=1:length(fields)
            fprintf('%s,',char(fields(i)))
            if rem(i,10)==0;fprintf('\n');end
        end
        fprintf('\n')
    end
else
    prj    = read_projection(traffile);
    fprintf('Reading large file: \n\t %s  ...',traffile)
    RLinks = shaperead(traffile);
    fprintf('done\n')
    fields = fieldnames(RLinks);
    fprintf('Fields:\n')
    for i=1:length(fields)
        fprintf('%s,',char(fields(i)))
        if rem(i,10)==0;fprintf('\n');end
    end
    fprintf('\n')
end
%--------------------------------------------------------------------------
rRLinks = Roads_clean_Input(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Remove_NoTrafficRoads(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Congestion_Parameters(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Calc_DISTANCE(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Calc_HBEFA_Slope(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Add_Municipality(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Find_StartMuncipalityFraction(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Fix_Roadfields(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Add_Urban(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Add_HBEFA_Parameters(RLinks);
%--------------------------------------------------------------------------
RLinks = Roads_Congestion_Parameters(RLinks);
%--------------------------------------------------------------------------
if use_temporary_files
    Save_shape(RLinks,tfiles.CleanRoads)
end

CLinks = Roads_Scale_Traffic_to_Year(RLinks);
CLinks = Roads_Clean_Annual(CLinks);
%--------------------------------------------------------------------------
if use_temporary_files
    Save_shape(CLinks,tfiles.RL)
end
fprintf('\n\n\n\n ROADS PROCESSED \n\n\n\n\n\n')
end
