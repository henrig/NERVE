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
function [Sn] = Roads_Add_Urban(RLinks)
%-------------------------------------------------------------------------
% Miljodirektoratet Traffic emission model NERVE:
%
%     FUNCTION :: HEDGE_Add_ROAD_URBAN_RURAL ::
%
%
% -INPUT  : this function reads the road link shapes per municipality to the
%           temporary folder (tpath).
% -OUTPUT : this function adds  properties to the road link shapes per
%           municipality in the temporary folder (tpath).
%
% 20.09.2018 -Henrik Grythe
% Kjeller NILU
%-------------------------------------------------------------------------
global tettfile debug_mode

fprintf('* in Add_ROAD_URBAN_RURAL   *\n')
fprintf('Adding URBAN field based on geo file: \n%s\n',tettfile)
Tshape  = shaperead(tettfile);
prj     = read_projection(tettfile);

% Adds a Roads Midpoint based on the starting point and ending point of
% Road
if ~isfield(RLinks,'MIDPOINT_X') || ~isfield(RLinks,'MIDPOINT_Y')
    for i=1:length(RLinks)
        x=RLinks(i).X;
        y=RLinks(i).Y;
        MIDPOINT_X(i)=(x(1)+x(end-1))/2;
        MIDPOINT_Y(i)=(y(1)+y(end-1))/2;
    end
else
     MIDPOINT_X = extractfield(RLinks,'MIDPOINT_X');
     MIDPOINT_Y = extractfield(RLinks,'MIDPOINT_Y');
end
% sped up calculations, check that they give same result
% Only checking Raods close by (closest 20 km of tettsted geographic center)
max_dist = 27;
I  = [];
ISURBAN = zeros(length(RLinks),1); 
for i=1:length(Tshape)
     mx=nanmean(Tshape(i).X);    
     my=nanmean(Tshape(i).Y);    
     dst =sqrt((MIDPOINT_X-mx).^2 + (MIDPOINT_Y-my).^2)*1e-3;    
     idx = find(dst<max_dist);
     fprintf('Checking %6i roads in tettsted %04i_%-18s',length(idx),i,char(Tshape(i).Tettstedsn))
     IN = inpolygon(MIDPOINT_X(idx), MIDPOINT_Y(idx),Tshape(i).X,Tshape(i).Y);
     I  = cat(2,I,idx(IN));
     fprintf('... Found %6i roads \n',length(idx(IN)))
end
ISURBAN(I)  = 1; 

fprintf('Of %6i Roads \n We forund\n %6i URBAN\n %6i RURAL roads\n',length(ISURBAN),sum(ISURBAN),sum(ISURBAN==0))
fprintf('Adding field to Table ...\n')
T        = struct2table(RLinks);
T.URBAN  = ISURBAN;
fprintf('Structuring Table ...\n')
Sn       = table2struct(T);

fprintf('Setfields\n  -- URBAN:   [YES=1,NO=0] ...on all links \n')
end
