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
function [Sn] = Roads_Add_Municipality(RLinks)
% This function looks for and checks the national traffic file. It makes
% sure that it has the required fields and adds (some) fields if they are
% deducable from the input-traffic. As there is a lot of uneccesary data
% and variable number and properties of the roads, properties and values
% with each iteration of the traffic input file, this is spesifically
% designed to treat the ones we have recieved. As the shapefiles recieved
% are unique in design, and different every time this function mainly deals
% with known problems. Field names for example are treated by name and
% may not be same in future extractions, and a change will have to be made
% in the script.
%
%-------------------------------------------------------------------------
% 20.09.2018 -Henrik Grythe
% Kjeller NILU
%-------------------------------------------------------------------------
global Komm_shape
fprintf('---------------------------------------------------------------\n')
fprintf('in Roads_Add_Municipality   *\n')
fprintf('---------------------------------------------------------------\n')

fprintf('Adding kommune field based on geo file: \n%s\n',Komm_shape)
Ks      = shaperead(Komm_shape);

% Adds a Roads Midpoint based on the starting point and ending point of
% Road
% if ~isfield(RLinks,'MIDPOINT_X') || ~isfield(RLinks,'MIDPOINT_Y')
%     for i = 1:length(RLinks)
%         x            = RLinks(i).X;
%         y            = RLinks(i).Y;
%         MIDPOINT_X(i)=(x(1)+x(end-1))/2;
%         MIDPOINT_Y(i)=(y(1)+y(end-1))/2;
%     end
% else
%     MIDPOINT_X = extractfield(RLinks,'MIDPOINT_X');
%     MIDPOINT_Y = extractfield(RLinks,'MIDPOINT_Y');
% end

for i = 1:length(RLinks)
    x            = RLinks(i).X;
    y            = RLinks(i).Y;
    ix           = find(~isnan(x));
    iy           = find(~isnan(y));
    START_X(i)   = x(ix(1));
    START_Y(i)   = y(iy(1));
    END_X(i)     = x(ix(end));
    END_Y(i)     = y(iy(end));
    MIDPOINT_X(i) = mean(x(ix));
    MIDPOINT_Y(i) = mean(y(iy));        
end

KOMMS = nan(size(START_X));
KOMME = nan(size(START_X));
KOMM  = nan(size(START_X));
for i = 1:length(Ks)
    % FIND STARTING MUNICIPALITY
    inS = inpolygon(START_X,START_Y,Ks(i).X,Ks(i).Y);
    % FIND ENDING MUNICIPALITY
    inE = inpolygon(END_X,END_Y,Ks(i).X,Ks(i).Y);
    % FIND MIDPOINT MUNICIPALITY
    inM = inpolygon(MIDPOINT_X,MIDPOINT_Y,Ks(i).X,Ks(i).Y);
    KOMMS(inS) = extractfield(Ks(i),'KOMMUNENUM');
    KOMME(inE) = extractfield(Ks(i),'KOMMUNENUM');
    KOMM(inM)  = extractfield(Ks(i),'KOMMUNENUM');
    hasS(inS)  = 1;
    hasE(inE)  = 1;
    hasM(inM)  = 1;
    fprintf('%03i Kommune %04i_%-24s has %5i roadS \t',i,extractfield(Ks(i),'KOMMUNENUM'),char(extractfield(Ks(i),'NAVN')),sum(inS))
    fprintf('(%i/%i = %4.1f%%) \n',sum(hasS),length(RLinks), 100*sum(hasS)/length(RLinks))
    clear inS inE inM
end
fprintf('Have: %7i of %7i Roadlinks "KOMM START", missing: %i \n',sum(hasS),length(RLinks),length(RLinks)-sum(hasS))
fprintf('Have: %7i of %7i Roadlinks "KOMM END  ", missing: %i \n',sum(hasE),length(RLinks),length(RLinks)-sum(hasE))
fprintf('Have: %7i of %7i Roadlinks "KOMM MIDP.", missing: %i \n',sum(hasM),length(RLinks),length(RLinks)-sum(hasM))

% teller = 0;
% try midpoint start point and end point of road, which commune it
% will be decided.
missing = find(hasM==0);
if ~isempty(missing)
    % If midpoint doesnt work, do it on the starting point
    for i=1:length(missing)
        x=extractfield(RLinks(missing(i)),'X');
        y=extractfield(RLinks(missing(i)),'Y');
        for j=1:length(Ks)
            in = inpolygon(x(1),y(1),Ks(j).X,Ks(j).Y);
            if ~isempty(in)
                KOMM(missing(i)) = extractfield(Ks(j),'KOMMUNENUM');
                hasM(missing(i))  = 1;
            end
        end
    end
end
fprintf('Have: %7i of %7i Roadlinks "KOMM", missing: %i \n',sum(hasM),length(RLinks),length(RLinks)-sum(hasM))

missing = find(hasM==0);
if ~isempty(missing)
    % If startpoint doesnt work, do it on the ending point
    for i=1:length(missing)
        x=extractfield(RLinks(missing(i)),'X');
        y=extractfield(RLinks(missing(i)),'Y');
        for j=1:length(Ks)
            in = inpolygon(x(end-1),y(end-1),Ks(j).X,Ks(j).Y);
            if ~isempty(in)
                KOMM(missing(i)) = extractfield(Ks(j),'KOMMUNENUM');
                hasM(missing(i))  = 1;
            end
        end
    end
end

fprintf('Have: %7i of %7i Roadlinks "KOMM", missing: %i \n',sum(hasM),length(RLinks),length(RLinks)-sum(hasM))

% Convert briefly to a table variable to make some calculations
% that take more time in structures.
T=struct2table(RLinks);
T.KOMM   = KOMM';
T.KOMMS  = KOMMS';
T.KOMME  = KOMME';
T.MIDPOINT_X = MIDPOINT_X';
T.MIDPOINT_Y = MIDPOINT_Y';
fprintf('Setfields\n  -- KOMM  Kommune (Mid)\n  -- KOMMS Start Kommune \n')
fprintf('  -- KOMME Endekommune \n  -- MIDPOINT_X\n  -- MIDPOINT_Y ...on all links \n')
Sn  = table2struct(T);
end
