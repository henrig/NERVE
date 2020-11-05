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
function [Sn] = Roads_Scale_Traffic_to_Year(RLinks)
% Scale_Road_Traffic_Year function takes a traffic file and scales it
% from the year (Ryear) to a new year (Tyear) based on the csv file. It
% does scaling at the municipal basis and treats each municipality
% independently.

% It also makes a lighter version of the Roads that only keep the
% information that is required to do emission calculations. 

global Ryear Tyear Light_traff_years Heavy_traff_years Buses_traff_years
global Listfields


fprintf('in Scale_Road_Traffic_Year  \n')


fprintf('To scale Light Traffic from %i to %i \n',Ryear,Tyear)

% Read the Light csv scaling file into table
fprintf('Reading Annual Scaling files (traff_years):\n %s \n',Light_traff_years)
T       = readtable(Light_traff_years,'ReadVariableNames',1);
fldList = T.Properties.VariableNames;
f       = find(ismember(fldList,sprintf('x%i',Ryear)));
fy      = find(ismember(fldList,sprintf('x%i',Tyear)));

% Extract the scale 
Scale         = table;
Scale.KOMM    = T.Kommune;
Scale.ToYearL = table2array(T(:,fy))./table2array(T(:,f));
k = unique(extractfield(RLinks,'KOMM'));

% Read the Heavy csv scaling file into table
fprintf(' %s \n',Heavy_traff_years)
TH       = readtable(Heavy_traff_years,'ReadVariableNames',1);
fldListH = T.Properties.VariableNames;
fH       = find(ismember(fldList,sprintf('x%i',Ryear)));
fyH      = find(ismember(fldList,sprintf('x%i',Tyear)));
Scale.ToYearH = table2array(TH(:,fyH))./table2array(TH(:,fH));

% Read the Heavy csv scaling file into table
fprintf(' %s \n',Buses_traff_years)
TB       = readtable(Buses_traff_years,'ReadVariableNames',1);
fldListB = T.Properties.VariableNames;
fB       = find(ismember(fldList,sprintf('x%i',Ryear)));
fyB      = find(ismember(fldList,sprintf('x%i',Tyear)));
Scale.ToYearB = table2array(TB(:,fyB))./table2array(TB(:,fB));

fprintf('Flat Average Scale\nFrom %i to %i:\n',Ryear,Tyear)
fprintf('Light=%4.1f%%\n',100*(mean(Scale.ToYearL)-1))
fprintf('Heavy=%4.1f%%\n',100*(mean(Scale.ToYearH)-1))
fprintf('Buses=%4.1f%%\n',100*(mean(Scale.ToYearB)-1))

% Tabelize Roads
Tlink   = struct2table(RLinks);

found_err = 0;
L_adt = zeros(length(RLinks),1);
H_adt = zeros(length(RLinks),1);
B_adt = zeros(length(RLinks),1);
for i = 1:length(k)
    idk = find(Scale.KOMM==k(i));
    if ~isempty(idk)
        idr = find(Tlink.KOMM == k(i));
        if ~isempty(idr)
            L_adt(idr) = Tlink.LETTE_BILE(idr)*Scale.ToYearL(idk);
            H_adt(idr) = Tlink.GODS_ADT(idr)  *Scale.ToYearH(idk);
            B_adt(idr) = Tlink.KOLL_ADT(idr)  *Scale.ToYearB(idk);
        else
            fprintf('### Missinig roads for Municipality %i ###\n', k(i))
            found_err = 1;
        end
    else
        fprintf('### Missinig scale for Municipality %i ###\n', k(i))
        found_err = 1;
    end
end
if found_err
    warning('### could not complete calculations for all roads ###')
end 
fprintf('Scaled all Municipalities \n')
Tlink.L_ADT = L_adt;
Tlink.H_ADT = H_adt;
Tlink.B_ADT = B_adt;



fprintf('Structuring Roads \n')
light_adt = find(ismember(Tlink.Properties.VariableNames,'L_ADT'));
heavy_adt = find(ismember(Tlink.Properties.VariableNames,'H_ADT'));
bus_adt   = find(ismember(Tlink.Properties.VariableNames,'B_ADT'));

Tlink.Properties.VariableNames(light_adt) = {sprintf('L_ADT%04i',Tyear)};
Tlink.Properties.VariableNames(heavy_adt) = {sprintf('H_ADT%04i',Tyear)};
Tlink.Properties.VariableNames(bus_adt)   = {sprintf('B_ADT%04i',Tyear)};

Sn = table2struct(Tlink);

% clear f
% for i=1:length(Listfields)
%    f(i)  = find(ismember(Tlink.Properties.VariableNames,Listfields(i)));
% end
% 
% light_adt = find(ismember(Tlink.Properties.VariableNames,'L_ADT'));
% heavy_adt = find(ismember(Tlink.Properties.VariableNames,'H_ADT'));
% bus_adt   = find(ismember(Tlink.Properties.VariableNames,'B_ADT'));
% 
% Tlink.Properties.VariableNames(light_adt) = {sprintf('L_ADT%04i',Tyear)};
% Tlink.Properties.VariableNames(heavy_adt) = {sprintf('H_ADT%04i',Tyear)};
% Tlink.Properties.VariableNames(bus_adt)   = {sprintf('B_ADT%04i',Tyear)};
% 
% Sn = table2struct(Tlink(:,f));
end