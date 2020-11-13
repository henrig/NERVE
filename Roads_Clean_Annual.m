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
function [Sn] = Roads_Clean_Annual(RLinks)
global Listfields Tyear
fprintf('---------------------------------------------------------------\n')
fprintf('in Roads_clean_Annual  *\n')
fprintf('---------------------------------------------------------------\n')

% Tabelize Roads
Tlink   = struct2table(RLinks);
fprintf('--- Structuring Roads --- \n')

ladt = find(ismember(Tlink.Properties.VariableNames,sprintf('L_ADT%i',Tyear)));
hadt = find(ismember(Tlink.Properties.VariableNames,sprintf('H_ADT%i',Tyear)));
badt = find(ismember(Tlink.Properties.VariableNames,sprintf('B_ADT%i',Tyear)));

Outfields = [Listfields, Tlink.Properties.VariableNames(ladt), ...
    Tlink.Properties.VariableNames(hadt), Tlink.Properties.VariableNames(badt)];

clear f
for i=1:length(Outfields)
    try
        f(i)  = find(ismember(Tlink.Properties.VariableNames,Outfields(i)));
    catch
        fprintf('### !! MISSING OUTPUT FIELD !! : %s\n',char(Outfields(i)))
    end
    
end
f = f(find(f));
Sn = table2struct(Tlink(:,f));

fprintf('--->\n')
end

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