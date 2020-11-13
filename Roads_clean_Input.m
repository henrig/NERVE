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
function [Sn] = Roads_clean_Input(RLinks)
global Inputfields  Listfields
fprintf('---------------------------------------------------------------\n')
fprintf('in Roads_clean_Input  *\n')
fprintf('---------------------------------------------------------------\n')

% Tabelize Roads
Tlink   = struct2table(RLinks);
fields  = unique([Inputfields,Listfields]);
f = zeros(length(fields),1)
for i=1:length(fields)
    % fprintf('%s\n',char(fields(i)))
    try
        f(i)  = find(ismember(Tlink.Properties.VariableNames,fields(i)));
    catch
        fprintf('### !! MISSING FIELD !! : %s\n',char(fields(i)))
    end
    
end
f = f(find(f));

Tlink.Properties.VariableNames(f);

Sn = table2struct(Tlink(:,f));
fprintf('--->\n')
end

