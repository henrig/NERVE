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
function OLinks =  Roads_Remove_NoTrafficRoads(RLinks)
fprintf('---------------------------------------------------------------\n')
fprintf('\nin Roads_Remove_NoTrafficRoads *\n...\n')
fprintf('---------------------------------------------------------------\n')

global remove_NoTrafficRoads
if remove_NoTrafficRoads == 0
    return
end
   
fieldN = fieldnames(RLinks);
t1 = contains(upper(fieldN),'ADT');
t2 = contains(upper(fieldN),'LETTE');

t = find(t1|t2);
if isempty(t)
    fprintf('### No ADT fields found ###\n')
    return
end

traff = zeros(1,length(RLinks));
for i = 1:length(t)
    if t(i)
        f = extractfield(RLinks,char(fieldN(t(i))));
        traff = traff+f;
        fprintf('%10s  %6i : %6i : %6i \n',char(fieldN(t(i))),length((find(f>0))),length((find(traff>0))),length(RLinks))
    end
end

OLinks = RLinks(traff>0);

end