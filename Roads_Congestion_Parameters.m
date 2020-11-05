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
function Sn = Roads_Congestion_Parameters(RLinks)

fprintf('\n in Roads_Congestion_Parameters \n...\n')

mrg = zeros(size(RLinks));
% Morgen
if isfield(RLinks,'KAPTID_06_')
    a(:,1) = extractfield(RLinks,'KAPTID_06_');
    mrg = max([mrg,a],[],2);
end
if isfield(RLinks,'KAPTID_07_')
    a(:,1) = extractfield(RLinks,'KAPTID_07_');
    mrg = max([mrg,a],[],2);
end

if isfield(RLinks,'KAPTID_08_')
    a(:,1) = extractfield(RLinks,'KAPTID_08_');
    mrg = max([mrg,a],[],2);
end

% Ettermiddag
kvl = zeros(size(RLinks));
if isfield(RLinks,'KAPTID_15_')
    a(:,1) = extractfield(RLinks,'KAPTID_15_');
    kvl = max([kvl,a],[],2);
end

if isfield(RLinks,'KAPTID_06_')
    a(:,1) = extractfield(RLinks,'KAPTID_16_');
    kvl = max([kvl,a],[],2);
end

if isfield(RLinks,'KAPTID_06_')
    a(:,1) = extractfield(RLinks,'KAPTID_17_');
    kvl = max([kvl,a],[],2);
end
T = struct2table(RLinks);
T.KO_MORGEN = mrg;
T.KO_ETTERM = kvl;
fprintf('Structuring Roads...\n')
Sn = table2struct(T);
fprintf('Setfields\n')
fprintf('  -- KO_ETTERM\n')
fprintf('  -- KO_MORGEN ...on all links\n') 
end


