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
fprintf('---------------------------------------------------------------\n')
fprintf('in Roads_Congestion_Parameters *\n')
fprintf('---------------------------------------------------------------\n')

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

if isfield(RLinks,'KAPTID_16_')
    a(:,1) = extractfield(RLinks,'KAPTID_16_');
    kvl = max([kvl,a],[],2);
end

if isfield(RLinks,'KAPTID_17_')
    a(:,1) = extractfield(RLinks,'KAPTID_17_');
    kvl = max([kvl,a],[],2);
end

if sum(mrg)>0 && sum(kvl)>0
    T = struct2table(RLinks);
    T.KO_MORGEN = mrg;
    T.KO_ETTERM = kvl;
    fprintf('Structuring Roads...\n')
    RLinks = table2struct(T);
    fprintf('Setfields\n')
    fprintf('  -- KO_ETTERM\n')
    fprintf('  -- KO_MORGEN ...on all links\n')
else

di = extractfield(RLinks,'DISTANCE');
kt = extractfield(RLinks,'FM_TIME');
km = extractfield(RLinks,'KO_MORGEN');
ke = extractfield(RLinks,'KO_ETTERM');
ks = extractfield(RLinks,'FM_SPEED');

% The time traffic that is delayed is calculated by the mean of the
% morning and evening congestion
RUSH_DELAY                    = round(100*(kt+mean([ke;km]))./kt-100);
RUSH_DELAY(isnan(RUSH_DELAY)) = 0;
RUSH_SPEED                    = round(di./((kt.*(1+RUSH_DELAY/100)/60)));
RUSH_SPEED(isnan(RUSH_SPEED)) = 0;
RUSH_SPEED(RUSH_SPEED>kt)     = kt(RUSH_SPEED>kt);

T = struct2table(RLinks);
T.RUSH_DELAY = RUSH_DELAY';
T.RUSH_SPEED = RUSH_SPEED';
Sn = table2struct(T);
end


