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
function [Sn] = Roads_Calc_HBEFA_Slope(RLinks)
    % Convert Precentage into HBEAF interpretable 0,2,4,6 pct slopes
    st = min(abs(round(50*extractfield(RLinks,'STIGNING_P'))*2),6);
    T = struct2table(RLinks);
    T.SLOPE      = st';
    Sn = table2struct(T);
    fprintf('Setfields\n  -- SLOPE ...on all links \n')
end