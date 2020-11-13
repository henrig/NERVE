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
function  PreProcess_Emission_Factors_HBEFA()
%--------------------------------------------------------------------------
global use_temporary_files do_preProcessing_HBEFA comps tfold
% 22.10.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------

if use_temporary_files && ~do_preProcessing_HBEFA
    try
        for com = 1: length(comps)
            ifile = sprintf('%sEFA_Table_MODEL_%s.mat',tfold,char(comps(com)));
            if ~exist(ifile,'file')
                error(sprintf('### file not found\n%s\n will make new',ifile))
            end
        end
        fprintf('Temporary files found\n\t\tCONTINUE\n\n\n')
        return
    catch
    end
end

Emission_Factor_Process_HBEFA_Matrix_Raw()

Emission_Factor_Model_group_HBEFA()

end

