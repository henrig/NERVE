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
function Vehicle_Distribution_Statistical_output(Vehicle_data,T)
%--------------------------------------------------------------------------
% National Statistical output.

% National Number of vehicles
global ofiles
National = squeeze(sum(Vehicle_data.YNV,2));
k=1;
National2D = zeros(size(National,2)*size(National,3),size(National,1));
Ts = table;
for k = 1:size(National,1)
    t = 1;
    Ttemp          = table;
    for j = 1:size(National,3)
        for i = 1:size(National,2)
            if k==1
                Ttemp          = table;
                Ttemp.Num      = t;
                Ttemp.Euro     = Vehicle_data.D4_euro(j);
                Ttemp.Kategori = Vehicle_data.D3_nybiltype(i);
                Ttemp.Var1     = National(k,i,j);
                Ts             = [Ts;Ttemp];
                t              = t+1;
            else
                Ts.Var1(t)  = National(k,i,j);
                t              = t+1;
            end
        end
    end
    vn = find(ismember(Ts.Properties.VariableNames,'Var1'));
    Ts.Properties.VariableNames(vn)= {sprintf('x%i',Vehicle_data.D1_yrs(k))};
end
writetable(T,ofiles.SSB_Stat,'Sheet','NationalNumberVehicles')
fprintf('Wrote National Vehicle stats file :\n%s\n',ofiles.SSB_Stat)
% National Driving Distance
National = squeeze(sum(Vehicle_data.YTD,2));
k=1;
National2D = zeros(size(National,2)*size(National,3),size(National,1));
Ts2 = table;
for k = 1:size(National,1)
    t = 1;
    Ttemp          = table;
    for j = 1:size(National,3)
        for i = 1:size(National,2)
            if k==1
                Ttemp          = table;
                Ttemp.Num      = t;
                Ttemp.Euro     = Vehicle_data.D4_euro(j);
                Ttemp.Kategori = Vehicle_data.D3_nybiltype(i);
                Ttemp.Var1     = National(k,i,j);
                Ts2 = [Ts2;Ttemp];
                t=t+1;
            else
                Ts2.Var1(t)  = National(k,i,j);
                t              = t+1;
            end
        end
    end
    vn = find(ismember(Ts2.Properties.VariableNames,'Var1'));
    Ts2.Properties.VariableNames(vn)= {sprintf('x%i',Vehicle_data.D1_yrs(k))};
end
writetable(Ts2,ofiles.SSB_Stat,'Sheet','NationalTotalDistance')
end
