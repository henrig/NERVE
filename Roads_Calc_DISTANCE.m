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
function [Sn] = Roads_Calc_DISTANCE(RLinks)
fprintf('in Roads_Calc_DISTANCE.\n')

if isfield(RLinks,'DISTANCE')
    old_dst = extractfield(RLinks,'DISTANCE');
    mean_length =  nanmean(old_dst);
elseif isfield(RLinks,'LEN')
    old_dst = extractfield(RLinks,'LEN');
    mean_length =  nanmean(old_dst);
else
    fprintf('### NO Distance field found \n')
    old_dst = nan(size(RLinks));
    mean_length = NaN;
end

if mean_length < 10
    unit = 1000;
    fprintf('Old Unit Km."mean_length < 10" \n')
else
    unit = 1;
    fprintf('Old Unit meter "mean_length > 10" \n')
    old_dst = old_dst/1000;
end

dst = zeros(size(RLinks));
k=0;
for i=1:length(RLinks)
    xt = extractfield(RLinks(i),'X');
    yt = extractfield(RLinks(i),'Y');
    for j = 2:length(xt)
        dst_x = abs(xt(j-1)-xt(j));
        dst_y = abs(yt(j-1)-yt(j));
        if ~isnan(dst_x) && ~isnan(dst_y)
            dst(i)= dst(i)+(sqrt(dst_x^2 + dst_y^2))/1000;            
        end
    end
    if ((old_dst(i)-dst(i))/dst(i)> 0.5 || (old_dst(i)-dst(i))/dst(i)< -0.5) && ~isnan(mean_length)
       k = k+1;
       fprintf('i = %6i Old distance / New distance / Diff %fkm/%fkm/%f%% --- Segments %i \n',i,old_dst(i),dst(i),100*(old_dst(i)-dst(i))/dst(i),length(xt)-2)
    end
end


fprintf('MEAN  distance old / new \n\t %f %f change = %5.1f%%\n',mean(old_dst),mean(dst),100*(mean(old_dst)-mean(dst))/mean(old_dst))
fprintf('TOTAL distance old / new \n\t %f %f change = %5.1f%%\n',sum(old_dst),sum(dst),100*(sum(old_dst)-sum(dst))/sum(old_dst))
fprintf('MIN   distance old / new \n\t %f %f change = %5.1f%%\n',min(old_dst),min(dst),100*(min(old_dst)-min(dst))/min(old_dst))
fprintf('MAX   distance old / new \n\t %f %f change = %5.1f%%\n',max(old_dst),max(dst),100*(max(old_dst)-max(dst))/max(old_dst))

TLinks = struct2table(RLinks);
TLinks.DISTANCE = dst;
Sn = table2struct(TLinks);

end

