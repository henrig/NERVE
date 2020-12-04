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
function [fra] = Roads_inDomain_Roads(S,E)
global debug_mode
% S = t_RLinks(i)
% E = Ks(pla)
fra = NaN;
for s=1:length(S)
    % This does not necessarily find all the intersection
    in = inpolygon(S(s).X,S(s).Y,E.X,E.Y);
    nn = sum(isnan(S(s).X));
    DST(s) = S(s).DISTANCE;
    if sum(in)   == 0                   % No vertices inside the domain: Road is outside and will be cut
        fra = 0;
    elseif sum(in) == length(S(s).X)-nn % All vertices inside the domain Road is inside and will be used
        fra = 1;
    else                                % There are intersections, partly in domain. Calculate the part inside the domain.
        xt = S(s).X(in);
        yt = S(s).Y(in);
        
        % This test sometimes fail to find all the intersections
        if length(S(s).X) < 500
            [xi,yi] = polyxpoly(S(s).X,S(s).Y,E.X,E.Y);
        else
            if debug_mode
                warning(sprintf('Straightened Road with %i subnodes',length(S(s).X)))
            end
            [xi,yi] = polyxpoly(S(s).X([1,end-1]),S(s).Y([1,end-1]),E.X,E.Y);
        end
        
        if ~isempty(xi)                 % At least one intersection with the kommune border is found.            
            if length(xi)==1            % Exactly one intersection is found
                if length(xt)==1        % Its only one point inside.
                    tx = [xt,xi,NaN];
                    ty = [yt,yi,NaN];
                    nd = sqrt((xt-xi)^2+(yt-yi)^2);
                else                    % There are several points inside domain.                    
                    % find the nearesty point to the intersection
                    clear d1
                    for nx = 1:length(xt)
                        d1(nx) = sqrt((xi-xt(nx))^2 +(yi-yt(nx))^2);
                    end
                    itp = find(d1==min(d1));
                    % Replace the closest point inside municipality by the
                    % intersection point. For distance calculations.
                    xt(itp) = xi;
                    yt(itp) = yi;    
                    for nx = 1:length(xt)-1
                        d(nx) = sqrt((xt(nx)-xt(nx+1))^2 +(yt(nx)-yt(nx+1))^2);
                    end
                    tx = [xt,NaN];
                    ty = [yt,NaN];
                    nd = nansum(d); 
                    clear d
                end
            else
                % multiple intersections:
                % Straighten shit out!
                % Find the maximum distance of intersections and link subnodes.
                xst = [xi;xt'];
                yst = [yi;yt'];
                nd  = sqrt((min(xst)-max(xst)).^2 +(min(yst)-max(yst)).^2);
            end
            
        else                           % Did not find intersection cant deal with this right now            
            warning(sprintf('\nNo intersection found for RL: Border?\n'))
            nd = (S(s).DISTANCE*1000);
        end
        fra = nd/(S(s).DISTANCE*1000);
        if fra >1.01
            fra = nd/(S(s).DISTANCE)
        end
        if debug_mode
            fprintf('%04.4f / %04.4f (%04.1f%%) in StartGrid \n',DST(s)*fra,DST(s),fra*100)
        end
        clear xi yi tx ty xt yt
    end
    clear in
end

end

