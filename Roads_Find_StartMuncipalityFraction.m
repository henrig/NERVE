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
function Sn = Roads_Find_StartMuncipalityFraction(RLinks)
global  Komm_shape 
fprintf('---------------------------------------------------------------\n')
fprintf('in Roads_Find_StartMuncipalityFraction   *\n')
fprintf('---------------------------------------------------------------\n')

Ks = shaperead(Komm_shape);

startkommune = extractfield(RLinks,'KOMMS');
endekommune  = extractfield(RLinks,'KOMME');
midtkommune  = extractfield(RLinks,'KOMM');


komVeg       = nan(size(startkommune));
komVeg(startkommune == endekommune&(startkommune == midtkommune))=1;
komVeg(startkommune ~= endekommune) = 0;
komVeg(startkommune ~= midtkommune) = 0;

fprintf('----------------------------------------\n')
fprintf('Roads internal in    municipalities %7i \n',nansum(komVeg))
fprintf('Roads shared between municipalities %7i \n',length(find(komVeg==0)))
fprintf('Roads   without      municipality   %7i \n',length(find(isnan(komVeg))))
fprintf('----------------------------------------\n')


Keep_RLinks = RLinks(find(komVeg));
t_RLinks    = RLinks(komVeg==0);

Komms       = extractfield(Ks,'KOMMUNENUM');



% First deal with roads that are continouus.  
for i=1:length(t_RLinks)
    startkommune = extractfield(t_RLinks(i),'KOMMS');
    endekommune  = extractfield(t_RLinks(i),'KOMME');
    midtkommune  = extractfield(t_RLinks(i),'KOMM');
    xs  = extractfield(t_RLinks(i),'X');
    individual(i)  = length(find(isnan(xs)));
    subsegments(i) = length(find(~isnan(xs)))-1;     
end

idx = individual==1 & subsegments==1;
Straight_single_roads = find(idx);
Curvy_single_roads    = find(individual==1 & ~idx);
Complex_roads         = find(individual>1 & subsegments~=1);

fprintf('----------------------------------------\n')
fprintf('Straight_single_roads municipalities %7i \n',length(Straight_single_roads))
fprintf('Curvy_single_roads    municipalities %7i \n',length(Curvy_single_roads))
fprintf('Complex_roads         municipality   %7i \n',length(Complex_roads))
fprintf('Unexplained roads     municipality   %7i \n',length(find(komVeg==0))-length(Complex_roads)-length(Curvy_single_roads)-length(Straight_single_roads) )
fprintf('----------------------------------------\n')

RDcomplexity(Straight_single_roads) = 1;
RDcomplexity(Curvy_single_roads)    = 2;
RDcomplexity(Complex_roads)         = 3;




% keep a backup ht = t_RLinks;

% trim some of the subnodes off the roads a bit so they are not too long,
% as the test can run out of memory.
for i = 1:length(t_RLinks)    
    if(RDcomplexity(i)==2 && subsegments(i)>500)
        xs  = extractfield(t_RLinks(i),'X');
        ys  = extractfield(t_RLinks(i),'Y');
        lxy = length(xs);
        startx = xs(1);
        starty = ys(1);
        stopx  = xs(end-1);
        stopy  = ys(end-1);
        while lxy > 450
            segD(1) = sqrt((xs(1)-xs(2))^2 +(ys(1)-ys(2))^2);
            for r = 2:length(xs)
                segD(r) = sqrt((xs(r-1)-xs(r))^2 +(ys(r-1)-ys(r))^2);
            end
            idx = find(segD <= min(segD)+1e-0);
            segD(idx) = NaN;
            idx2 = ~isnan(segD);
            xs = xs(idx2);
            ys = ys(idx2);
            lxy = length(xs);
            clear segD
        end
        newx = [startx,xs,stopx,NaN];
        newy = [starty,ys,stopy,NaN];
        t_RLinks(i).X = newx;
        t_RLinks(i).Y = newy;
    end
end

fail = 1;
Fraction    = zeros(size(t_RLinks));
for i=1:length(t_RLinks)
    startkommune = extractfield(t_RLinks(i),'KOMMS');
    endekommune  = extractfield(t_RLinks(i),'KOMME');
    midtkommune  = extractfield(t_RLinks(i),'KOMM');
    fprintf('%-5i',i)
    
    if ~isnan(startkommune)
        pla = find(Komms == startkommune);
        for p=1:length(pla)
            fra = Roads_inDomain_Roads(t_RLinks(i),Ks(pla(p)));
            Fraction(i) =max([Fraction(i),fra]);
        end
    else
        uKomm = unique([startkommune, endekommune, midtkommune]);
        ik = find(uKomm>0);        
        if length(ik)==1
            t_RLinks(i).KOMM  = uKomm(ik);
            t_RLinks(i).KOMMS = uKomm(ik);
            t_RLinks(i).KOMME = uKomm(ik);
            Fraction(i) = 1;
            fprintf(' ###%i  \n',fail)
            fail = fail+1;
        else
            fprintf(' ###%i  erraneious road @ %i\n',fail,i)
            t_RLinks(i)
            Fraction(i) = 1;
            fail = fail+1;
        end
    end
    if rem(i,50)==0; fprintf('\n'); end
end
KeepT = struct2table(Keep_RLinks);    
t_T   = struct2table(t_RLinks);    

KeepT.IDO(:) = 1; 
t_T.IDO(:)   = Fraction; 
Tble = [KeepT;t_T];
Tble = sortrows(Tble,{'KOMMS'});
Sn = table2struct(Tble);
end

% startkommune = extractfield(RLinks,'KOMMS');
% RLinks  = shaperead('HEDGE_RL_2019_Norway02Dec2019_temp');
% EMLinks = shaperead('HEDGE_Emissions_2019');
% Ks    = shaperead('Kommuner2020');
% 
% Komms   = extractfield(Komm,'KOMMUNENUM');
% uKomm   = unique(Komms);
% for k = 1:length(uKomm)
%     idx    = find(Komms == uKomm(k));
%     Kshape = Komm(idx);        
% end




% pla = find(Komms == 301)
%    close all
%     mapshow(Ks(pla))
%     hold on
%     mapshow(t_RLinks(i))
%     scatter(t_RLinks(i).X,t_RLinks(i).Y)
% mapshow(t_RLinks(33))





