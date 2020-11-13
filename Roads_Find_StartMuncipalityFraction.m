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
function wRLinks = Roads_Find_StartMuncipalityFraction(RLinks)
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


Fraction = zeros(size(t_RLinks));
fail = 1;
for i=1:length(t_RLinks)
    startkommune = extractfield(t_RLinks(i),'KOMMS');
    endekommune  = extractfield(t_RLinks(i),'KOMME');
    midtkommune  = extractfield(t_RLinks(i),'KOMM');
    fprintf('%-4i',i)
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
end
KeepT = struct2table(Keep_RLinks);    
t_T   = struct2table(t_RLinks);    

KeepT.IDO(:) = 1; 
t_T.IDO(:)   = Fraction; 
Tble = [KeepT;t_T];
Tble = sortrows(Tble,{'KOMMS'});
wRLinks = table2struct(Tble);
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





