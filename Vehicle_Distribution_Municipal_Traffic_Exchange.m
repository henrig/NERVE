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
function [kmne,TrafficIN,TrafficFROM] = Vehicle_Distribution_Municipal_Traffic_Exchange()
%--------------------------------------------------------------------------
% -INPUT  : .xlsx Exchange matrix file (RTMpexch) sheet='Biofuels'. Emssion
%  factors per municipality (KEF).
% -OUTPUT : KEF adjusted for biofuels.
%
% 09.03.2018 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global tfold Tyear KmneNr KmneNavn debug_mode text_div

global traff_exchange traff_exchange_sh

fprintf('%s\n',text_div)
fprintf('in Vehicle_Distribution_Municipal_Traffic_Exchange *\n')
fprintf('%s\n',text_div)
fprintf('---  Finding Origin of Traffic in each municipality \n')

%RTMsheet = '2020_Utveksling';

fprintf('Municipal_Traffic_Exchange_NERVE \n ... reading sheet "%s" from file \n   %s\n',traff_exchange_sh, traff_exchange)
Mex  = readtable(traff_exchange,'Sheet',traff_exchange_sh);
Mex.Properties.VariableNames(1)={'KommNr'};
kmne = Mex.KommNr;
exchange= table2array(Mex(:,2:end));

% Debug the exchange Matrix.

% the assumption is that all kmnes must have at least some traffic in their
% own municipality. Set the lower threshold to 1km. The model will npt
% allow holes in the data
Test=diag(exchange);
I=find(Test==0);
if ~isempty(I)
    warning('There are Municipalities that have no traffic on their own roads.')
    for i=1:length(I)
    fprintf('%04i_%s\n',KmneNr(I(i)),char(KmneNavn(I(i))))
    end
    warning('Setting Driving Distance to 1 km for these kommunes')    
    exchange(I,I)=1;
end

% Make the caluclations
for komm=1:length(kmne)
      DD_FROM(komm)       = sum(exchange(:,komm));
      DD_IN(komm)         = sum(exchange(komm,:));
      TrafficIN(komm,:)   = 100*exchange(komm,:)/DD_IN(komm);
      TrafficFROM(:,komm) = 100*exchange(:,komm)/DD_FROM(komm);
end
InFrom = diag(TrafficFROM);
InIn   = diag(TrafficIN);

if debug_mode
     If = find(InFrom>99.9);
    fprintf('Municipality with 100%% Traffic From their own muncipality\n')
    for i=1:length(If)
     fprintf('%04i_%s\n',kmne(If(i)),char(KmneNavn(If(i))))
    end
    Ii = find(InIn>99.9);
    fprintf('Municipality with 100%% Traffic In their own muncipality\n')
    for i=1:length(Ii)
     fprintf('%04i_%s\n',KmneNr(Ii(i)),char(KmneNavn(Ii(i))))
    end
    fprintf('***** Test Oslo \n\n')
    fprintf('%04i_%s %% Traffic %4.1f %4.1f\n',KmneNr(1),char(KmneNavn(1)),InIn(1),InFrom(1))
    fprintf('***** Test Oslo \n\n')
    fprintf('%04i_%s  %% Traffic On roads in Oslo With Oslo Registered Cars %4.1f %%\n',KmneNr(41),char(KmneNavn(41)),InIn(41))
    fprintf('%04i_%s  Oslo Registered Cars %%  On roads in Oslo      %4.1f %%\n',KmneNr(41),char(KmneNavn(41)),InFrom(41))
end


save(sprintf('%sMunicipal_Traffic_Exchange_%04i.mat',tfold,Tyear),'kmne','TrafficIN','TrafficFROM','InIn','InFrom')

end

