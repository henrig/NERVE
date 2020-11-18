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
function TMout = Emissions_Factor_Mix_in_biofuels(TFout)
%--------------------------------------------------------------------------
% Miljodirektoratet Traffic emission model NERVE:
%
%     FUNCTION :: mix_in_biofuels_NERVE ::
%
%  Function to adjust emission factor for biofuels. Only applied if comp == CO2.
%
% -INPUT  : .xlsx biofules file (biofuelsf) sheet='Biofuels'. Emssion
%  factors per municipality (KEF).
% -OUTPUT : KEF adjusted for biofuels.
%
% 09.03.2018 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global Tyear debug_mode text_div exchfile input
fprintf('\t---------------------------------------------------------------\n')
fprintf('\tin Emissions_Factor_Mix_in_biofuels *\n')
fprintf('\t---------------------------------------------------------------\n')
fprintf('\t---- Subtracting biofuel share from CO2 Emission Factors \n')

sheet           = 'Biofuels';
BioInnblanding  = readtable(input.files.Bio_mix_file,'Sheet',sheet,'ReadVariableNames',true);
fprintf('\tRead Biofuels from %s\n\tSheet ::: %s\n',input.files.Bio_mix_file,sheet)
pos             = find(ismember(BioInnblanding.Properties.VariableNames,sprintf('x%i',Tyear)));
Bio             = table2array(BioInnblanding(:,pos));
fprintf('\t-----------\n')
fprintf('\tBioInnblanding i år %-20i var: \n\t %5.2f %% Totalt \n\t% 5.2f %% i Bensin\n\t %5.2f %% i Diesel\n\t %5.2f %% i Gass\n',Tyear,Bio*100)
fprintf('\t-----------\n')
fprintf('\tAdjusting Emission Factors *\n')

% load exchange file sheet model ID tags SSB_Vehicle_dist
TM       = readtable(input.files.SSB_Vehicle_dist,'Sheet','MODEL');
fuels    = unique(TM.FuelNum);
Vehnames = TFout.Properties.VariableNames;

EFs      = table2array(TFout(:,8:end));
EFheader = TFout.Properties.VariableNames(8:end);

%__________________________________________________________________________
ZeroE = find(TM.FuelNum == 0 | TM.FuelNum == 10| TM.FuelNum == 4);
for i = 1:length(ZeroE)
    idE(i) = find(ismember(EFheader,TM.Name(ZeroE(i))));
    if debug_mode
        fprintf('\tZero EMission, No BIO : %i %s \n',TM.ModelNumber(ZeroE(i)),char(TM.Name(ZeroE(i))))
    end
end
fprintf('\tFound %i Zero Emission Vehicle Groups\n',length(idE))

%__________________________________________________________________________
Petrol   = find(TM.FuelNum == 1 | TM.FuelNum == 6 | TM.FuelNum == 8| TM.FuelNum == 11 |TM.FuelNum == 12| TM.FuelNum == 13);
for i = 1:length(Petrol)
    idP(i) = find(ismember(EFheader,TM.Name(Petrol(i))));
    if isempty(idP)
        fprintf('### MODEL NOT FOUND %s\n',char(TM.Name(Petrol(i))))
    else
        if debug_mode
            fprintf('Petrol MODEL : %i %s \n',TM.ModelNumber(Petrol(i)),char(TM.Name(Petrol(i))));
        end
        EFs(:,idP(i)) = EFs(:,idP(i))*(1-Bio(2));
    end
end
fprintf('\tFound %i Petrol Vehicle Groups\n',length(idP))

%__________________________________________________________________________

Diesel   = find(TM.FuelNum == 2 | TM.FuelNum == 5 | TM.FuelNum == 7 | TM.FuelNum == 9 );
for i = 1:length(Diesel)
    idD(i) = find(ismember(EFheader,TM.Name(Diesel(i))));
    if isempty(idD)
        fprintf('### MODEL NOT FOUND %s\n',char(TM.Name(Diesel(i))))
    else
        if debug_mode
            fprintf('Diesel MODEL : %i %s \n',TM.ModelNumber(Diesel(i)),char(TM.Name(Diesel(i))));
        end
        EFs(:,idD(i)) = EFs(:,idD(i))*(1-Bio(3));
    end
end
fprintf('\tFound %i Diesel Vehicle Groups\n',length(idD))

%__________________________________________________________________________
Gas      = find(TM.FuelNum == 3);
for i = 1:length(Gas)
    idG(i) = find(ismember(EFheader,TM.Name(Gas(i))));
    if isempty(idG)
        fprintf('### MODEL NOT FOUND %s\n',char(TM.Name(Gas(i))))
    else
        if debug_mode
            fprintf('Diesel MODEL : %i %s \n',TM.ModelNumber(Gas(i)),char(TM.Name(Gas(i))));
        end
        EFs(:,idG(i)) = EFs(:,idG(i))*(1-Bio(4));
    end
end
fprintf('\tFound %i Gas Vehicle Groups\n',length(idG))


%__________________________________________________________________________
% 0 1 2 3 4 x 6 7 8 9 10 11 12 13   
unspecified = find( TM.FuelNum >=14);
idU = [];
for i = 1:length(unspecified)
    idU(i) = find(ismember(EFheader,TM.Name(unspecified(i))));
    if isempty(idU)
        fprintf('### MODEL NOT FOUND %s\n',char(TM.Name(Gas(i))))
    else
        if debug_mode
        fprintf('unspecified fuel by MODEL : %i %s idU(i)=%i\n',TM.ModelNumber(unspecified(i)),char(TM.Name(unspecified(i))),idU(i));
        end
        %EFs(:,idU(i)) = EFs(:,idU(i))*(1-Bio(4));
    end
end
fprintf('\tFound %i unspecified Vehicle Groups\n',length(idU))
fprintf('\tTotal %i Vehicle Groups\n',length(idE)+length(idP)+length(idD)+length(idG)+length(idU))


EFold = nanmean(table2array(TFout(:,8:end)));
EFnew = nanmean(EFs);
fprintf('\tAdjusted Petrol by -%4.1f%% of %4.1f%%\n',100*(1-nanmean(EFnew(Petrol))/nanmean(EFold(Petrol))),100*Bio(2))
fprintf('\tAdjusted Diesel by -%4.1f%% of %4.1f%%\n',100*(1-nanmean(EFnew(Diesel))/nanmean(EFold(Diesel))),100*Bio(3))
fprintf('\tAdjusted Gas    by -%4.1f%% of %4.1f%%\n',100*(1-nanmean(EFnew(Gas))/nanmean(EFold(Gas))),100*Bio(4))    
fprintf('\t CONTINUING\n-----------------------------\n')    
TFout(:,8:end) = array2table(EFs);
TMout = TFout;
end

