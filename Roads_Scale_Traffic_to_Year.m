function [Sn] = Roads_Scale_Traffic_to_Year(RLinks)
% Scale_Road_Traffic_Year function takes a traffic file and scales it
% from the year (Ryear) to a new year (Tyear) based on the csv file. It
% does scaling at the municipal basis and treats each municipality
% independently.

% It also makes a lighter version of the Roads that only keep the
% information that is required to do emission calculations. 

global Ryear Tyear traff_years traffile tfold

a = strfind(traffile,'/');
if isempty(a)
    ofile = traffile;
else
    fname = traffile(a(end)+1:end);
    ofile = strcat(tfold,fname,sprintf('HEDGE_%i',Tyear));
end

try
    Sn = shaperead(ofile);
    return
catch
end

fprintf('in Scale_Road_Traffic_Year  \n')

fprintf('Reading Annual file :\n %s \n',traff_years)
fprintf('To scale Traffic from %i to %i \n',Ryear,Tyear)

T       = readtable(traff_years,'ReadVariableNames',1);
fldList = T.Properties.VariableNames;
f       = find(ismember(fldList,sprintf('x%i',Ryear)));
fy      = find(ismember(fldList,sprintf('x%i',Tyear)));


Scale         = table;
Scale.KOMM    = T.Kommune;
Scale.ToYearL = table2array(T(:,fy))./table2array(T(:,f));


k = unique(extractfield(RLinks,'KOMM'));


Tlink   = struct2table(RLinks);


found_err = 0;
L_adt = zeros(length(RLinks),1);
H_adt = zeros(length(RLinks),1);
B_adt = zeros(length(RLinks),1);
for i = 1:length(k)
    idk = find(Scale.KOMM==k(i));
    if ~isempty(idk)
        idr = find(Tlink.KOMM == k(i));
        if ~isempty(idr)
            L_adt(idr) = Tlink.LETTE_BILE(idr)*Scale.ToYearL(idk);
            H_adt(idr) = Tlink.GODS_ADT(idr)*Scale.ToYearL(idk);
            B_adt(idr) = Tlink.KOLL_ADT(idr)*Scale.ToYearL(idk);
        else
            fprintf('Missinig roads for Municipality %i\n', k(i))
            found_err = 1;
        end
        
    else
        fprintf('Missinig scale for Municipality %i\n', k(i))
        found_err = 1;
    end
end
if found_err
    warning('could not complete calculations for all roads')
end 
fprintf('Scaled all Municipalities \n')

Tlink.L_ADT = L_adt;
Tlink.H_ADT = H_adt;
Tlink.B_ADT = B_adt;

Listfields = [{'Geometry'},{'X'},{'Y'},{'BoundingBox'},{'DISTANCE'},{'KOMM'},{'SPEED'},{'SLOPE'},{'URBAN'},{'RUSH_DELAY'},{'HBEFA_EQIV'},{'WIDTH'},{'N_LANES'},{'L_ADT'},{'H_ADT'},{'B_ADT'}];
for i=1:length(Listfields)
   f(i)  = find(ismember(Tlink.Properties.VariableNames,Listfields(i)));
   %Listfields(i)
end
Tlink.Properties.VariableNames(f(14)) = {sprintf('L_ADT%04i',Tyear)};
Tlink.Properties.VariableNames(f(15)) = {sprintf('H_ADT%04i',Tyear)};
Tlink.Properties.VariableNames(f(16)) = {sprintf('B_ADT%04i',Tyear)};


fprintf('Structuring Roads \n')
Sn = table2struct(Tlink(:,f));
end







