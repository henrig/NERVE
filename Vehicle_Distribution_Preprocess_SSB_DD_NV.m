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
function Vehicle_data = Vehicle_Distribution_Preprocess_SSB_DD_NV()
% This function is made to read and transform the Driving Distance (DD)
% and Number of Vehicles per Municipality supplied by SSB. With
% Norwegian administrative borders changing over time, the timeseries need
% adjustment to the current administrative borders. In addition the format
% is changed into a matrix.
% New data from SSB:
% This preprocessing is primarily to convert the SSB municipalities 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
% 546 possible vehicle types
%  78 nybiltype
%   7 regaar

global SSB_vehicle Kartverket use_temporary_files tfiles

fprintf('\t in Vehicle_Distribution_Preprocess_SSB_DD_NV\n')

if use_temporary_files
    try
        fprintf('Looking for file %s \n',tfiles.CarPark)
        load(tfiles.CarPark)
        return
    catch
        fprintf('#### No temporary file found \n')
    end
end

nybiltype = 78;
regaar    = 7;

% import the SSB data 
T = readtable(SSB_vehicle);


% Read Kartverket data for kommunesammenslåing, two sheets
k1 = readtable(Kartverket,'Sheet','Kartverket_kommuner_2020');
k2 = readtable(Kartverket,'Sheet','Kommuneendring2000_2019');

% Define all the kommunenumbers that appear in the period
kommuner         = k1.Kommunenr_2020;
kommuneNavn      = k1.Kommunenavn_2020;
kommuner2020     = unique(k1.Kommunenr_2020);
for i = 1:length(kommuner2020)
    idx = find(k1.Kommunenr_2020 ==kommuner2020(i));
    kommuneNavn2020(i)     = {k1.Kommunenavn_2020(idx(1))};
end
kommuner2019 = k1.Kommunenr_2019;
kommuner_old = k2.Old_ID;

% First step: Find and replace kommunenumbers that are were changed in or
% before 2019:
komm =  unique(T.kommkode);
for i =1:length(komm)
    idx = T.kommkode==komm(i);
    % First test if the kommune_nr prior to 2020
    idx2 = kommuner_old==komm(i);
    if sum(idx2)==1
        fprintf('Kommune: %04i changed to %04i %s\n',komm(i),kommuner_old(idx2),char(k2.Date(idx2)))
        T.kommkode(idx) = k2.New_ID(idx2);
    elseif sum(idx2) > 1
        fprintf('Kommune: %04i Multiple Name changes to %04i in %s\n',komm(i),kommuner_old(idx2),char(k2.Date(idx2)))
    else
        fprintf('Kommune: %04i Number / Name changed \n',komm(i))
    end
end

% Second step: make the kommunemergers that happened in 2020 
komm  =  unique(T.kommkode);
Ttemp = table;
Tadd  = table;
for i =1:length(komm)
    idx = T.kommkode==komm(i);
    % First test if the kommune_nr exist in 2020
    idx2 = kommuner2019==komm(i);
    if sum(idx2)==1
        fprintf('Kommune: %04i changed to %04i \n',komm(i),kommuner(idx2))
        T.kommkode(idx) = kommuner(idx2);
    elseif sum(idx2) > 1
        a = find(idx2);
        % Split kommuner by even parts
        Ttemp = T(idx,:);
        Ttemp.antallkomm = Ttemp.antallkomm/sum(idx);
        for j =1:length(a)
            fprintf('\tKommune: %04i is split in 2020 \n',kommuner2019(a(j)))
            Ttemp.kommkode = repmat(kommuner(a(j)),height(Ttemp),1);
            Tadd = [Tadd;Ttemp];
        end        
        T = [T(~idx,:);Tadd];
        %Ttemp = repmat(T(idx,:),sum(idx),1);
    end
end

k = 0;
for i =1:length(kommuner)
    idx = T.kommkode==kommuner(i);
    if sum(idx)==0
        fprintf('No match for kommune %04i\n',komm(i))
        k = k+1;
    end
end

komm =  unique(T.kommkode);
for i =1:length(komm)
    idx = kommuner==komm(i);
    if sum(idx)==0
        fprintf('No match for kommune %04i\n',komm(i))
        k = k+1;
    end
end

% % % id1 = T.kommkode==4601;
% % % id2 = T.nybiltype>26;
% % % id3 = T.nybiltype<42;
% % % sum(id1&id2&id3)
% % % Tbergenbuss = T(id1&id2&id3,:)
% % % writetable(Tbergenbuss,'BergenBusesforTØI.xlsx')
% end Kommunesammenlåing




ST = [];
yrs = unique(T.tid);
st = 1;
YNV     = zeros(length(yrs),length(kommuner2020),nybiltype,regaar);
YNV     = zeros(length(yrs),length(kommuner2020),nybiltype,regaar);

c       = zeros(4,1); % counter to keep track of what driving distance is applied to a car
for y = 1:length(yrs)
    Tyear = table;
    idx = T.tid== yrs(y);
    Data = T(idx,:);
    fprintf('Year %i\n',yrs(y))
    
    NV      = zeros(length(kommuner2020),nybiltype,regaar);
    TL      = zeros(size(NV));
    VL      = zeros(size(NV));
    SSB_type= NaN(size(NV));
    SSB_EU  = NaN(size(NV));

    for komm = 1:length(kommuner2020)
         fprintf('\tKommune %i %s\n',kommuner2020(komm),char(kommuneNavn2020{komm}))
        I1 = Data.kommkode==kommuner2020(komm);
        BP = Data(I1,:);
        ST = [ST;[BP.kommkode(1),  BP.fylke(1),  BP.tid(1),  nansum(BP.antallkomm)]];           
        % Remove the extras as a result of merging municipalities where we
        % have double entries.
        remo = zeros(height(BP),1);
        for bp=1:nybiltype
            Ibt = BP.nybiltype==bp;
            for k=1:regaar
                Ieu = BP.regaar==k;
                I=Ibt&Ieu;
                if sum(I)>1
                    idx = find(I);
                    % kommkode;fylke;nybiltype;regaar;snittkomm;antallkomm;snittfylke;antallfylke;snittland;antalland;tid
                    %    5054    17      51      3         NaN       1       36191       17         36297      447     2016
                    if isnan(sum(BP.snittkomm(idx)))
                        gx = ~isnan(BP.snittkomm);
                        idf = find(gx&I);
                        if ~isempty(idf)
                            BP.snittkomm(idx(1))  = nansum(BP.antallkomm(idf).*BP.snittkomm(idf)) / nansum(BP.antallkomm(idf));
                        end
                    else
                        BP.snittkomm(idx(1))  = nansum(BP.antallkomm(idx).*BP.snittkomm(idx)) / nansum(BP.antallkomm(idx));
                    end
                    BP.antallkomm(idx(1)) = nansum(BP.antallkomm(idx));
                    for am = 2:length(idx)
                        remo(idx(am)) = 1;
                    end
                end
            end
        end
        BP = BP(~remo,:);
        
        for bp=1:nybiltype
            Ibt=BP.nybiltype==bp;
            % loop "produksjonsår"
            for k=1:regaar
                Ieu=BP.regaar==k;
                I=Ibt&Ieu;
                % test that a car was found that need a matching
                if sum(I)==1
                    SSB_type(komm,bp,k)= bp;
                    SSB_EU(komm,bp,k)  = 7-k;
                    % Fills in NaN or car number registered value
                    if NV(komm,bp,k)==0 || isnan(NV(komm,bp,k))
                        NV(komm,bp,k)=BP.antallkomm(I);
                    elseif NV(i,bp,k)>=0
                        NV(komm,bp,k)=NV(i,bp,k)+BP.antallkomm(I);
                    else %IF here: Something has gone wrong, in here, only one car per kommune cat.
                        warning('KMNE %i EV %i FUE %i EUR %i N %i BP %i\n',i,k,NV(komm,bp,k), bp)
                    end
                    % Add the driving distance and total driving distance
                    % of the group of cars.
                    % 1) Try first Kommune DD if exist
                    % 2) Apply Fylke DD if exist
                    % 3) Apply National DD if exist
                    % 4) final try give it an average of cars in kommune
                    if ~isnan(BP.snittkomm(I))      %1)
                        TL(komm,bp,k) = BP.antallkomm(I).*BP.snittkomm(I);
                        VL(komm,bp,k) = BP.snittkomm(I);
                        c(1)          = c(1)+BP.antallkomm(I);
                    elseif ~isnan(BP.snittfylke(I))  %2)
                        TL(komm,bp,k) = BP.antallkomm(I).*BP.snittfylke(I);
                        VL(komm,bp,k) = BP.snittfylke(I);
                        c(2)          = c(2)+BP.antallkomm(I);
                    elseif ~isnan(BP.snittland(I))  %3)
                        TL(komm,bp,k) = BP.antallkomm(I).*BP.snittland(I);
                        VL(komm,bp,k) = BP.snittland(I);
                        c(3)          = c(3)+BP.antallkomm(I);
                    else                    %4)
                        %warning(sprintf('Could not match Vehicle: %i %s',bp,char(regaar(k))))
                        c(4)=c(4)+BP.antallkomm(I);
                        TL(komm,bp,k)=BP.antallkomm(I).*nanmedian(BP.snittland(I));
                        VL(komm,bp,k)=nanmedian(BP.snittland(I));
                    end
                elseif sum(I)>1
                    NV(komm,bp,k)=sum(BP.antallkomm(I),1);
                    error(sprintf('found %i vehicles of %i,%i in municipality %04i \n',length(find(I)),bp,k,BP(1,1)))
                else
                    NV(komm,bp,k)=0;
                end
            end
        end
        Numkomm(komm)  = squeeze(nansum(nansum(NV(komm,:,:),2),3));
        Distkomm(komm) = squeeze(nansum(nansum(TL(komm,:,:),2),3));
        
        st = st+1;
    end
        YNV(y,:,:,:) = NV;
        YTD(y,:,:,:) = TL;
        clear NV TL
end
Vehicle_data.YNV = YNV;
Vehicle_data.YTD = YTD;
Vehicle_data.D1_yrs = yrs;
Vehicle_data.D2_komm = kommuner2020;
Vehicle_data.D3_nybiltype = [1:nybiltype];
Vehicle_data.D4_euro      = [regaar-1:-1:0];

if use_temporary_files
   save(tfiles.CarPark, 'Vehicle_data') 
end

end


