function [Vehicle_dist] = Vehicle_Distribution_per_Municipality_SSB()
% Module for processing SSB DD. 
global Tyear SSB_Vehicle_dist

fprintf('\t in Vehicle_Distribution_per_Municipality_SSB\n')


% OPT TODO: ADD a test to check for data from municipality and year
Vehicle_data = Vehicle_Distribution_Preprocess_SSB_DD_NV();

% Extract only the traffic from the traffic year and reshape euro and SSB
% individual car type to one class. 
idx    = find(Vehicle_data.D1_yrs == Tyear);
tempNV = reshape(Vehicle_data.YNV(idx,:,:,:),size(Vehicle_data.YNV,2),size(Vehicle_data.YNV,3)*size(Vehicle_data.YNV,4));
tempTD = reshape(Vehicle_data.YTD(idx,:,:,:),size(Vehicle_data.YTD,2),size(Vehicle_data.YTD,3)*size(Vehicle_data.YTD,4));

SSBkomm = Vehicle_data.D2_komm;
% Now use the table(s) to merge the vehicles into the MODEL vehicle types. 
T       = readtable(SSB_Vehicle_dist,'Sheet','MODEL');


% Find the right vehicles for each MODEL group. For DD and NV this is a
% plain sum.
for i=1:height(T)
    idx = str2num(char(split(T.VehicleNumbers_SSB(i),';')));
    modelNV(:,i) = nansum(tempNV(:,idx),2);   
    modelTD(:,i) = nansum(tempTD(:,idx),2);
end

% simple test for data, if seems ok proceed.
a  = sum(modelNV,1);
for i=1:height(T)
   if T.VehicleYears_SSB(i)<a(i)
    fprintf('### error? SSBtot %7i < newMODEL %i\n',round(T.VehicleYears_SSB(i)),round(a(i)))
   end
end

% Recalculate as a fraction of each traffic group L,B,H
idL = T.ClassNum<3;
idH = T.ClassNum>4;
idB = ~idL&~idH;
modelNVfrac = nan(size(modelNV));
modelTDfrac = nan(size(modelTD));
for i=1:size(modelNV,1)
    modelNVfrac(i,idL) = modelNV(i,idL)/nansum(modelNV(i,idL));
    modelNVfrac(i,idH) = modelNV(i,idH)/nansum(modelNV(i,idH));
    modelNVfrac(i,idB) = modelNV(i,idB)/nansum(modelNV(i,idB));

    modelTDfrac(i,idL) = modelTD(i,idL)/nansum(modelTD(i,idL));
    modelTDfrac(i,idH) = modelTD(i,idH)/nansum(modelTD(i,idH));
    modelTDfrac(i,idB) = modelTD(i,idB)/nansum(modelTD(i,idB));
end

% Call the exchange matrix:
[EXkmneNr,TrafficIN,TrafficFROM] = Vehicle_Distribution_Municipal_Traffic_Exchange();

% Calculate based on the traffic exchange what the actual mix of vehicles
% driving in each municipality
modelVdistIN = zeros(size(modelTDfrac));
for i=1:size(TrafficIN)
    % add test to check that i match the right SSB municipality (replace i).
    idkSSB = find(SSBkomm==EXkmneNr(i));
    in = find(TrafficIN(i,:)>0);
    for j = 1:length(in)
           kmm = EXkmneNr(in(j));
           modelVdistIN(idkSSB,:) = modelVdistIN(idkSSB,:) +(modelTDfrac(in(j),:)*TrafficIN(i,in(j))/100);
    end
end

Vehicle_dist.Vdist         = modelVdistIN;
Vehicle_dist.D1_KommNr     = SSBkomm;
Vehicle_dist.D2_Vehicle    = T.Name;
Vehicle_dist.D2_VehicleNum = T.ModelNumber;
Vehicle_dist.TrafficIN     = TrafficIN;
Vehicle_dist.TrafficFROM   = TrafficFROM;
Vehicle_dist.D12_Traffic   = EXkmneNr;


%--------------------------------------------------------------------------
% Statistical output.
    Vehicle_Distribution_Statistical_output(Vehicle_data,T)
end
