function [Tn] = Vehicle_Distribution_National_HBEFA()
% Function to make 778_sorted, and extract the relavant years distribution
% of vehicles.

global vehicles Tyear

% 778 Vehicles
try % Try this if it exist and reads
    Tout  = readtable(vehicles,'Sheet','778_sorted');
catch
    load('/storage/nilu/Inby/Emission_Group/Emission_Models/HEDGE/Temp/EFA_matrix41_RAW_CO2.mat')
    HasEF = roads.Vehicles;
    T  = readtable(vehicles,'Sheet','all_sorted');
    YW = table2array(T(:,2:end-2));
    
    nYW = NaN(length(HasEF),size(YW,2));
    
    teller = 0;
    found = zeros(size(HasEF));
    pos   = zeros(size(T,1),1);
    for i = 1:height(T)
        idx = find(ismember(HasEF,T.Vehicle(i)));
        if isempty(idx)
            teller = teller+1;
        else
            nYW(idx,:) = YW(i,:);
            found(idx) = i;
            pos(i)     = idx;
        end
    end
    
    if teller > 0
        fprintf('Found Weighted Vehicle with no Emission Factor \n')
    end
    
    
    Tnew = table;
    Tnew.Vehciles = HasEF;
    Tt = array2table(nYW);
    Tout =[Tnew,Tt];
    Tout.Properties.VariableNames = T.Properties.VariableNames(1:end-2);
    %     classSTR = [{'UBus'}, {'Bike'}, {'Scooter'}, {'TT/AT'}, {'RT'}, {'RigidTruck'}, {'PC'},{'MC'},{'LCV'},{'HGV'}.{'Coach'}];
    %     classNUM = [3       , 0       , 0          , 5        , 5     ,5              , 1     ,0     ,2      ,5      .4];
    writetable(Tout,vehicles,'Sheet','778_sorted');
end
yn = sprintf('x%i',Tyear);
Tn = Tout(:,[1,find(ismember(Tout.Properties.VariableNames,yn)),width(Tout)]);

end
