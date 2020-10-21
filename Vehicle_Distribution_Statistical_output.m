function Vehicle_Distribution_Statistical_output(Vehicle_data,T)
%--------------------------------------------------------------------------
% Statistical output.
% Tb = table;
% Tb.ModNum  = T.ModelNumber;
% Tb.ModName = T.Name;
% Tb.Class   = str2num(char(T.ClassNum));
% 
% Tb.OsloNV    = modelNV(1,:)';
% Tb.OsloADD   = modelTD(1,:)';
% Tb.OsloPct   = modelTDfrac(1,:)'*100;
% Tb.OsloPctEX = modelVdistIN(1,:)'*100;
% writetable(T,'National_SSB_Vehicle_Number_Distribution.xlsx','Sheet','MunicipalStats_VehicleDistribution')
% National Number of vehicles
National = squeeze(sum(Vehicle_data.YNV,2));
k=1;
National2D = zeros(size(National,2)*size(National,3),size(National,1));
Ts = table;
for k = 1:size(National,1)
    t = 1;
    for j = 1:size(National,3)
        for i = 1:size(National,2)
            Ts.Num(t)  = t;
            Ts.Euro(t) = Vehicle_data.D4_euro(j);
            Ts.Kategori(t)  = Vehicle_data.D3_nybiltype(i);
            Ts.Var1(t) = National(k,i,j);
            t=t+1;
        end
    end
    vn = find(ismember(Ts.Properties.VariableNames,'Var1'));
    Ts.Properties.VariableNames(vn)= {sprintf('x%i',Vehicle_data.D1_yrs(k))};
end
writetable(T,'National_SSB_Vehicle_Number_Distribution.xlsx','Sheet','NationalNumberVehicles')
% National Driving Distance
National = squeeze(sum(Vehicle_data.YTD,2));
k=1;
National2D = zeros(size(National,2)*size(National,3),size(National,1));
Ts2 = table;
for k = 1:size(National,1)
    t = 1;
    for j = 1:size(National,3)
        for i = 1:size(National,2)
            Ts2.Num(t)  = t;
            Ts2.Euro(t) = Vehicle_data.D4_euro(j);
            Ts2.Kategori(t)  = Vehicle_data.D3_nybiltype(i);
            Ts2.Var1(t) = National(k,i,j);
            t=t+1;
        end
    end
    vn = find(ismember(Ts2.Properties.VariableNames,'Var1'));
    Ts2.Properties.VariableNames(vn)= {sprintf('x%i',Vehicle_data.D1_yrs(k))};
end
writetable(T,'National_SSB_Vehicle_Number_Distribution.xlsx','Sheet','NationalTotalDistance')
end
