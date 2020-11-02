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
    Ttemp          = table;
    for j = 1:size(National,3)
        for i = 1:size(National,2)
            if k==1
                Ttemp          = table;
                Ttemp.Num      = t;
                Ttemp.Euro     = Vehicle_data.D4_euro(j);
                Ttemp.Kategori = Vehicle_data.D3_nybiltype(i);
                Ttemp.Var1     = National(k,i,j);
                Ts             = [Ts;Ttemp];
                t              = t+1;
            else
                Ts.Var1(t)  = National(k,i,j);
                t              = t+1;
            end
        end
    end
    vn = find(ismember(Ts.Properties.VariableNames,'Var1'));
    Ts.Properties.VariableNames(vn)= {sprintf('x%i',Vehicle_data.D1_yrs(k))};
end
file = 'National_SSB_Vehicle_Number_Distribution.xlsx';
writetable(T,file,'Sheet','NationalNumberVehicles')
fprintf('Wrote National Vehicle stats file :\n%s\n',file)
% National Driving Distance
National = squeeze(sum(Vehicle_data.YTD,2));
k=1;
National2D = zeros(size(National,2)*size(National,3),size(National,1));
Ts2 = table;
for k = 1:size(National,1)
    t = 1;
    Ttemp          = table;
    for j = 1:size(National,3)
        for i = 1:size(National,2)
            if k==1
                Ttemp          = table;
                Ttemp.Num      = t;
                Ttemp.Euro     = Vehicle_data.D4_euro(j);
                Ttemp.Kategori = Vehicle_data.D3_nybiltype(i);
                Ttemp.Var1     = National(k,i,j);
                Ts2 = [Ts2;Ttemp];
                t=t+1;
            else
                Ts2.Var1(t)  = National(k,i,j);
                t              = t+1;
            end
        end
    end
    vn = find(ismember(Ts2.Properties.VariableNames,'Var1'));
    Ts2.Properties.VariableNames(vn)= {sprintf('x%i',Vehicle_data.D1_yrs(k))};
end
writetable(T,'National_SSB_Vehicle_Number_Distribution.xlsx','Sheet','NationalTotalDistance')
end
