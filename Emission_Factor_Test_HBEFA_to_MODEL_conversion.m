function Emission_Factor_Test_HBEFA_to_MODEL_conversion(com,ConVfile)
global EFA EF_AVG comps roads
fprintf('\tin Emission_Factor_Test_HBEFA_to_MODEL_conversion\n')
Traff = table;
k = 1;
for roa = 1:size(EF_AVG,2)
    for spd = 1:size(EF_AVG,3)
        for dec = 1:size(EF_AVG,4)
            for cog = 1:size(EF_AVG,5)
                for urb = 1:size(EF_AVG,6)
                    Tb = table;
                    Tb.name        = {sprintf('%s/%s-%i/%i%%/%s',char(roads.RoadEnv(urb)),char(roads.RoadType(roa)),roads.RoadSpeeds(spd),roads.RoadGradient(dec),char(roads.Congestion(cog)))};
                    Tb.roadTypeID  = roa;
                    Tb.roadType    = roads.RoadType(roa);
                    Tb.roadSpeedID = spd;
                    Tb.roadSpeed   = roads.RoadSpeeds(spd);
                    Tb.roadGradID  = dec;
                    Tb.roadGrad    = roads.RoadGradient(dec);
                    Tb.roadCongID  = cog;
                    Tb.roadCong    = roads.Congestion(cog);
                    Tb.roadCongID  = urb;
                    Tb.roadCong    = roads.RoadEnv(urb);
                    VEH = zeros(1,size(EF_AVG,1));
                    VEHp = zeros(1,size(EF_AVG,1));
                    for veh = 1:size(EF_AVG,1)
                        if ~isnan(EF_AVG(veh,roa,spd,dec,cog,urb))
                            VEH(veh) = VEH(veh)+1;
                        end
                        if EF_AVG(veh,roa,spd,dec,cog,urb)>0
                            VEHp(veh) = VEHp(veh)+1;
                        end
                    end
                    nums1(k) = sum(VEH);
                    nums2(k) = sum(VEHp);
                    k = k+1;
                    Traff = [Traff;Tb];
                end
            end
        end
    end
end
Traff.NumVehEF_HBEFA  = nums1';
Traff.NumVehEFz_HBEFA = nums2';

Traff2 = table;
k = 1;
for roa = 1:size(EFA,2)
    for spd = 1:size(EFA,3)
        for dec = 1:size(EFA,4)
            for cog = 1:size(EFA,5)
                for urb = 1:size(EFA,6)
                    Traff2.name(k)        = {sprintf('%s/%s-%i/%i%%/%s',char(roads.RoadEnv(urb)),char(roads.RoadType(roa)),roads.RoadSpeeds(spd),roads.RoadGradient(dec),char(roads.Congestion(cog)))};
                    VEH = zeros(1,size(EFA,1));
                    VEHp = zeros(1,size(EFA,1));
                    for veh = 1:size(EFA,1)
                        if ~isnan(EFA(veh,roa,spd,dec,cog,urb))
                            VEH(veh) = VEH(veh)+1;
                        end
                        if EFA(veh,roa,spd,dec,cog,urb)>0
                            VEHp(veh) = VEHp(veh)+1;
                        end
                    end
                    nums1(k) = sum(VEH);
                    nums2(k) = sum(VEHp);
                    k = k+1;
                end; end; end; end; end;
Traff2.NumVehEF_MODEL  = nums1';
Traff2.NumVehEFz_MODEL = nums2';

errs = 0;
for i=1:height(Traff)
    idx = find(ismember(Traff2.name,Traff.name(i)));
    if (Traff.NumVehEF_HBEFA(i))>0 && (Traff2.NumVehEF_MODEL(idx))>0
        
    elseif (Traff.NumVehEF_HBEFA(i))>0
        errs = errs +1;
    else
    end
end
if errs >0
    fprintf('#### SHEET CONVERISION ERROR ####\n')
else
    fprintf('CONVERISION TESTED OK.\n')
end
Traff = [Traff,Traff2(:,2:end)];
fprintf('HBEFA : %i / %i / %i \n',round(min(Traff.NumVehEF_HBEFA)),round(mean(Traff.NumVehEF_HBEFA)),round(max(Traff.NumVehEF_HBEFA)) )
fprintf('MODEL : %i / %i / %i \n',round(min(Traff.NumVehEF_MODEL)),round(mean(Traff.NumVehEF_MODEL)),round(max(Traff.NumVehEF_MODEL)) )


Traff = Traff(Traff.NumVehEF_MODEL>0|Traff.NumVehEF_HBEFA>0,:);
writetable(Traff,ConVfile,'Sheet',sprintf('%s',char(comps(com))))
end