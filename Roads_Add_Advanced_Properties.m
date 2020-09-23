function [Sn] = Roads_Add_Advanced_Properties(RLinks)
    % Convert Precentage into HBEAF interpretable 0,2,4,6 pct slopes
    st = min(abs(round(50*extractfield(RLinks,'STIGNING_P'))*2),6);
    fprintf('* in Add_Road_Advanced_Properties  *\n')

    % Calculate rush hour traffic speed and delay
    di = extractfield(RLinks,'DISTANCE');
    kt = extractfield(RLinks,'FM_TIME');
    km = extractfield(RLinks,'KO_MORGEN');
    ke = extractfield(RLinks,'KO_ETTERM');
    ks = extractfield(RLinks,'FM_SPEED');
    
    % The time traffic that is delayed is calculated by the mean of the 
    % morning and evening congestion
    RUSH_DELAY                    = round(100*(kt+mean([ke;km]))./kt-100);
    RUSH_DELAY(isnan(RUSH_DELAY)) = 0;
    RUSH_SPEED                    = round(di./((kt.*(1+RUSH_DELAY/100)/60)));
    RUSH_SPEED(isnan(RUSH_SPEED)) = 0;
    RUSH_SPEED(RUSH_SPEED>kt)     = kt(RUSH_SPEED>kt);

    T = struct2table(RLinks);
    T.RUSH_DELAY = RUSH_DELAY';
    T.RUSH_SPEED = RUSH_SPEED';
    T.SLOPE      = st';
    Sn = table2struct(T);
    fprintf('Setfields\n  -- RUSH_DELAY\n  -- RUSH_SPEED\n  -- SLOPE ...on all links \n')
end
