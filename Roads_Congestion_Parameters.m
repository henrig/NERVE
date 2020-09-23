function Sout = Roads_Congestion_Parameters(RLinks)

fprintf('\n in fix_KAPTID_to_KOTID \n...\n')

mrg = zeros(size(RLinks));
% Morgen
if isfield(RLinks,'KAPTID_06_')
    a(:,1) = extractfield(RLinks,'KAPTID_06_');
    mrg = max([mrg,a],[],2);
end
if isfield(RLinks,'KAPTID_07_')
    a(:,1) = extractfield(RLinks,'KAPTID_07_');
    mrg = max([mrg,a],[],2);
end

if isfield(RLinks,'KAPTID_08_')
    a(:,1) = extractfield(RLinks,'KAPTID_08_');
    mrg = max([mrg,a],[],2);
end

% Ettermiddag
kvl = zeros(size(RLinks));
if isfield(RLinks,'KAPTID_15_')
    a(:,1) = extractfield(RLinks,'KAPTID_15_');
    kvl = max([kvl,a],[],2);
end

if isfield(RLinks,'KAPTID_06_')
    a(:,1) = extractfield(RLinks,'KAPTID_16_');
    kvl = max([kvl,a],[],2);
end

if isfield(RLinks,'KAPTID_06_')
    a(:,1) = extractfield(RLinks,'KAPTID_17_');
    kvl = max([kvl,a],[],2);
end
T = struct2table(RLinks);
T.KO_MORGEN = mrg;
T.KO_ETTERM = kvl;
fprintf('Structuring Roads...\n')
Sout = table2struct(T);
fprintf('Setfields\n')
fprintf('  -- KO_ETTERM\n')
fprintf('  -- KO_MORGEN ...on all links\n') 
end


