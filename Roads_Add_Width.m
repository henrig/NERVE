function Sout = Roads_Add_Width(RLinks)
% Dekkebredde is top in hieararchy (use if it is defined:)
DB   = extractfield(RLinks,'DEKKEBREDD');
has  = find(DB);

fprintf('found %i /%i RLinks with "DEKKEBREDDE>0"\n',length(has), length(RLinks))
need = find(~DB);

LA = extractfield(RLinks,'LANES');
a = find(isempty(LA));

if isempty(a)
    fprintf('found %i /%i RLinks with "LANES>0"\n',length(LA), length(RLinks))
end

for i=1:length(RLinks)
    nlanes(i) = length(strfind(char(LA(i)),','))+1;
end

% setting WIDTH
for i=1:length(RLinks)
    if ismember(i,has)
       width(i) = DB(i); 
    else        
        type = RLinks(i).VK;
        if ismember(type,'E') 
               width(i) = min([nlanes(i),8])*3.1; 
        elseif ismember(type,'R') || ismember(type,'F')
               width(i) = min([nlanes(i),8])*2.8;         
        elseif ismember(type,'P') || ismember(type,'K')
               width(i) = min([nlanes(i),8])*2.6; 
        else
               width(i) = min([nlanes(i),8])*2.6;             
        end
    end
end
T = struct2table(RLinks);
T.WIDTH   = width';
T.N_LANES = nlanes';
fprintf('Setfields\n')
fprintf('  -- WIDTH\n')
fprintf('  -- N_LANES ...on all links\n')
fprintf('Structuring Roads \n')
Sout = table2struct(T);
end
