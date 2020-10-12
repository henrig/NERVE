function OLinks =  Roads_Remove_NoTrafficRoads(RLinks)

global remove_NoTrafficRoads
if remove_NoTrafficRoads == 0
    return
end
   
fieldN = fieldnames(RLinks);
t1 = contains(upper(fieldN),'ADT');
t2 = contains(upper(fieldN),'LETTE');

t = find(t1|t2);
if isempty(t)
    fprintf('### No ADT fields found ###\n')
    return
end

traff = zeros(1,length(RLinks));
for i = 1:length(t)
    if t(i)
        f = extractfield(RLinks,char(fieldN(t(i))));
        traff = traff+f;
        fprintf('%10s  %6i : %6i : %6i \n',char(fieldN(t(i))),length((find(f>0))),length((find(traff>0))),length(RLinks))
    end
end

OLinks = RLinks(traff>0);

end