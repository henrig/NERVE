function [Sout] = Roads_Roadsplitter(S)
% The main aim is to reduce the amounts of subnodes necessary to describe
% each road. Straighten or split roads by their subnodes.

global maxSegments minLength angle_change


if mean(extractfield(S,'DISTANCE')) < 10
    Length_unit = 'km';
    minLength = minLength/1000;
    fprintf('Length Unit assumed = %s minLength = %i\n',Length_unit,minLength)
else
    Length_unit = 'm';
    fprintf('Length Unit assumed = %s minLength = %i\n',Length_unit,minLength)
    
end

segments  = zeros(size(S));
subsegm   = length(extractfield(S,'X'))-(length(S)*3);

potential_Roads = length(subsegm)-length(segments);
fprintf('Total roads: %i \n',length(S))
fprintf('Total subnodes: %i Potential Roads = \n',subsegm,potential_Roads )

Sout   = [];
teller =0;
for i=1:length(S)
    if rem(i,500)==0; fprintf('Splitting Roads: @ %i of %i added %i Links  to system\n',i,length(S),teller);end
    V = S(i);
    % splitted = 0;
    % while splitted < max_split
    
    segments(i) = length(V.X)-2;
    if segments(i) > maxSegments

        if S(i).DISTANCE < minLength
            % Segment is short and will be made straight
            V.X = [xt(1),xt(end-1:end)];
            V.Y = [yt(1),yt(end-1:end)];
            Sout = [Sout,V];
        else
            % Segment is long with multiple geographical details.
            xt = V.X;
            yt = V.Y;
            v1 = [xt(1),yt(1),0];
            v2 = [xt(end-1),yt(end-1),0];
            % Test here to remove subsegment xt and yt points.
            for m=1:length(xt)
                pt =  [xt(m),yt(m),0];
                a = v1 - v2;
                b = pt - v2;
                d(m) = norm(cross(a,b)) / norm(a);
            end
            delta_d(1) = 0;
            delta_d(2:length(d)) = d(1:end-1)-d(2:end);
            
            % Keep only the segments wich cause a large change in direction
            keepes  = find(abs(delta_d)>angle_change);
            
            % Keep at least one segment to display shape.
            if isempty(keepes) && nanmax(d)>angle_change
                keepes = find(d == nanmax(d));
            end
            potential_Roads = potential_Roads+length(keepes)-(length(xt)-2);
            xt = [xt(1),xt(keepes),xt(end-1:end)];
            yt = [yt(1),yt(keepes),yt(end-1:end)];
            clear d delta_d
            
            TOTdistance = sqrt((xt(1)-xt(end-1))^2+(yt(1)-yt(end-1))^2);
            % This part cacluates the new emissions fields.
            V2 =[];
            for nl = 2:length(xt)-1
                Segment_Length = sqrt((xt(nl-1)-xt(nl))^2+(yt(nl-1)-yt(nl))^2);
                NewLink = V;
                NewLink.X = [xt(nl-1),xt(nl),NaN];
                NewLink.Y = [yt(nl-1),yt(nl),NaN];
                NewLink.DISTANCE = Segment_Length;
                if ismember(Length_unit,'km')
                    NewLink.DISTANCE = Segment_Length/1000;
                else
                    NewLink.DISTANCE = Segment_Length;
                end
                % Should ideally make a test for each field.
                % if isfield(NewLink,'EM_FC')
                NewLink.EM_FC    = V.EM_FC*Segment_Length/TOTdistance;
                NewLink.EM_FC_MJ = V.EM_FC_MJ*Segment_Length/TOTdistance;
                NewLink.EM_CH4   = V.EM_CH4*Segment_Length/TOTdistance;
                NewLink.EM_BC    = V.EM_BC*Segment_Length/TOTdistance;
                NewLink.EM_PM    = V.EM_PM*Segment_Length/TOTdistance;
                NewLink.EM_HC    = V.EM_HC*Segment_Length/TOTdistance;
                NewLink.EM_CO    = V.EM_CO*Segment_Length/TOTdistance;
                NewLink.EM_NOx   = V.EM_NOx*Segment_Length/TOTdistance;
                NewLink.EM_Be    = V.EM_Be*Segment_Length/TOTdistance;
                NewLink.EM_NMHC  = V.EM_NMHC*Segment_Length/TOTdistance;
                NewLink.EM_NO2   = V.EM_NO2*Segment_Length/TOTdistance;
                NewLink.EM_PN    = V.EM_PN*Segment_Length/TOTdistance;
                NewLink.EM_CO2   = V.EM_CO2*Segment_Length/TOTdistance;
                NewLink.EM_TW10  = V.EM_TW10*Segment_Length/TOTdistance;
                V2 =[V2,NewLink];
            end
            Sout = [Sout,V2]; % Add new buch of roads V2 to roadnetwork instead.
            teller = teller+length(V2)-1;
        end
    else
        Sout = [Sout,V];
    end
end

subsegm2   = length(extractfield(Sout,'X'))-(length(Sout)*3);
fprintf('Total roads: %i \n',length(Sout))
fprintf('Total subnodes: %i \n',subsegm2)
end

