function [Sn] = Roads_Add_Municipality(RLinks,Komm_shape)
%-------------------------------------------------------------------------
% Miljodirektoratet Traffic emission model NERVE:
%
%     FUNCTION :: HEDGE_Add_ROAD_Municipality ::
%
% This function looks for and checks the national traffic file. It makes
% sure that it has the required fields and adds (some) fields if they are
% deducable from the input-traffic. As there is a lot of uneccesary data
% and variable number and properties of the roads, properties and values
% with each iteration of the traffic input file, this is spesifically
% designed to treat the ones we have recieved. As the shapefiles recieved
% are unique in design, and different every time this function mainly deals
% with known problems. Field names for example are treated by name and
% may not be same in future extractions, and a change will have to be made
% in the script.
%
% -INPUT  : this function reads the road link shapes per municipality to the
%           temporary folder (tpath).
% -OUTPUT : this function adds  properties to the road link shapes per
%           municipality in the temporary folder (tpath).
%
% 20.09.2018 -Henrik Grythe
% Kjeller NILU
%-------------------------------------------------------------------------

fprintf('* call Add_ROAD_Municipality   *\n')
fprintf('Adding kommune field based on geo file: \n%s\n',Komm_shape)
Ks      = shaperead(Komm_shape);
prj     = HEDGE_read_projection(Komm_shape);

% Adds a Roads Midpoint based on the starting point and ending point of
% Road
if ~isfield(RLinks,'MIDPOINT_X') || ~isfield(RLinks,'MIDPOINT_Y')
    for i=1:length(RLinks)
        x=RLinks(i).X;
        y=RLinks(i).Y;
        MIDPOINT_X(i)=(x(1)+x(end-1))/2;
        MIDPOINT_Y(i)=(y(1)+y(end-1))/2;
    end
else
     MIDPOINT_X = extractfield(RLinks,'MIDPOINT_X');
     MIDPOINT_Y = extractfield(RLinks,'MIDPOINT_Y');
end

teller = 0;
for i=1:length(Ks)
    in = inpolygon(MIDPOINT_X,MIDPOINT_Y,Ks(i).X,Ks(i).Y);
    KOMM(in) = extractfield(Ks(i),'KOMMUNENUM');
    fprintf('%03i Kommune %04i_%s has %i roads \n',i,extractfield(Ks(i),'KOMMUNENUM'),char(extractfield(Ks(i),'NAVN')),sum(in))
    has(in)=1; clear in
    fprintf('Placed %i/%i = %4.1f%% \n',sum(has),length(RLinks), 100*sum(has)/length(RLinks))
end


fprintf('Have: %i of %i Roadlinks "KOMM", missing: %i \n',sum(has),length(RLinks),length(RLinks)-sum(has))
% try midpoint start point and end point of road, which commune it
% will be decided.
missing=find(has==0);
if ~isempty(missing)
    % If midpoint doesnt work, do it on the starting point
    for i=1:length(missing)
        x=extractfield(RLinks(missing(i)),'X');
        y=extractfield(RLinks(missing(i)),'Y');
        for j=1:length(Ks)
            in = inpolygon(x(1),y(1),Ks(j).X,Ks(j).Y);
            if ~isempty(in)
                KOMM(missing(i)) = extractfield(Ks(j),'KOMMUNENUM');
                has(missing(i))  = 1;
            end
        end
    end
end
fprintf('Have: %i of %i Roadlinks "KOMM", missing: %i \n',sum(has),length(RLinks),length(RLinks)-sum(has))
missing=find(has==0);
if ~isempty(missing)
    % If startpoint doesnt work, do it on the ending point
    for i=1:length(missing)
        x=extractfield(RLinks(missing(i)),'X');
        y=extractfield(RLinks(missing(i)),'Y');
        for j=1:length(Ks)
            in = inpolygon(x(end-1),y(end-1),Ks(j).X,Ks(j).Y);
            if ~isempty(in)
                KOMM(missing(i)) = extractfield(Ks(j),'KOMM');
                has(missing(i))  = 1;
            end
        end
    end
end
fprintf('Have: %i of %i Roadlinks "KOMM", missing: %i \n',sum(has),length(RLinks),length(RLinks)-sum(has))

if ~isfield(RLinks,'DISTANCE')
    for i=1:length(RLinks)
        x=RLinks(i).X;
        y=RLinks(i).Y;
        DISTANCE(i)=sqrt((x(1)+x(end-1))^2+(y(1)+y(end-1))^2);
    end
end



% Convert briefly to a table variable to make some calculations
% that take more time in structures.
T=struct2table(RLinks);
T.KOMM=KOMM';
T.MIDPOINT_X=MIDPOINT_X';
T.MIDPOINT_Y=MIDPOINT_Y';
fprintf('Setfields\n  -- KOMM\n  -- MIDPOINT_X\n  -- MIDPOINT_Y ...on all links \n')
if ~isfield(RLinks,'DISTANCE')
    for i=1:length(RLinks)
        x=RLinks(i).X;
        y=RLinks(i).Y;
        DISTANCE(i)=sqrt((x(1)+x(end-1))^2+(y(1)+y(end-1))^2);
    end
    T.DISTANCE=DISTANCE';
    fprintf('Setfields\n  -- DISTANCE ...on all links \n')
end

Sn  = table2struct(T);
end
