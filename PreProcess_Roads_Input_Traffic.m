function [RLinks] = PreProcess_Roads_Input_Traffic()
%--------------------------------------------------------------------------
% Module to read in an original Road Network and traffic file.
% Processes data to calculate the required fields and splits it to
% municipal level data.
% Function is not generic and require input of specified format (.shp) and
% also require existing fieldnames specified.
%
% -INPUT  : this function reads road and traffic network from file (ESRI shape)
% -OUTPUT : this function adds the road link shapes per municipality to the
%           temporary folder (tpath).
%
% Functions also read ancillary files to define fields.
%   read_projection()
%   Roads_Add_Urban()
%   Roads_Add_HBEFA_Parameters()
%   Roads_Add_Advanced_Properties() 
%   Roads_Add_Width()
%   Roads_Scale_Traffic_to_Year()
%   Roads_Congestion_Parameters()
%
% 20.09.2020 -Henrik Grythe
% Kjeller NILU
%--------------------------------------------------------------------------
global tfold use_temporary_files
global traffile Komm_shape Listfields

fprintf('in PreProcess_Input_Traffic.\n')

% compare the list of field names to a list of requred fields
required_fields = [{'Geometry'},{'X'},{'Y'},{'KOLL_ADT'},{'LETTE_BILE' },{'GODS_ADT'},{'SUM_ADT'}];

% list of fields needed to do (if needed) matching on municipality grids:
location_fields = [{'KOMM'},{'MIDPOINT_X'},{'MIDPOINT_Y'},{'DISTANCE'}];

% list of fields erquired for advanced calculations:
advanced_fields = [{'FM_TIME'},{'KO_MORGEN'},{'KO_ETTERM'},{'FM_SPEED'},{'SHAPE_LENG'},{'VK'},{'STIGNING_P'},{'SPEED'}];

% list of fields generated by the model
fields_made = [{'URBAN'},{'RUSH_DELAY'},{'RUSH_SPEED'},{'SLOPE'},{'HBEFA_NUM'},{'HBEFA_EQIV'}];

% List of fields a finished file has:
Listfields = [{'Geometry'},{'X'},{'Y'},{'BoundingBox'},{'DISTANCE'},{'KOMM'},{'SPEED'},{'SLOPE'},...
    {'URBAN'},{'RUSH_DELAY'},{'HBEFA_EQIV'},{'WIDTH'},{'CAPACITY'},{'N_LANES'},{'L_ADT'},{'H_ADT'},{'B_ADT'}];

if use_temporary_files
    try
        % First do tests if there is a temporary file saved:
        fprintf('Trying to read temp file\n%s ...\n',tfiles.RL)
        RLinks = shaperead(tfiles.RL);
        fprintf('Found Pre-Processed file:\n %s \n',tfiles.RL)
        fields = fieldnames(RLinks);
        fprintf('Fields:\n')
        for i=1:length(fields)
            fprintf('%s\n',char(fields(i)))
        end
    catch
        fprintf('### No temporary file found ###\n')
        prj    = read_projection(traffile);
        fprintf('Reading large file: \n\t %s  ...',traffile)
        RLinks = shaperead(traffile);
        fprintf('done\n')
        fields = fieldnames(RLinks);
        fprintf('Fields:\n')
        for i=1:length(fields)
            fprintf('%s,',char(fields(i)))
            if rem(i,10)==0;fprintf('\n');end
        end
        fprintf('\n')
    end
else
    prj    = read_projection(traffile);
    fprintf('Reading large file: \n\t %s  ...',traffile)
    RLinks = shaperead(traffile);
    fprintf('done\n')
    fields = fieldnames(RLinks);
    fprintf('Fields:\n')
    for i=1:length(fields)
        fprintf('%s,',char(fields(i)))
        if rem(i,10)==0;fprintf('\n');end
    end
    fprintf('\n')
end


% Test if required fields are there:
for i=1:length(required_fields)
    if isfield(RLinks,required_fields(i))
    else
        required_fields(i)
        error('Missing required field (possibly others)')
    end
end

% Test if required fields are there:
for i=1:length(Listfields)
    if isfield(RLinks,Listfields(i))
    else
        fprintf('Output field not defined yet: %s \n',char(Listfields(i)))
    end
end

%--------------------------------------------------------------------------
miss_location = 0;
for i=1:length(location_fields)
    if isfield(RLinks,location_fields(i))
    else
        fprintf('%s\n',char(location_fields(i)))
        miss_location = 1;
    end
end

if miss_location
    RLinks = Roads_Add_Municipality(RLinks,Komm_shape);
else
    fprintf('Already have necessary fields Location fields \n')
end

%--------------------------------------------------------------------------
RLinks = Roads_Remove_NoTrafficRoads(RLinks);

RLinks = Fix_Roadfields(RLinks);

%--------------------------------------------------------------------------
if isfield(RLinks,'KAPTID_06_') && ~isfield(RLinks,'KOTID_MORGEN')
    RLinks = Roads_Congestion_Parameters(RLinks);
end

miss_advanced = 0;
for i=1:length(advanced_fields)
    if isfield(RLinks,advanced_fields(i))
    else
        miss_advanced = 1;
        fprintf('Missing field to do Advanced calculations: %s \n',char(advanced_fields(i)))
        fprintf('Maybe field %s is misnamed?\n',char(advanced_fields(i)))
        warning('or proceed with simplistic calculation?\n')
    end
end


%--------------------------------------------------------------------------
% ADD fields (if missing)
% URBAN
if ~isfield(RLinks,'URBAN')
    RLinks = Roads_Add_Urban(RLinks);
end

%--------------------------------------------------------------------------
% ADD fields (if missing)
% HBEFA
if ~isfield(RLinks,'HBEFA_NUM') || ~isfield(RLinks,'HBEFA_EQIV')
    RLinks = Roads_Add_HBEFA_Parameters(RLinks);
end
%--------------------------------------------------------------------------
% ADD fields (if missing)
% Advanced
if ~isfield(RLinks,'RUSH_DELAY') || ~isfield(RLinks,'RUSH_SPEED')|| ~isfield(RLinks,'SLOPE')
    RLinks = Roads_Add_Advanced_Properties(RLinks);
end

%--------------------------------------------------------------------------
% ADD fields (if missing)
if ~isfield(RLinks,'WIDTH') && isfield(RLinks,'DEKKEBREDD')|| isfield(RLinks,'LANES')
        RLinks = Roads_Add_Width(RLinks);
end

% Scale Traffic from road year to traffic year:
RnLinks = Roads_Scale_Traffic_to_Year(RLinks);

if use_temporary_files
    Save_shape(RLinks,tfiles.RL)
end

end
%
%     try
%         prj    = read_projection(traffile);
%         fprintf('Reading large file: \n\t %s  ...',fname)
%         RLinks = shaperead(traffile);
%         fprintf('done\n')
%         fields = fieldnames(RLinks);
%         fprintf('Fields:\n')
%         for i=1:length(fields)
%             fprintf('%s,',char(fields(i)))
%             if rem(i,10)==0;fprintf('\n');end
%         end
%         fprintf('\n')
%     catch
%         RLinks=[];
%         warning('Did not find road links at defined path')
%         fprintf('continuing...\n')
%         return
%     end
% 
% % 
% 
% a = strfind(traffile,'/');
% if isempty(a)
%     fname = [];
%     ofile = strcat(tfold,tfiles.RL);
% else
%     fname = traffile(a(end)+1:end);
%     ofile = strcat(tfold,tfiles.RL);
% end
% % Try to load a temporary file:
% try
%     fprintf('Trying to read temp file\n%s_temp...\n',tfiles.RL)
%     RLinks = shaperead(tfiles.RL);
%     fprintf('Found Pre-Processed file:\n %s \n',tfiles.RL)
%     fields = fieldnames(RLinks);
%     fprintf('Fields:\n')
%     for i=1:length(fields)
%         fprintf('%s\n',char(fields(i)))        
%     end
%     return
% catch
%     fprintf('### No temporary file found ###\n')
%     % try to read the road shape
%     try
%         prj    = read_projection(traffile);
%         fprintf('Reading large file: \n\t %s  ...',fname)
%         RLinks = shaperead(traffile);
%         fprintf('done\n')
%         fields = fieldnames(RLinks);
%         fprintf('Fields:\n')
%         for i=1:length(fields)
%             fprintf('%s,',char(fields(i)))
%             if rem(i,10)==0;fprintf('\n');end
%         end
%         fprintf('\n')
%     catch
%         RLinks=[];
%         warning('Did not find road links at defined path')
%         fprintf('continuing...\n')
%         return
%     end
% end