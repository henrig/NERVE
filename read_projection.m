function prj = read_projection(shapefileName)
% Summary of this function goes here
%   Detailed explanation goes here

if exist(sprintf('%s%s',shapefileName,'.prj'))
    fid  = fopen(sprintf('%s%s',shapefileName,'.prj'),'r');
    txt = textscan(fid,'%s');
    fclose(fid);
    
elseif strcmp(shapefileName(length(shapefileName)-3:length(shapefileName)),'.prj')
    fprintf('%s\n',shapefileName)
    fid  = fopen(shapefileName,'r');
    txt = textscan(fid,'%s');
    fclose(fid);
    
elseif exist(sprintf('%s%s',shapefileName(1:end-4),'.prj'))
    fprintf(sprintf('%s%s\n',shapefileName(1:end-4),'.prj'))
    fid  = fopen(sprintf('%s%s',shapefileName(1:end-4),'.prj'),'r');
    txt = textscan(fid,'%s');
    fclose(fid);
    
else
    fprintf('%s\n',shapefileName)
    error('file not found')
end

a=char(txt{1,1});
[l,r]=size(a);
prj=[];
for i=1:l
    prj=strcat(prj,a(i,:));
end

end

