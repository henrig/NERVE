function Save_shape(Sout,ofile,prj_file)
% Small function to save a Structure (or a suited table) as an ESRI
% shapefile with shapewrite. MatLab needs to copy a .prj file as it does
% not create these files by default.

% INPUT: Sout     --- Shape or Table variable for writing.
%        ofile    --- Filename of written file (no extension)
%        prj_file --- Name of file (without extension) of projection to be
%                      copied.
% -------------------------------------------------------------------------
% a = class(Sout);
if istable(Sout)
    Sout = table2struct(Sout);
end
dbf = makedbfspec(Sout);
shapewrite(Sout, ofile, 'DbfSpec', dbf)
pfilename = strcat(ofile,'.prj');
fprintf('Wrote Shape file:\n %s\n ',ofile)
try
    prj    = HEDGE_read_projection(prj_file);
    fid       = fopen(pfilename,'w+');
    fprintf(fid,'%s',prj);
    fclose(fid);
    fprintf('Wrote Projection file:\n %s\n ',ofile)
catch
    fprintf('No Valid Projection given\n ')
end
end