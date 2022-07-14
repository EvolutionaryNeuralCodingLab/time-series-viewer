function []=saveLayoutFiles(target,ext,chLayout)
extLen=numel(ext);
if iscell(target)
    d=target;
elseif isdir(target)
    d=dir([target filesep '*.' ext]);
    d={d.name};
end

for i=1:numel(d)
    fid=fopen([d{i}(1:(end-extLen-1)) '.chMap'],'w');
    fprintf(fid,chLayout);
    fclose(fid);
end
    
