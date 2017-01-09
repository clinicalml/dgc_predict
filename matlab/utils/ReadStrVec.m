function s = ReadStrVec(file, str)

if ~exist('str')
    str = '%s';
end

fid = fopen(file, 'r');
C = textscan(fid, str);
s = C{1};

end