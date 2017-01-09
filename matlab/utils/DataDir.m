function dataDir = DataDir(str)

if ~exist('str')
    str = '';
end

fid = fopen('config/config.txt', 'r');
tline = fgets(fid);
lineParsed = strread(tline,'%s','delimiter','=');

while( ~strcmp(lineParsed(1),'DATAPATH') && ischar(tline))
	tline = fgets(fid);
	lineParsed = strread(tline,'%s','delimiter','=');
	dataDir = lineParsed(2);
end

dataDir = [dataDir{1} '/' str];

fclose(fid);

end
