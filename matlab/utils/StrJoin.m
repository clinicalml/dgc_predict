function [str] = StrJoin(cellArray, delimiter)

if ~exist('delimiter')
  delimiter = '';
end

str = '';
for i = 1:length(cellArray)
	str = [str delimiter cellArray{i}];
end

if length(str) > length(delimiter)
  str = str(length(delimiter)+1:end);
end

end
