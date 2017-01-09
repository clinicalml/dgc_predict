d = DataDir();
assert(exist(d)==7);
a = strsplit(d,'/');
if (length(a{end})==0)
    a = a(1:end-1);
end
assert(strcmp(a{end},'data'))
