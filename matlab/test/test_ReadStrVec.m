file = DataDir('metadata/test_ids.txt');
s = ReadStrVec(file);
assert(length(s)==14);
assert(strcmp(s(13), 'BRD-A00520476'));
