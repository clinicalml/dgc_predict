file = '/Users/rhodos/Desktop/Dropbox/LINCS/data/annot/lincs_annot/pertAnnot_columns/test_ids.txt';
s = ReadStrVec(file);
assert(length(s)==14);
assert(strcmp(s(13), 'BRD-A00520476'));