function [T,pertIds,geneIds,cellIds] = GetTensor(tensorName, printFlag)
if ~exist('tensorName')
    tensorName = 'toy_dense';
end
                
if ~exist('printFlag')
    printFlag = false;
end

if(printFlag) disp(sprintf('loading %s...', tensorName)); end

out = load(tensorName);
T = out.T;
pertIds = out.pertIds;
cellIds = out.cellIds;
geneIds = out.geneIds;

if(printFlag) disp('...done!'); end
    
end
