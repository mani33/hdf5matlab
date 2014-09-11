function rs = getRecordSize(br)
% Get record size of data. Typicall Neuralynx saves the data in 512 sample
% chunks. A record simply refers to a 512 sample block.
rs = br.recordSize;