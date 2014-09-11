function br = baseReaderNeuralynx(fileName,channels)
% Base reader for electrophysiology recordings
%   br = baseReaderElectrophysiology(fileName) opens a base reader for the
%   file given in fileName.
%
%   br = baseReaderElectrophysiology(fileName, channels) opens a reader for
%   only the given channels, where channels is either a numerical vector of
%   channel indices, a string containing a channel name or a cell array of
%   stings containig multiple channel names.
%
%   br = baseReaderElectrophysiology(fileName, pattern) opens a reader for
%   a group of channels matching the given pattern. Channel groups can be
%   for instance tetrodes. In this case the pattern would be 't10c*'.
%
% Example
% br = baseReaderNeuralynx('y:\ephys\2014-09-09_11-10-12\tt*.ncs) - This
% will read all ncs files starting with tt.
% br = baseReaderNeuralynx('y:\ephys\2014-09-09_11-10-12\tt_1_2.ncs) - This
% will read only the channel named tt_1_2.
% br = baseReaderNeuralynx('y:\ephys\2014-09-09_11-10-12\tt_1_*.ncs) - This
% will read all 4 channels of tetrode 1.
% MS 2014-09-09

br.fileName = fileName;

basePath = fileparts(fileName);
% Get all file names specified by the input fileName
d = dir(fileName);
allFileNames = arrayfun(@(x) fullfile(basePath,x.name),d,'uni',false);

oneFile = allFileNames{1};
br.fp = []; % file pointer



% Read just one channel and find total samples
br.nbSamples = getNsamples(oneFile);

% % Find true channel indices. These are physical AD channel indices assigned by
% % Neuralynx hardware.
% br.chIndices = getChannelIndices(fileNames);
allADchanIndices = getChannelIndices(allFileNames);

allChanNames = getChanNames(allFileNames);
selFileNames = cell(1,1);
if nargin < 2
    br.chIndices = allADchanIndices;
    selFileNames = allFileNames;
    br.chNames = allChanNames;
else
    
    % If channels are numbers, then those numbers must be the AD Channel
    % numbers assigned by Neuralynx hardware
    if isnumeric(channels)
        % Find all available AD Channels
        br.chIndices = channels;
        [~,fileIndices] = intersect(allADchanIndices,channels);
        selFileNames = allFileNames(fileIndices);
        br.chNames = allChanNames(fileIndices);
    else % pattern such as t10c1 or t10c* or t*c*
        % one channel string or cell array of strings?
        if ischar(channels)
            channels = {channels};
        end
        nChanGrp = length(channels);
        for i = 1:nChanGrp
            ch = strrep(channels{i},'*','\d+'); % to make it easy for regexp
            % Find which file name(s) have these channels           
            fileInd = cellfun(@(x) ~isempty(regexp(x,ch,'start')),allFileNames);
            selFileNames{i} = allFileNames(fileInd);
            br.chIndices{i} = allADchanIndices(fileInd);
            br.chNames{i} = allChanNames(fileInd);
        end
        br.chIndices = [br.chIndices{:}];
        br.chNames = [br.chNames{:}];
        selFileNames = [selFileNames{:}];
    end
end
br.fileList = selFileNames;
% Find number of channels
br.nbChannels = length(br.chNames);


[br.t0,br.Fs,header] = Nlx2MatCSC(oneFile,[1 0 1 0 0],1,3,1);
scaleInd = cellfun(@(x) ~isempty(strfind(x,'ADBitVolts')), header);
hh = header{scaleInd};
sc = regexp(hh,'ADBitVolts.+(\d+\.\d+\w\-\d+)','tokens');
br.scale = str2double(sc{:});
br.recordSize = size(Nlx2MatCSC(oneFile,[0 0 0 0 1],0,3,1),1);
br = class(br, 'baseReaderNeuralynx');

function nSamples = getNsamples(fileName)
% You must verify that this function returns the correct number of total
% samples by reading a whole file using Nlx2MatCSC function. Here, we will
% use some hard coded information such as headersize to determine total
% samples.
header = Nlx2MatCSC(fileName,[0 0 0 0 0],1,3,1);
hInd = cellfun(@(x) ~isempty(regexp(x,'-RecordSize\s\d+','match')),header);
hh = header{hInd};
rs = regexp(hh,'-RecordSize\s(\d+)','tokens');
recordsize = str2double(rs{:});
% the file starts with a 16*1024 bytes header in ascii, followed by a number of records
fid = fopen(fileName, 'rb', 'ieee-le');

% determine the length of the file
fseek(fid, 0, 'eof'); % move pointer to the end of file
headersize = 16384; % DANGER! - you got to make sure that Neuralynx doesn't change these numbers when they
% change their acquisition software.
NRecords   = floor((ftell(fid) - headersize)/recordsize);
dataRecordSize = size(Nlx2MatCSC(fileName,[0 0 0 0 1],0,3,1),1);
nSamples = NRecords*dataRecordSize;
fclose(fid);


function chNames = getChanNames(fileNames)
nFiles = length(fileNames);
chNames = cell(1,nFiles);
for i = 1:nFiles
    [~,chNames{i}] = fileparts(fileNames{i});
end

function chIndices = getChannelIndices(fileNames)

% Get entire header info first
nFiles = length(fileNames);
chIndices = nan(1,nFiles);
for i = 1:nFiles
    fn = fileNames{i};
    h = Nlx2MatCSC(fn,[0 0 0 0 0],1,3,1);
    adInd = cellfun(@(x) ~isempty(regexp(x,'-ADChannel\s','start')), h);
    hh = h{adInd};
    ci = regexp(hh,'-ADChannel\s+(\d+)\s+','tokens');
    ci = ci{:};
    chIndices(i) = str2double(ci);
end