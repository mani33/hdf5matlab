function x = subsref(br, s)
% Subscripting.
%   x = br(samples, channels). channels can be either channel indices or
%   't' for the timestamps in milliseconds.
%
% MS 2015-04-22


fname = br.fileName;
[~,~,ext] = fileparts(fname);

if ~strcmp(ext,'.ncs')
    x = useHdf5SubsRef(br,s);
    return
end

% make sure subscripting has the right form
assert(numel(s) == 1 && strcmp(s.type, '()') && numel(s.subs) == 2, ...
    'MATLAB:badsubscript', 'Only subscripting of the form (samples, channels) is allowed!')

% samples and channels
samples = s(1).subs{1};
channels = s(1).subs{2};

% Neuralynx data are in blocks of 512 samples. So we first find
% the block range in which the requested data resides

if ischar(channels) && strcmp(channels,'t') % Time stamps requested for given sample indices
    
    assert(br.t0 > 0, 't0 has not been updated in this file!')
    if iscolon(samples)
        x = br.t0 + 1e6 * (0:br.nbSamples-1)' / br.Fs;
    else
        x = br.t0 + 1e6 * (samples(:)-1)' / br.Fs;
    end
else
    
    if ischar(channels) && strcmp(channels,'t_range')
        % When samples are requested corresponding to a time interval
        t_range = samples;
        x = getSamplesForTimeInt(fname,t_range,br);
        
    elseif (ischar(channels) && ~strcmp(channels,'t')) || (isnumeric(channels) || strcmp(channels, ':'))% sample index based data request
        % We assume that the sample indices requested correspond to
        % uninterruped samples i.e, no sample lost during acquisition and if
        % there was any lost samples, we will replace them with NaN's.
        % Since sometimes part or whole of the 512 sample records are lost
        % during acquisition, we will first find the timestamp of the beginning
        % of the requested sample based on t0.
        t_range(1) = br.t0 + 1e6*((samples(1)-1) * (1/br.Fs));
        t_range(2) = br.t0 + 1e6*((samples(end)-1) * (1/br.Fs));
        x = getSamplesForTimeInt(fname,t_range,br);
    end
    
    % scale to volts
    order = numel(br.scale);
    if order == 1
        x = x * br.scale;
    else
        y = 0;
        for i = 1:order
            y = y + x.^(i - 1) * br.scale(i);
        end
        x = y;
    end
    if br.input_inverted
        x = -x;
    end
end

function x = getSamplesForTimeInt(fname,t_range,br)
[timeStamps,nValidSamples,records] =  Nlx2MatCSC(fname,[1 0 0 1 1],0,4,t_range);
[rr, tt] = get_voltage_and_time(timeStamps, nValidSamples, records, br);
% Get the created time stamp indices nearest to the requested
% timestamp indices.
[~, start] = min(abs(tt-t_range(1)));
[~,stop] = min(abs(tt-t_range(2)));

x(:,1) = rr(start:stop);

function [rr, tt] = get_voltage_and_time(timeStamps, nValidSamples, records, br)
% First check the missing part of the data
rs = getRecordSize(br);
badRecords = find(nValidSamples < rs);
nBad = length(badRecords);
tfs = 1e6*(1/br.Fs);
bd = struct;
nRecords = length(timeStamps);
if nBad > 0
    for iRec = 1:nRecords
        nv = nValidSamples(iRec);
        iTs = timeStamps(iRec);
        if nv < rs
            % The code below won't work if the last record had lost
            % some data because of iRec+1. Let's hope that never
            % happens.
            tEnd = timeStamps(iRec) + (nv-1)*tfs;
            nFillSamples = round((((timeStamps(iRec+1))-tEnd)/tfs)) - 1;
            rr = cat(1,records(1:nv,iRec),nan(nFillSamples,1));
            tt = (0:(length(rr)-1))*tfs + iTs;
        else
            rr = records(:,iRec);
            tt = (0:(rs-1))*tfs + iTs;
        end
        bd.rr{iRec} = rr;
        bd.tt{iRec} = tt(:);
    end
    rr = cat(1,bd.rr{:});
    tt = cat(1,bd.tt{:});
else
    f = tfs*repmat((0:(rs-1))',1,nRecords); % Neuralynx time stamps are in microseconds so used 1e6
    tsm = repmat(timeStamps,rs,1);
    tt = f + tsm;
    tt = tt(:);
    rr = records(:);
end

function b = iscolon(x)
b = ischar(x) && isscalar(x) && x == ':';