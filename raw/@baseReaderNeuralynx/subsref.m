function x = subsref(br, s)
% Subscripting.
%   x = br(samples, channels). channels can be either channel indices or
%   't' for the timestamps in milliseconds.
%
% AE 2011-04-11

% make sure subscripting has the right form
assert(numel(s) == 1 && strcmp(s.type, '()') && numel(s.subs) == 2, ...
    'MATLAB:badsubscript', 'Only subscripting of the form (samples, channels) is allowed!')

% samples and channels
samples = s(1).subs{1};
channels = s(1).subs{2};

% Neuralynx data are in blocks of 512 samples. So we first find
% the block range in which the requested data resides
rs = getRecordSize(br);

% all samples requested?
if iscolon(samples)
    nSamples = br.nbSamples;
else
    % Check for valid range of samples
    assert(all(samples <= br.nbSamples & samples > 0), 'MATLAB:badsubscript', ...
        'Sample index out of range [1 %d]', br.nbSamples);
    nSamples = numel(samples);
end

% time channel requested?
if ischar(channels) && channels == 't'
    assert(br.t0 > 0, 't0 has not been updated in this file!')
    if iscolon(samples)
        x = br.t0 + 1000 * (0:br.nbSamples-1)' / br.Fs;
    else
        x = br.t0 + 1000 * (samples(:)-1)' / br.Fs;
    end
else
    
    % all channels requested?
    if iscolon(channels)
        channels = 1:(br.nbChannels);
    else
        % Check for valid range of channels
        assert(all(channels <= br.nbChannels & channels > 0), ...
            'MATLAB:badsubscript', 'Channel index out of range [1 %d]', br.nbChannels);
    end
    nChannels = numel(channels);
    
    % Convert to actual channel numbers in the recording file
    %     channels = br.chIndices(channels);
    
    x = zeros(nSamples, nChannels);
    
    if iscolon(samples)
        % reading all samples
        for i = 1:nChannels
            fname = br.fileList{channels(i)};
            records =  Nlx2MatCSC(fname,[0 0 0 0 1],0,1,[]);
            x(:,i) = records(:);
        end
    elseif length(samples) > 2 && samples(end) - samples(1) == length(samples) - 1 && all(diff(samples) == 1)
        % reading continuous block of samples. Add an extra block at the
        % end of the block range so that Nlx2MatCSC does not complain when
        % you need to get samples that reside within a single block.
        block_range = [ceil(samples(1)/rs) ceil(samples(end)/rs)+1];
        for i = 1:nChannels
            fname = br.fileList{channels(i)};
            records =  Nlx2MatCSC(fname,[0 0 0 0 1],0,2,block_range);
            rr = records(:);
            % Get starting and ending points within the blocks of samples
            % retrieved
            start = samples(1)- (block_range(1)-1)*rs;
            stop = start + nSamples - 1;
            x(:,i) = rr(start:stop);
        end
    else
        % reading arbitrary set of samples
        error('reading aribitray samples is not implemented yet')
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
end



function b = iscolon(x)
b = ischar(x) && isscalar(x) && x == ':';
