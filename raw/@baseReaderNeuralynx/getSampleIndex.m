function idx = getSampleIndex(br, t)
% idx = getSampleIndex(br, t)
% Returns the sample indices 'idx' for the vector of timestamps 't'.
% Neurlynx timestamps are in microseconds. Otherwise, you need to update
% this function
% Mani Subramaniyan 2015-04-17
idx = round(1e-6 * (t - br.t0) * br.Fs) + 1;
idx(idx < 1) = nan;
idx(idx > br.nbSamples) = nan;
