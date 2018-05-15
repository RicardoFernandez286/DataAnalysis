function [Ids, DiffVals] = findClosestId2Val(X,Vals)
% Finds the index of the most similar value in an array.
%
% @Input:
%   X (required)
%     A nx1 column vector with the data.
%
%   Vals (required)
%     A 1xk finite row vector with values for which the indices should be found.
%     If 'Vals' contains non-finite values a 'DJM:NonFiniteValues'-warning will
%     be thrown. The according ids will be 1 for -INF, length(X) for INF and NaN
%     for NaN. So check the ids in case of warnings!
%
% @Output:
%   Ids
%     A 1xk vector with the index of the closest values of 'Vals' in 'X'. Ids =
%     1 or Ids = length(X) indicates, that there is no good match to 'Vals' in
%     'X' of some values where infinite.
%
%   DiffVals
%     A 1xk vector with the absolute differences between 'Vals' and the
%     closest value in 'X' (e.g. abs(X(Ids)-Vals)).
%
% @Remarks:
%   - No checking is done on the dimensions of input arguments, so choose your
%     inputs wisely!
%
% @Dependancies
%   none
%
% @Changelog
%   2011-11-13 (DJM): Release.
%   2015-09-29 (DJM): Added support for multiple values.
%   2015-10-28 (DJM): Added some improvements for speed.
%   2016-04-19 (DJM): Added the handling of non-finite values in 'Vals'.
%   2016-08-08 (DJM):
%     - Ids will now be 1 for -INF and length(X) for INF values in 'Vals' as
%       this is closer to the expected behavior. NaNs will still produce NaNs in
%       the ids.
%     - The according warning is been split into a NaN and an INF version. The
%       identifier is still the same for backwards compatibility.

%% Check input
%Assure column vector.
if isrow(X)
	X = X';
end

%Assure row vector.
if iscolumn(Vals)
	Vals = Vals';
end

%Check all-finite.
Ids = NaN(size(Vals));
IsFinite = isfinite(Vals);
if ~all(IsFinite)
	IsInf = isinf(Vals);
	if any(IsInf)
		warning('DJM:NonFiniteValues', ['''Vals'' contained %d infinite ' ...
			'value(s)! The according id(s) will be 1 for -INF or length(X) ' ...
			'for INF!'], nnz(IsInf));
		Ids(Vals == -inf) = 1;
		Ids(Vals ==  inf) = length(X);
	else %NaNs
		warning('DJM:NonFiniteValues', ['''Vals'' contained %d NaN value(s)! ' ...
			'The according id(s) will also be set to NaN!'], nnz(~IsFinite));
	end
end

%% Compute
if any(IsFinite)
	%Get the diffs.
	if isscalar(Vals)
		AbsDiff = X - Vals(IsFinite); %Slightly faster than BSXFUN with scalar.
	else
		AbsDiff = bsxfun(@minus,X,Vals(IsFinite));
	end
	AbsDiff = abs(AbsDiff);

	%Compute the ids.
	[MinDiffVals, Ids(IsFinite)] = min(AbsDiff);
end

%Get the diff values if requested.
if nargout() == 2
	DiffVals = NaN(size(Vals));
	DiffVals(IsFinite) = MinDiffVals;
	if ~all(IsFinite)
		DiffVals(IsInf) = inf; %Since abs we only have positive inf vals.
	end
end