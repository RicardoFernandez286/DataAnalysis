% anyEq - Fast check if two arrays have any common element
% R = anyEq(X, Y)
% INPUT:
%   X, Y: Arrays of any size. Complex or sparse array are rejected.
%         Types: DOUBLE, SINGLE, (U)INT8/16/32/64, CHAR, LOGICAL.
% OUTPUT:
%   R:    TRUE is replied if any element of X occurs in Y, FALSE otherwise.
%         NaNs are treated as equal.
%
% NOTES:
% - This is equivalent to:
%     R = any(X(:) == Y(1)) || any(X(:) == Y(2)) || ...
% - This MEX version is faster than the Matlab method, because the creation of
%   the intermediate logical array is avoided:
%     Worst case (no matching element):  25% to 60% faster
%     Best case (first element matches): 99.99% faster for 1e6 elements
%     (Matlab 2011b/64, MSVC 2008).
% - For small LOGICAL arrays (< 5'000 elements) anyEq is up to 65% slower than
%   any(X) or ~all(X). For larger arrays the speedup depends on the position of
%   the first match and can reach a factor of 7.
% - X and Y are ordered automatically such, that the elements of the smaller
%   array are searched in the larger one.
%
% EXAMPLES:
%   anyEq(0:0.1:1, 0.3)   % FALSE: Effect of limited precision
%   anyEq(1:4, [5,2])     % TRUE:  2 is found in 1:4
%
% COMPILATION: See anyEq.c for instructions how to compile the C-file.
%
% TEST: Run uTest_anyEq to check validity and speed of the Mex function.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
%         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008/2010
% Assumed Compatibility: higher Matlab versions, Mac, Linux
% Author: Jan Simon, Heidelberg, (C) 2006-2013 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also ANY, ALL, anyExceed.

% $JRev: R-h V:007 Sum:0eKiyncQHl49 Date:11-Sep-2013 01:15:04 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_anyEq $
% $File: Tools\GLSets\anyEq.m $

%#mex
% MEX file! This is a dummy to feed HELP.
