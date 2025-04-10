function [h, crit_p, adj_p, ind_p_sig] = fdr_bh(pvals, q)
% FDR_BH performs false discovery rate correction using the Benjamini-Hochberg procedure.
%
%   [h, crit_p, adj_p] = fdr_bh(pvals, q)
%
%   Inputs:
%       pvals - a vector of p-values from multiple hypothesis tests.
%       q     - the desired false discovery rate (e.g., 0.05). Default is 0.05.
%
%   Outputs:
%       h      - a logical vector (same size as pvals) indicating which hypotheses
%                are rejected (true = reject).
%       crit_p - the critical p-value threshold; p-values less than or equal to this
%                are considered significant.
%       adj_p  - the adjusted p-values based on the Benjamini-Hochberg procedure.
%
%   This implementation follows the approach described in:
%   Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate:
%   a practical and powerful approach to multiple testing. Journal of the Royal
%   Statistical Society, Series B, 57(1), 289–300.

if nargin < 2
    q = 0.05;
end

% Ensure pvals is a column vector
pvals = pvals(:);
m = length(pvals);

% Sort the p-values and get the sorting order
[sorted_p, sort_index] = sort(pvals);
[~, unsort_index] = sort(sort_index);

% Compute thresholds for each sorted p-value
thresholds = ( (1:m)' / m ) * q;

% Find the largest index where the p-value is less than or equal to its threshold
w = find(sorted_p <= thresholds);
if isempty(w)
    crit_p = 0;
    h = false(m,1);
else
    max_k = max(w);
    crit_p = sorted_p(max_k);
    h = pvals <= crit_p;
    ind_p_sig = find(pvals==crit_p);
end

% Compute adjusted p-values
adj_p = zeros(m,1);
min_adj_p = 1;
for i = m:-1:1
    min_adj_p = min(min_adj_p, m * sorted_p(i) / i);
    adj_p(i) = min_adj_p;
end
adj_p = adj_p(unsort_index);

end
