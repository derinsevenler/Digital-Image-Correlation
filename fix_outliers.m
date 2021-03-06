function A = fix_outliers(A)
% fix_outliers 

% Identify outliers as values greater than 3 standard deviations

Aelements = A(isfinite(A));
[r c] = size(A);
threshold= 3*std(Aelements);
outliervals = Aelements(abs(Aelements-mean(Aelements))>threshold);
numpoints = length(outliervals);
% outliers are replaced with the average of their 8 nearest neighbors. If
% the outlier is on an edge, it is set to NaN.
for i = 1:numpoints
    [outr, outc] = find(A == outliervals(i),1);
    if (outr>1 && outr<r && outc>1 && outc<c)
        A(outr,outc) = (A(outr-1,outc-1) + A(outr-1,outc) + ...
            A(outr-1,outc+1) + A(outr,outc-1) + A(outr,outc+1) + ...
            A(outr+1,outc-1) + A(outr+1,outc) + A(outr+1,outc+1))/8;
    else
        A(outr,outc) = NaN;
    end
end

