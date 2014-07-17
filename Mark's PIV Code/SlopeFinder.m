% SlopeFinder
%
% use: [out]=getstrain_lsq2(y,u)
%
% Calculates the strain for a given displacement field in the x-direction
% (u) as a function of vertical distance(y) 
%
% y and u are both one-dimensional arrays 
%
% This function is just like getstrain_5pt, except it does a 2 point least
% squares fit on the first point (i.e., a forward difference calculation), a 
% three point least squares fit on the second point, and a three point least 
% squares fit on the (n-1)th point.  The strain at the other n-4 points is
% calculated via a five point least squares fit. 
% 
% Unlike getstrain_lsq, this function can successfully deal with data with
% NaNs
%
function [out] = SlopeFinder(y,u)
strain=zeros(length(y),1);
nanornot=isnan(u);
for n=1:1:length(y)-2
    if (nanornot(n) ~= 1) 
       strain(n)= (u(n+1)-u(n))/(y(n+1)-y(n));
       p=polyfit(y(n:n+2),u(n:n+2),1);
       strain(n+1)=p(1);
       for i=n+2:1:length(y)-2
          p=polyfit(y(i-2:i+2),u(i-2:i+2),1);
          strain(i)=p(1);
       end
    break
    else
       strain(n)=NaN;
    end
end
p=polyfit(y(length(y)-2:length(y)),u(length(y)-2:length(y)),1);
strain(length(y)-1)=p(1);
strain(length(y))= NaN;
out=strain;

