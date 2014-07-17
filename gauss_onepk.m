function F = gauss_onepk(x,xdata)
        F = (x(1)*exp((-(xdata-x(2)).^2)./(2.*x(3).^2)))  + x(4);