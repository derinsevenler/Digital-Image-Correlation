function F = gauss_twopk(x,xdata)
        F = (x(1)*exp((-(xdata-x(2)).^2)./(2.*x(3).^2))) + (x(4)*exp((-(xdata-x(5)).^2)./(2.*x(6).^2))) + x(7);