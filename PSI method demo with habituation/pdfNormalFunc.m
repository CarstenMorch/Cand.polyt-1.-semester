function y = pdfNormalFunc(x, Mean, SD)

%y = exp(-1.*(x-Mean).^2./(2*SD.^2))./((2*pi).^.5*SD);
y = (1/(SD*sqrt(2*pi))).*exp(-1.*((x-Mean).^2)./(2*SD^2)) 
