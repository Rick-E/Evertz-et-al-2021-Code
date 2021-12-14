function [Peak_alpha, damper] = alpha_peak(X,Y)

% [~, xal] = min(abs(X-8));
% [~, xau] = min(abs(X-20));
% [~,MxI] = max(Y(xal:xau));
% ID = xal+MxI-1;
% Peak_alpha = X(ID);
% damper = 1;




[~, xal] = min(abs(X-8));
[~, xau] = min(abs(X-15));
x = X(xal:xau);
y = Y(xal:xau)';
% y = y/trapz(x,y);
% 
x0 = [8,0,0.01];
xz = [15,15,1000];
sp = [10,1,0.01];
fitfun = fittype( @(a,b,c,x) (1/pi)*(((pi*b)/2)./((x-a).^2+((pi*b)/2).^2))*c );
% fitfun = fittype( @(a,b,c,x) 1./( (x-a).^2 + (b/(2*pi)).^2)*c );
[fitted_curve] = fit(x,y,fitfun,'Lower',x0,'Upper',xz,'StartPoint', sp);

% [fitted_curve] = fit(x,y,fitfun,'StartPoint',x0);

alpha_fit = coeffvalues(fitted_curve);
% disp(alpha_fit);

Peak_alpha = alpha_fit(1,1);
damper = alpha_fit(1,2);
% disp(Peak_alpha);
% figure
% plot(x,y); hold on;
% plot(x,fitted_curve(x)); hold on;


end