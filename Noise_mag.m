function [Noise_magnitude, Noise_beta] = Noise_mag(X,Y)

[~, xal] = min(abs(X-15)); % 15 for exp data 
[~, xau] = min(abs(X-40)); % 40 for exp data 
x = X(xal:xau);
y = Y(xal:xau);


x0 = [1,2];
% xz = [1e-1,2.5];
% sp = [1e-4,2];

fitfun = fittype( @(a,b,x) a./x.^b);
% [fitted_curve] = fit(x,y,fitfun,'Lower',x0,'Upper',xz,'StartPoint', sp);
[fitted_curve] = fit(x,y,fitfun,'StartPoint',x0);
noise_fit = coeffvalues(fitted_curve);
% disp(noise_fit);
Noise_magnitude = noise_fit(1,1);
Noise_beta = noise_fit(1,2);
% 
% figure
% loglog(x,y); hold on;
% loglog(x,fitted_curve(x),'r--'); hold on;


% plot(x,y); hold on;
% Noise_magnitude = mean(y);
% disp(Noise_magnitude);
end