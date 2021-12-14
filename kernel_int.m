function K = kernel_int(x,y,F_A)
% Fredholm integral equation kernel - approximate lorentzian form

% Passed arguments [ frequency vector (x), damping vector (y) and peak alpha frequency (F_A) ]
% Function returns power (K) for each value of frequency (x) and damping (y)

K = (1./(4*pi).^2)*(1./( (y./(2*pi)).^2 + (x - F_A).^2));

end
