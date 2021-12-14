function DIV = jesh(x,y,f)
% Jensen shannon divergence between two spectra

X = x/(trapz(f,x));
Y = y/(trapz(f,y));

Z = (1./2)*(X+Y);
JENX = X.*log(X./Z);
JENY = Y.*log(Y./Z);
DIVX = trapz(f,JENX);
DIVY = trapz(f,JENY);
DIV = 1/2*DIVX + 1/2*DIVY;


end