fun = @(F)0.5831 +(1-F).^x2*(0.81883-0.5841);
F = fminsearch(fun,0);