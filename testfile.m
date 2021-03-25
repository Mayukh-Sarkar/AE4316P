%rng default % for reproducibility
tdata = 0:0.1:10;
ydata = 40*exp(-0.5*tdata) + randn(size(tdata));

fun = @(x)sseval(x,tdata,ydata);
x0 = rand(2,1);
bestx = fminsearch(fun,x0)