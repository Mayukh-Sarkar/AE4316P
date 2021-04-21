
res = [-1,0,1,2,3,4] ;
y = zeros(length(res),1);
for i = 1 : length(res)
    y(i) = 4 - 2*res(i) + 1.5 * res(i)^2 - 0.75 * res(i)^3 ;
    
end

f = @(x) 0 ;
% options = optimset('Display','iter','PlotFcns',@optimplotfval);

for i = 1 : length(res)
    g = @(x) (y(i) - (x(1) - x(2) * res(i) + x(3)* res(i)^2 + x(4) * res(i)^3 ))^2 ;
    f = @(x) f(x) + g(x) ;
end

x0 = [3, 1, 3] ;
% [x,fval,exitflag,output] = fminsearch(f, x0,options);
x = fminsearch(f,x0) ;
