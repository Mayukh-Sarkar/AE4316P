x1 = [1:1:60];
y1 = Med_M_e;
y = Med_M_e(end)+((1-y1).^x1)*(Med_M_e(1)-Med_M_e(end));
function sse = sseval(x2,x1,y)

y1 = x(2);
sse = sum((y1 - y).^2);
end 
fun = @(x)sseval(x,x1,y);
 
x0 = [Med_M_e(1);Med_M_e(end)];
bestx = fminsearch(fun,x0);
F = bestx(2);
yfit = Med_M_e(end)+((1-0.0475).^x1)*(Med_M_e(1)-Med_M_e(end));
 
plot(x1,y1,'*');
hold on
plot(x1,yfit,'r');
hold off
