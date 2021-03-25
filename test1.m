x1 = [1:1:60];
y1 = Med_M_e;
constant = lsqcurvefit(@f,[0;0;0],x1,y1);
pa = constant(3)
p0 = constant(1)
F = constant(2)
yfit =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3));
 
plot(x1,y1,'*');
hold on
plot(x1,yfit,'r');
hold off
