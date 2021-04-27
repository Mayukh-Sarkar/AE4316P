%% motion group
x = [1:1:60];
%plotting error
figure(1)
scatter(x,Med_M_e,'k','o')
hold on
plot(x,yfit_e,'r','LineStyle','--');
hold on
plot(x,yfit_e+0.03,'color',[0 0.5 0],'LineStyle','--');
hold on
plot(x,yfit_e-0.03);
hold off
xlabel('Trial runs')
ylabel('RMS(e)')
legend('RMS(e)','Learning curve')
%plotting stick input
figure(2)
scatter(x,med_u,'*')
hold on
plot(x,yfit_u,'r');
hold off
xlabel('Trial runs')
ylabel('RMS(u)')
legend('RMS(u)','Learning curve')
%plotting stick deflection
figure(3)
scatter(x,med_x,'*')
hold on
plot(x,yfit_x,'r');
hold off
xlabel('Trial runs')
ylabel('RMS(x)')
legend('RMS(x)','Learning curve')
%% non motion
%plotting error
figure(4)
scatter(x,y_nm_e,'k','o')
hold on
plot(x,yfit_nm_e,'r','LineStyle','--');
hold on
plot(x,yfit_nm_e+0.03,'color',[0 0.5 0],'LineStyle','--');
hold on
plot(x,yfit_nm_e-0.03);
hold off
xlabel('Trial runs')
ylabel('RMS(e)')
legend('RMS(e)','Learning curve')
%plotting stick input
figure(5)
scatter(x,y_nm_u,'*')
hold on
plot(x,yfit_nm_u,'r');
hold off
xlabel('Trial runs')
ylabel('RMS(u)')
legend('RMS(u)','Learning curve')
%plotting stick deflection
figure(6)
scatter(x,y_NM_x,'*')
hold on
plot(x,yfit_nm_x,'r');
hold off
xlabel('Trial runs')
ylabel('RMS(x)')
legend('RMS(x)','Learning curve')
