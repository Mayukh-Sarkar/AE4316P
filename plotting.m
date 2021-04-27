% motion group
x = [1:1:60];
%plotting error
figure(1)
scatter(x,Med_M_e,'k','o')
hold on
plot(x,yfit_e,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_e+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_e-0.03,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS(e)')
legend({'RMS(e)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
legend('location','northeast')
%%plotting stick input
figure(2)
scatter(x,med_u,'k','o')
hold on
plot(x,yfit_u,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_u+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_u-0.03,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS(u)')
legend({'RMS(u)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plotting stick deflection
figure(3)
scatter(x,med_x,'k','o')
hold on
plot(x,yfit_x,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_x+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_x-0.03,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS(x)')
legend({'RMS(u)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%% non motion
%plotting error
figure(4)
scatter(x,y_nm_e,'k','o')
hold on
plot(x,yfit_nm_e,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_e+0.02,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_e-0.02,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS(e)')
legend({'RMS(e)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plotting stick input
figure(5)
scatter(x,y_nm_u,'k','o')
hold on
plot(x,yfit_nm_u,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_u+0.05,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_u-0.05,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS(u)')
legend({'RMS(u)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plotting stick deflection
figure(6)
scatter(x,y_NM_x,'k','o')
hold on
plot(x,yfit_nm_x,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_x+0.01,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_x-0.01,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS(x)')
legend({'RMS(x)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
