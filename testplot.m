%% error plot
x = [1:1:60];
figure(1)
scatter(x,Med_M_e,'k','o')
hold on
plot(x,yfit_e,'k','LineWidth',1.2);
hold on
scatter(x,y_nm_e,'k','o','filled')
hold on
plot(x,yfit_nm_e,'color',[0 0 0],'LineStyle','--','LineWidth',1.2);
xlabel('Trial runs')
ylabel('RMS(e), deg')
legend({'RMS(e) M group','Learning curve  M group','RMS(e) NM group','Learning curve NM group'},'Location','northeast','Orientation','vertical')
legend('location','northeast')
%% stick input
 figure(2)
scatter(x,med_u,'k','o')
hold on
plot(x,yfit_u,'k','LineWidth',1.2);
hold on
scatter(x,y_nm_u,'k','o','filled')
hold on
plot(x,yfit_nm_u,'color',[0 0 0],'LineStyle','--','LineWidth',1.2);
xlabel('Trial runs')
ylabel('RMS(u), deg')
legend({'RMS(u) M group','Learning curve  M group','RMS(u) NM group','Learning curve NM group'},'Location','northeast','Orientation','vertical')
legend('location','northwest')

%% stick deflection
figure(3)
scatter(x,med_x,'k','o')
hold on
plot(x,yfit_x,'k','LineWidth',1.2);
hold on
scatter(x,y_NM_x,'k','o','filled')
hold on
plot(x,yfit_nm_x,'color',[0 0 0],'LineStyle','--','LineWidth',1.2);
xlabel('Trial runs')
ylabel('RMS(x), deg')
legend({'RMS(x) M group','Learning curve  M group','RMS(x) NM group','Learning curve NM group'},'Location','northeast','Orientation','vertical')
legend('location','northeast')