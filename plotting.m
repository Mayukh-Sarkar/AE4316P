% motion group
x = [1:1:60];
%plottinf error
figure(1)
scatter(x,Med_M_e,'k','o')
hold on
plot(x,yfit_e,'r','LineStyle','--','LineWidth',1.2);
%hold on
%plot(x,yfit_e+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
%hold on
%plot(x,yfit_e-0.03,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS_m(e)')
legend({'RMS(e)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
legend('location','northeast')
%%%plottinf stick input
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
ylabel('RMS_m(u)')
legend({'RMS(u)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%%plottinf stick deflection
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
ylabel('RMS_m(x)')
legend({'RMS(x)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')

figure(4)
%plottinf pilot gain
scatter(x,kp_m,'k','o')
hold on
plot(x,yfit_m_kp,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_m_kp+0.1,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_m_kp-0.1,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('K_p_m')
legend({'K_p','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
legend('location','northeast')
%%plottinf time lead
figure(5)
scatter(x,Tl_m,'k','o')
hold on
plot(x,yfit_m_Tl,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_m_Tl+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_m_Tl-0.03,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\tau_L_m')
legend({'\tau_L','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%%plottinf pilot time
figure(6)
scatter(x,Tp_m,'k','o')
hold on
plot(x,yfit_m_Tp,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_m_Tp+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_m_Tp-0.03,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\tau_p_m')
legend({'\tau_p','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%%plottinf damping
x_new = [2:1:60];
figure(7)
scatter(x_new,zeta_m(2:end,:),'k','o')
hold on
plot(x_new,yfit_m_zeta,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x_new,yfit_m_zeta+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x_new,(yfit_m_zeta-0.03),'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\zeta_{m}')
legend({'\zeta_{m}','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%%plottinf omemga
figure(8)
scatter(x_new,omega_m(2:end,:),'k','o')
hold on
plot(x_new,yfit_m_omega,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x_new,yfit_m_omega+0.1,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x_new,yfit_m_omega-0.1,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\omega_{m}')
legend({'\omega_{m}','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')

%% non motion
%plottinf error
figure(10)
scatter(x,y_nm_e,'k','o')
hold on
plot(x,yfit_nm_e,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_e+0.02,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_e-0.02,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS_{nm}(e)')
legend({'RMS(e)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%%plottinf stick input
figure(11)
scatter(x,y_nm_u,'k','o')
hold on
plot(x,yfit_nm_u,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_u+0.05,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_u-0.05,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS_{nm}(u)')
legend({'RMS(u)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%%plottinf stick deflection
figure(12)
scatter(x,y_NM_x,'k','o')
hold on
plot(x,yfit_nm_x,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_x+0.01,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_x-0.01,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('RMS_{nm}(x)')
legend({'RMS(x)','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')

%plottinf pilot gain
figure(13)
scatter(x,kp_nm,'k','o')
hold on
plot(x,yfit_nm_kp,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_kp+0.1,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_kp-0.1,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('K_p_{nm}')
legend({'K_p','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
legend('location','northeast')
%plottinf time lead
figure(14)
scatter(x,Tl_nm,'k','o')
hold on
plot(x,yfit_nm_Tl,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_Tl+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_Tl-0.03,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\tau_L_{nm}')
legend({'\tau_L','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plottinf pilot time
figure(15)
scatter(x,Tp_nm,'k','o')
hold on
plot(x,yfit_nm_Tp,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_Tp+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x,yfit_nm_Tp-0.03,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\tau_p_{nm}')
legend({'\tau_p','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plottinf damping
x_new = [2:1:60];
figure(16)
scatter(x_new,zeta_nm(2:end,:),'k','o')
hold on
plot(x_new,yfit_nm_zeta,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x_new,yfit_nm_zeta+0.03,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x_new,(yfit_nm_zeta-0.03),'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\zeta_{nm}')
legend({'\zeta_{nm}','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plottinf omemga
figure(17)
scatter(x_new,omega_nm(2:end,:),'k','o')
hold on
plot(x_new,yfit_nm_omega,'r','LineStyle','--','LineWidth',1.2);
hold on
plot(x_new,yfit_nm_omega+0.1,'color',[0 0.5 0],'LineStyle','--','LineWidth',1.2);
hold on
plot(x_new,yfit_nm_omega-0.1,'color',[0 0 0.5],'LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\omega_{nm}')
legend({'\omega_{nm}','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
