%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% AE4316P Advanced Aerospace Human-Machine Interaction %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Group 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('ae4316P_2021_data_group1.mat')

%%%%%%% Plot RMS(e) %%%%%%%
%% error motion
RMS_M1_e = rms(groupM_data_subj1.e);

RMS_M2_e = rms(groupM_data_subj2.e);

RMS_M3_e = rms(groupM_data_subj3.e);

RMS_M4_e = rms(groupM_data_subj4.e);
RMS_M5_e = rms(groupM_data_subj5.e);
RMS_M6_e = rms(groupM_data_subj6.e);
RMS_M7_e = rms(groupM_data_subj7.e);
RMS_M8_e = rms(groupM_data_subj8.e);

x = [1:1:60];
RMS_M_e = [RMS_M1_e;RMS_M2_e;RMS_M3_e;RMS_M4_e;RMS_M5_e;RMS_M6_e;RMS_M7_e;RMS_M8_e];
Med_M_e = median(RMS_M_e);
figure(1)
scatter(x,Med_M_e,'*')
hold on
x1 = [1:1:60];
y1 = Med_M_e;
constant = lsqcurvefit(@f,[0;0;0],x1,y1);
pae = constant(3)
p0e = constant(1)
Fe = constant(2)
yfite =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3));

hold on
plot(x,yfit,'r');
hold off


%% input motion
%%%%%%% Plot RMS(u) %%%%%%%
RMS_M1_u = rms(groupM_data_subj1.u);
RMS_M2_u = rms(groupM_data_subj2.u);
RMS_M3_u = rms(groupM_data_subj3.u);
RMS_M4_u = rms(groupM_data_subj4.u);
RMS_M5_u = rms(groupM_data_subj5.u);
RMS_M6_u = rms(groupM_data_subj6.u);
RMS_M7_u = rms(groupM_data_subj7.u);
RMS_M8_u = rms(groupM_data_subj8.u);
RMS_M_u = [RMS_M1_u;RMS_M2_u;RMS_M3_u;RMS_M4_u;RMS_M5_u;RMS_M6_u;RMS_M7_u;RMS_M8_u];
x1 = [1:1:60];
y2 = median(RMS_M_u);
constant = lsqcurvefit(@f,[0;0;0],x1,y2);
pau = constant(3)
p0u = constant(1)
Fu = constant(2)
yfitu =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3));

figure(2);
boxplot(RMS_M_u)
hold on
plot(x1,yfitu,'r');
hold off
 %% x motion
%%%%%%% Plot RMS(x) %%%%%%%
RMS_M1_x = rms(groupM_data_subj1.x);
RMS_M2_x = rms(groupM_data_subj2.x);
RMS_M3_x = rms(groupM_data_subj3.x);
RMS_M4_x = rms(groupM_data_subj4.x);
RMS_M5_x = rms(groupM_data_subj5.x);
RMS_M6_x = rms(groupM_data_subj6.x);
RMS_M7_x = rms(groupM_data_subj7.x);
RMS_M8_x = rms(groupM_data_subj8.x);
RMS_M_x = [RMS_M1_x;RMS_M2_x;RMS_M3_x;RMS_M4_x;RMS_M5_x;RMS_M6_x;RMS_M7_x;RMS_M8_x];
x1 = [1:1:60];
y3 = median(RMS_M_x);
constant = lsqcurvefit(@f,[0;0;0],x1,y3);
pax = constant(3)
p0x = constant(1)
Fux = constant(2)
yfitx =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3));
 
figure(3);
boxplot(RMS_M_x)
hold on
plot(x1,yfitx,'r');
hold off

%% nonmotion error
%%%%%%% Plot RMS(e) %%%%%%%
RMS_NM1_e = rms(groupNM_data_subj1.e);
RMS_NM2_e = rms(groupNM_data_subj2.e);
RMS_NM3_e = rms(groupNM_data_subj3.e);
RMS_NM4_e = rms(groupNM_data_subj4.e);
RMS_NM5_e = rms(groupNM_data_subj5.e);
RMS_NM6_e = rms(groupNM_data_subj6.e);
RMS_NM7_e = rms(groupNM_data_subj7.e);
RMS_NM8_e = rms(groupNM_data_subj8.e);

x = [1:1:60];
RMS_NM_e = [RMS_NM1_e;RMS_NM2_e;RMS_NM3_e;RMS_NM4_e;RMS_NM5_e;RMS_NM6_e;RMS_NM7_e;RMS_NM8_e];
ynme = median(RMS_NM_e);
constant = lsqcurvefit(@f,[0;0;0],x1,ynme);
yfitnme =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3));
panme = constant(3)
p0nme = constant(1)
Funme = constant(2)
figure(4);
boxplot(RMS_NM_e)
hold on
plot(x1,yfitnme,'r');
hold off
%% nm inpuut
%%%%%%% Plot RMS(u) %%%%%%%
RMS_NM1_u = rms(groupNM_data_subj1.u);
RMS_NM2_u = rms(groupNM_data_subj2.u);
RMS_NM3_u = rms(groupNM_data_subj3.u);
RMS_NM4_u = rms(groupNM_data_subj4.u);
RMS_NM5_u = rms(groupNM_data_subj5.u);
RMS_NM6_u = rms(groupNM_data_subj6.u);
RMS_NM7_u = rms(groupNM_data_subj7.u);
RMS_NM8_u = rms(groupNM_data_subj8.u);
RMS_NM_u = [RMS_NM1_u;RMS_NM2_u;RMS_NM3_u;RMS_NM4_u;RMS_NM5_u;RMS_NM6_u;RMS_NM7_u;RMS_NM8_u];
y2NM = median(RMS_NM_u);
constant = lsqcurvefit(@f,[0;0;0],x1,y2NM);
yfitnmu =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3));
panmu = constant(3)
p0nmu = constant(1)
Funmu = constant(2)
figure(5);
boxplot(RMS_NM_u)
hold on
plot(x1,yfitnmu,'r')
hold off
%%%%%%% Plot RMS(u) %%%%%%%
%% nm x
RMS_NM1_x = rms(groupNM_data_subj1.x);
RMS_NM2_x = rms(groupNM_data_subj2.x);
RMS_NM3_x = rms(groupNM_data_subj3.x);
RMS_NM4_x = rms(groupNM_data_subj4.x);
RMS_NM5_x = rms(groupNM_data_subj5.x);
RMS_NM6_x = rms(groupNM_data_subj6.x);
RMS_NM7_x = rms(groupNM_data_subj7.x);
RMS_NM8_x = rms(groupNM_data_subj8.x);
RMS_NM_x = [RMS_NM1_x;RMS_NM2_x;RMS_NM3_x;RMS_NM4_x;RMS_NM5_x;RMS_NM6_x;RMS_NM7_x;RMS_NM8_x];
y3NM = median(RMS_NM_x);
constant = lsqcurvefit(@f,[0;0;0],x1,y3NM);
yfitnmx =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3));
panmx = constant(3)
p0nmx = constant(1)
Funmx = constant(2)
figure(6);
boxplot(RMS_NM_x)
hold on
plot(x1,yfitnmx,'r')
hold off



%% Bode Plots
phase_M_1 = groupM_data_subj1.phase_Hp(:,1) ;
phase_M_2 = groupM_data_subj2.phase_Hp(:,1) ;
phase_M_3 = groupM_data_subj3.phase_Hp(:,1) ;
phase_M_4 = groupM_data_subj4.phase_Hp(:,1) ;
phase_M_5 = groupM_data_subj5.phase_Hp(:,1) ;
phase_M_6 = groupM_data_subj6.phase_Hp(:,1) ;
phase_M_7 = groupM_data_subj7.phase_Hp(:,1) ;
phase_M_8 = groupM_data_subj8.phase_Hp(:,1) ;
phase_M_run1 = [phase_M_1 phase_M_2 phase_M_3 phase_M_4 phase_M_5 phase_M_6 phase_M_7 phase_M_8];
phase_mean_M_1 = mean(phase_M_run1, 2) ;

mag_M_1 = groupM_data_subj1.mag_Hp(:,1) ;
mag_M_2 = groupM_data_subj2.mag_Hp(:,1) ;
mag_M_3 = groupM_data_subj3.mag_Hp(:,1) ;
mag_M_4 = groupM_data_subj4.mag_Hp(:,1) ;
mag_M_5 = groupM_data_subj5.mag_Hp(:,1) ;
mag_M_6 = groupM_data_subj6.mag_Hp(:,1) ;
mag_M_7 = groupM_data_subj7.mag_Hp(:,1) ;
mag_M_8 = groupM_data_subj8.mag_Hp(:,1) ;
mag_M_run1 = [mag_M_1 mag_M_2 mag_M_3 mag_M_4 mag_M_5 mag_M_6 mag_M_7 mag_M_8];
mag_mean_M_1 = mean(mag_M_run1, 2) ;

phase_M_1 = groupM_data_subj1.phase_Hp(:,60) ;
phase_M_2 = groupM_data_subj2.phase_Hp(:,60) ;
phase_M_3 = groupM_data_subj3.phase_Hp(:,60) ;
phase_M_4 = groupM_data_subj4.phase_Hp(:,60) ;
phase_M_5 = groupM_data_subj5.phase_Hp(:,60) ;
phase_M_6 = groupM_data_subj6.phase_Hp(:,60) ;
phase_M_7 = groupM_data_subj7.phase_Hp(:,60) ;
phase_M_8 = groupM_data_subj8.phase_Hp(:,60) ;
phase_M_60 = [phase_M_1 phase_M_2 phase_M_3 phase_M_4 phase_M_5 phase_M_6 phase_M_7 phase_M_8];
phase_mean_M_60 = mean(phase_M_60, 2) ;

mag_M_1 = groupM_data_subj1.mag_Hp(:,60) ;
mag_M_2 = groupM_data_subj2.mag_Hp(:,60) ;
mag_M_3 = groupM_data_subj3.mag_Hp(:,60) ;
mag_M_4 = groupM_data_subj4.mag_Hp(:,60) ;
mag_M_5 = groupM_data_subj5.mag_Hp(:,60) ;
mag_M_6 = groupM_data_subj6.mag_Hp(:,60) ;
mag_M_7 = groupM_data_subj7.mag_Hp(:,60) ;
mag_M_8 = groupM_data_subj8.mag_Hp(:,60) ;
mag_M_run60 = [mag_M_1 mag_M_2 mag_M_3 mag_M_4 mag_M_5 mag_M_6 mag_M_7 mag_M_8];
mag_mean_M_60 = mean(mag_M_run60, 2) ;
omega = w ;


figure(7)
subplot(2,1,1)
loglog(omega, mag_mean_M_1, 'o')
hold on
loglog(omega, mag_mean_M_60, '*')
hold off
axis(10.^[-1 1.5 -1 2])
legend('Run 1','Run 60','Location','southwest')



subplot(2,1,2)
semilogx(omega, phase_mean_M_1, 'o')
hold on
semilogx(omega, phase_mean_M_60, '*')
hold off
axis([10.^-1 10.^1.5 -360 180])
legend('Run 1','Run 60','Location','southwest')