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
% figure(1)
% scatter(x,Med_M_e,'*')
% hold on
% x = [1:1:60];
y_m_e = Med_M_e;
constant = lsqcurvefit(@f,[0;0;0],x,y_m_e); % curve fit for the mean error across 60 runs
pa_e = constant(3); % asymptotic values
p0_e = constant(1); % initial error vaue
Fe = constant(2); % learning curve slope

yfit_e =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3)); % curve fitting values acrros 60 runs

% hold on
% plot(x,yfit_e,'r');
% hold off


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
x = [1:1:60];
med_u = median(RMS_M_u);
constant = lsqcurvefit(@f,[0;0;0],x,med_u); % curve fitting function for input u
pa_u = constant(3); % asymptotic value of the input curve
p0_u = constant(1); % initial error value for the input
F_u = constant(2);% learning slope 
yfit_u =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3)); % value for the fit

figure(2);
boxplot(RMS_M_u)
hold on
plot(x,yfit_u,'r');
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
x = [1:1:60];
med_x = median(RMS_M_x);
constant = lsqcurvefit(@f,[0;0;0],x,med_x); % leanring curve model
pa_x = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_x = constant(1); % initial value for the x which I dont know what it is...Hail lucifer
Fu_x = constant(2); % learning curve slope for the x which I dont know what it is...Hail lucifer
yfit_x =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3)); % fittig curve value for each trial
 
figure(3);
boxplot(RMS_M_x)
hold on
plot(x,yfit_x,'r');
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
y_nm_e = median(RMS_NM_e);

constant = lsqcurvefit(@f,[0;0;0],x,y_nm_e); % fitting curve
yfit_nm_e =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3)); % valyes of the fitting curve
pa_nm_e = constant(3); % asymptotic values for non motion error
p0_nm_e = constant(1); % initial value for non motion error
Fu_nm_e = constant(2); % learning curve rate
figure(4);
boxplot(RMS_NM_e)
hold on

plot(x,yfit_nm_e,'r');
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
y_nm_u = median(RMS_NM_u);
constant = lsqcurvefit(@f,[0;0;0],x,y_nm_u); % learning curve
yfit_nm_u =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3)); % learning curve values for the pilot input
pa_nm_u = constant(3); % asymptotic valuue for thge input curve
p0_nm_u = constant(1); % inital value of the input curve
Fu_nm_u = constant(2);% learning rate
figure(5);
boxplot(RMS_NM_u)
hold on
plot(x,yfit_nm_u,'r')
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
y_NM_x = median(RMS_NM_x);
constant = lsqcurvefit(@f,[0;0;0],x,y_NM_x);
yfit_nm_x =  constant(3) +(1-constant(2)).^x*(constant(1)-constant(3));
pa_nm_x = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_nm_x = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_nm_x = constant(2); % leanring rate  for the x which I dont know what it is...Hail lucifer
figure(6);
boxplot(RMS_NM_x)
hold on
plot(x,yfit_nm_x,'r')
hold off



%% Bode Plots
% phase_M_1 = groupM_data_subj1.phase_Hp(:,1) ;
% phase_M_2 = groupM_data_subj2.phase_Hp(:,1) ;
% phase_M_3 = groupM_data_subj3.phase_Hp(:,1) ;
% phase_M_4 = groupM_data_subj4.phase_Hp(:,1) ;
% phase_M_5 = groupM_data_subj5.phase_Hp(:,1) ;
% phase_M_6 = groupM_data_subj6.phase_Hp(:,1) ;
% phase_M_7 = groupM_data_subj7.phase_Hp(:,1) ;
% phase_M_8 = groupM_data_subj8.phase_Hp(:,1) ;
% phase_M_run1 = [phase_M_1 phase_M_2 phase_M_3 phase_M_4 phase_M_5 phase_M_6 phase_M_7 phase_M_8];
% phase_mean_M_1 = mean(phase_M_run1, 2) ;
% 
% mag_M_1 = groupM_data_subj1.mag_Hp(:,1) ;
% mag_M_2 = groupM_data_subj2.mag_Hp(:,1) ;
% mag_M_3 = groupM_data_subj3.mag_Hp(:,1) ;
% mag_M_4 = groupM_data_subj4.mag_Hp(:,1) ;
% mag_M_5 = groupM_data_subj5.mag_Hp(:,1) ;
% mag_M_6 = groupM_data_subj6.mag_Hp(:,1) ;
% mag_M_7 = groupM_data_subj7.mag_Hp(:,1) ;
% mag_M_8 = groupM_data_subj8.mag_Hp(:,1) ;
% mag_M_run1 = [mag_M_1 mag_M_2 mag_M_3 mag_M_4 mag_M_5 mag_M_6 mag_M_7 mag_M_8];
% mag_mean_M_1 = mean(mag_M_run1, 2) ;
% 
% phase_M_1 = groupM_data_subj1.phase_Hp(:,60) ;
% phase_M_2 = groupM_data_subj2.phase_Hp(:,60) ;
% phase_M_3 = groupM_data_subj3.phase_Hp(:,60) ;
% phase_M_4 = groupM_data_subj4.phase_Hp(:,60) ;
% phase_M_5 = groupM_data_subj5.phase_Hp(:,60) ;
% phase_M_6 = groupM_data_subj6.phase_Hp(:,60) ;
% phase_M_7 = groupM_data_subj7.phase_Hp(:,60) ;
% phase_M_8 = groupM_data_subj8.phase_Hp(:,60) ;
% phase_M_60 = [phase_M_1 phase_M_2 phase_M_3 phase_M_4 phase_M_5 phase_M_6 phase_M_7 phase_M_8];
% phase_mean_M_60 = mean(phase_M_60, 2) ;
% 
% mag_M_1 = groupM_data_subj1.mag_Hp(:,60) ;
% mag_M_2 = groupM_data_subj2.mag_Hp(:,60) ;
% mag_M_3 = groupM_data_subj3.mag_Hp(:,60) ;
% mag_M_4 = groupM_data_subj4.mag_Hp(:,60) ;
% mag_M_5 = groupM_data_subj5.mag_Hp(:,60) ;
% mag_M_6 = groupM_data_subj6.mag_Hp(:,60) ;
% mag_M_7 = groupM_data_subj7.mag_Hp(:,60) ;
% mag_M_8 = groupM_data_subj8.mag_Hp(:,60) ;
% mag_M_run60 = [mag_M_1 mag_M_2 mag_M_3 mag_M_4 mag_M_5 mag_M_6 mag_M_7 mag_M_8];
% mag_mean_M_60 = mean(mag_M_run60, 2) ;
% omega = w ;
% 
% 
% figure(7)
% subplot(2,1,1)
% loglog(omega, mag_mean_M_1, 'o')
% hold on
% loglog(omega, mag_mean_M_60, '*')
% hold on
% axis(10.^[-1 1.5 -1 2])
% legend('Run 1','Run 60','Location','southwest')
% constant = lsqcurvefit(@f,[0;0;0],omega,mag_mean_M_60);
% %test fit
% yfitM =  constant(3) +(1-constant(2)).^omega*(constant(1)-constant(3));
% paM = constant(3)
% p0M = constant(1)
% FuM = constant(2)
% 
% loglog(omega,yfitM,'r')
% hold off
% 
% subplot(2,1,2)
% semilogx(omega, phase_mean_M_1, 'o')
% hold on
% semilogx(omega, phase_mean_M_60, '*')
% hold on
% axis([10.^-1 10.^1.5 -360 180])
% legend('Run 1','Run 60','Location','southwest')
% 
% %test fit
% yfitM =  constant(3) +(1-constant(2)).^omega*(constant(1)-constant(3));
% paM = constant(3);
% p0M = constant(1);
% FuM = constant(2);
% 
% loglog(omega,yfitM,'r')
% hold off


%val % creating the table for the  fitting values; update the parameters in val.m


%% Human Pilot M model parameters

% options = optimset('Display','iter','PlotFcns',@optimplotfval);

% N_participants = 5 ;
N_participants = 60 ;
N_par = 5 ;
human_par_M = zeros(N_participants, N_par) ;

for i = 1 : N_participants
    phase_M_1 = groupM_data_subj1.phase_Hp(:,i) ;
    phase_M_2 = groupM_data_subj2.phase_Hp(:,i) ;
    phase_M_3 = groupM_data_subj3.phase_Hp(:,i) ;
    phase_M_4 = groupM_data_subj4.phase_Hp(:,i) ;
    phase_M_5 = groupM_data_subj5.phase_Hp(:,i) ;
    phase_M_6 = groupM_data_subj6.phase_Hp(:,i) ;
    phase_M_7 = groupM_data_subj7.phase_Hp(:,i) ;
    phase_M_8 = groupM_data_subj8.phase_Hp(:,i) ;
    phase_M = [phase_M_1 phase_M_2 phase_M_3 phase_M_4 phase_M_5 phase_M_6 phase_M_7 phase_M_8];
    phase_mean_M = mean(phase_M, 2) ;

    mag_M_1 = groupM_data_subj1.mag_Hp(:,i) ;
    mag_M_2 = groupM_data_subj2.mag_Hp(:,i) ;
    mag_M_3 = groupM_data_subj3.mag_Hp(:,i) ;
    mag_M_4 = groupM_data_subj4.mag_Hp(:,i) ;
    mag_M_5 = groupM_data_subj5.mag_Hp(:,i) ;
    mag_M_6 = groupM_data_subj6.mag_Hp(:,i) ;
    mag_M_7 = groupM_data_subj7.mag_Hp(:,i) ;
    mag_M_8 = groupM_data_subj8.mag_Hp(:,i) ;
    mag_M = [mag_M_1 mag_M_2 mag_M_3 mag_M_4 mag_M_5 mag_M_6 mag_M_7 mag_M_8];
    mag_mean_M = mean(mag_M, 2) ;

    H_pe = zeros(length(w),1) ; 
    f = @(x) 0 ;

    Phase = phase_mean_M;
    Mag = mag_mean_M;
    for k = 1 : length(w)
        H_pe(k) = Mag(k) * (cos(Phase(k)*pi/180) + 1j * sin(Phase(k)*pi/180)) ;
        g = @(x) (abs(H_pe(k) - x(1)*(1 + x(2)*(1j*w(k))) * exp(-1j*w(k)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*w(k) + (1j*w(k))^2))))^2;
        f = @(x) f(x) + g(x) ;
    end

    x0 = [3, 1, 0.35, 0.5, 15] ;
    % [x,fval,exitflag,output] = fminsearch(f, x0, options);
    x = fminsearch(f,x0);
    
    human_par_M(i,1) = x(1) ;
    human_par_M(i,2) = x(2) ;
    human_par_M(i,3) = x(3) ;
    human_par_M(i,4) = x(4) ;
    human_par_M(i,5) = x(5) ;
    
    
    omega_test = logspace(-1, 1.5, 200) ;
    phase_out = zeros(length(omega_test),1) ;
    mag_out = zeros(length(omega_test),1) ;

    for l = 1:length(omega_test)
        Hpe_model = x(1)*(1 + x(2)*(1j*omega_test(l))) * exp(-1j*omega_test(l)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*omega_test(l) + (1j*omega_test(l))^2));
        mag_out(l) = abs(Hpe_model) ;
        phase_out(l) = angle(Hpe_model)*180/pi ;
    end

%     title = "Bode plot for run #" + num2str(i) ;
%     figure(7+i)
%     subplot(2,1,1)
%     loglog(omega, mag_mean_M, 'o')
%     hold on
%     loglog(omega_test, mag_out)
%     hold off
%     axis(10.^[-1 1.5 -1 2])
%     legend('Experiment','Parameter estimation','Location','southwest')
%     
%     subplot(2,1,2)
%     semilogx(omega, phase_mean_M, 'o')
%     hold on
%     semilogx(omega_test, phase_out)
%     hold off
%     axis([10.^-1 10.^1.5 -360 180])
%     legend('Experiment','Parameter estimation','Location','southwest')
%     sgtitle(title)

end


%% Human Pilot NM model parameters

% options = optimset('Display','iter','PlotFcns',@optimplotfval);

% N_participants = 5 ;
N_participants = 60 ;
N_par = 5 ;
human_par_NM = zeros(N_participants, N_par) ;

for i = 1 : N_participants
    phase_NM_1 = groupNM_data_subj1.phase_Hp(:,i) ;
    phase_NM_2 = groupNM_data_subj2.phase_Hp(:,i) ;
    phase_NM_3 = groupNM_data_subj3.phase_Hp(:,i) ;
    phase_NM_4 = groupNM_data_subj4.phase_Hp(:,i) ;
    phase_NM_5 = groupNM_data_subj5.phase_Hp(:,i) ;
    phase_NM_6 = groupNM_data_subj6.phase_Hp(:,i) ;
    phase_NM_7 = groupNM_data_subj7.phase_Hp(:,i) ;
    phase_NM_8 = groupNM_data_subj8.phase_Hp(:,i) ;
    phase_NM = [phase_NM_1 phase_NM_2 phase_NM_3 phase_NM_4 phase_NM_5 phase_NM_6 phase_NM_7 phase_NM_8];
    phase_mean_NM = mean(phase_NM, 2) ;

    mag_NM_1 = groupNM_data_subj1.mag_Hp(:,i) ;
    mag_NM_2 = groupNM_data_subj2.mag_Hp(:,i) ;
    mag_NM_3 = groupNM_data_subj3.mag_Hp(:,i) ;
    mag_NM_4 = groupNM_data_subj4.mag_Hp(:,i) ;
    mag_NM_5 = groupNM_data_subj5.mag_Hp(:,i) ;
    mag_NM_6 = groupNM_data_subj6.mag_Hp(:,i) ;
    mag_NM_7 = groupNM_data_subj7.mag_Hp(:,i) ;
    mag_NM_8 = groupNM_data_subj8.mag_Hp(:,i) ;
    mag_NM = [mag_NM_1 mag_NM_2 mag_NM_3 mag_NM_4 mag_NM_5 mag_NM_6 mag_NM_7 mag_NM_8];
    mag_mean_NM = mean(mag_NM, 2) ;

    H_pe = zeros(length(w),1) ; 
    f = @(x) 0 ;

    Phase = phase_mean_NM;
    Mag = mag_mean_NM;
    for k = 1 : length(w)
        H_pe(k) = Mag(k) * (cos(Phase(k)*pi/180) + 1j * sin(Phase(k)*pi/180)) ;
        g = @(x) (abs(H_pe(k) - x(1)*(1 + x(2)*(1j*w(k))) * exp(-1j*w(k)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*w(k) + (1j*w(k))^2))))^2;
        f = @(x) f(x) + g(x) ;
    end

    x0 = [3, 1, 0.35, 0.5, 15] ;
    % [x,fval,exitflag,output] = fminsearch(f, x0, options);
    x = fminsearch(f,x0);
    
    human_par_NM(i,1) = x(1) ;
    human_par_NM(i,2) = x(2) ;
    human_par_NM(i,3) = x(3) ;
    human_par_NM(i,4) = x(4) ;
    human_par_NM(i,5) = x(5) ;
    
    
    omega_test = logspace(-1, 1.5, 200) ;
    phase_out = zeros(length(omega_test),1) ;
    mag_out = zeros(length(omega_test),1) ;

    for l = 1:length(omega_test)
        Hpe_model = x(1)*(1 + x(2)*(1j*omega_test(l))) * exp(-1j*omega_test(l)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*omega_test(l) + (1j*omega_test(l))^2));
        mag_out(l) = abs(Hpe_model) ;
        phase_out(l) = angle(Hpe_model)*180/pi ;
    end

%     title = "Bode plot for run #" + num2str(i) ;
%     figure(7+i)
%     subplot(2,1,1)
%     loglog(omega, mag_mean_M, 'o')
%     hold on
%     loglog(omega_test, mag_out)
%     hold off
%     axis(10.^[-1 1.5 -1 2])
%     legend('Experiment','Parameter estimation','Location','southwest')
%     
%     subplot(2,1,2)
%     semilogx(omega, phase_mean_M, 'o')
%     hold on
%     semilogx(omega_test, phase_out)
%     hold off
%     axis([10.^-1 10.^1.5 -360 180])
%     legend('Experiment','Parameter estimation','Location','southwest')
%     sgtitle(title)

end

%% nomotion group
x1 = [1:1:60]';
kp_nm = human_par_NM(:,1);
Tl_nm = human_par_NM(:,2);
Tp_nm  = human_par_NM(:,3);
zeta_nm = human_par_NM(:,4);
omega_nm = human_par_NM(:,5);

constant = lsqcurvefit(@f,[0;0;0],x1,kp_nm);
yfit_nm_kp =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_nm_kp = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_nm_kp = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_nm_kp = constant(2);

figure(7);
plot(x1,kp_nm,'*')
hold on
plot(x1,yfit_nm_kp,'r')
hold off

constant = lsqcurvefit(@f,[0;0;0],x1,Tl_nm);
yfit_nm_Tl =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_nm_Tl = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_nm_Tl = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_nm_Tl = constant(2);

figure(8);
plot(x1,Tl_nm,'*')
hold on
plot(x1,yfit_nm_Tl,'r')
hold off

figure(9);
constant = lsqcurvefit(@f,[0;0;0],x1,Tp_nm);
yfit_nm_Tp =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_nm_Tp = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_nm_Tp = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_nm_Tp = constant(2);
plot(x1,Tp_nm,'*')
hold on
plot(x1,yfit_nm_Tp,'r')
hold off
x2 = [2:1:60]';
figure(10);
constant = lsqcurvefit(@f,[0;0;0],x2,zeta_nm(2:end,:));
yfit_nm_zeta =  constant(3) +(1-constant(2)).^x2*(constant(1)-constant(3));
pa_nm_zeta = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_nm_zeta = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_nm_zeta = constant(2);
plot(x2,zeta_nm(2:end,:),'*')
hold on
plot(x2,yfit_nm_zeta,'r')
hold off

figure(11);
constant = lsqcurvefit(@f,[0;0;0],x2,omega_nm(2:end,:));
yfit_nm_omega =  constant(3) +(1-constant(2)).^x2*(constant(1)-constant(3));
pa_nm_omega = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_nm_omega = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_nm_omega = constant(2);
plot(x2,omega_nm(2:end,:),'*')
hold on
plot(x2,yfit_nm_omega,'r')
hold off

kp_m = human_par_M(:,1);
Tl_m = human_par_M(:,2);
Tp_m  = human_par_M(:,3);
zeta_m = human_par_M(:,4);
omega_m = human_par_M(:,5);

constant = lsqcurvefit(@f,[0;0;0],x1,kp_m);
yfit_m_kp =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_m_kp = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_m_kp = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_m_kp = constant(2);

figure(12);
plot(x1,kp_m,'*')
hold on
plot(x1,yfit_m_kp,'r')
hold off

constant = lsqcurvefit(@f,[0;0;0],x1,Tl_m);
yfit_m_Tl =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_m_Tl = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_m_Tl = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_m_Tl = constant(2);

figure(13);
plot(x1,Tl_m,'*')
hold on
plot(x1,yfit_m_Tl,'r')
hold off

figure(14);
constant = lsqcurvefit(@f,[0;0;0],x1,Tp_m);
yfit_m_Tp =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_m_Tp = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_m_Tp = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_m_Tp = constant(2);
plot(x1,Tp_m,'*')
hold on
plot(x1,yfit_m_Tp,'r')
hold off

figure(15);
constant = lsqcurvefit(@f,[0;0;0],x2,zeta_m(2:end,:));
yfit_m_zeta =  constant(3) +(1-constant(2)).^x2*(constant(1)-constant(3));
pa_m_zeta = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_m_zeta = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_m_zeta = constant(2);
plot(x2,zeta_m(2:end,:),'*')
hold on
plot(x2,yfit_m_zeta,'r')
hold off

figure(16);
constant = lsqcurvefit(@f,[0;0;0],x2,omega_m(2:end,:));
yfit_m_omega =  constant(3) +(1-constant(2)).^x2*(constant(1)-constant(3));
pa_m_omega = constant(3); % asymptotic value for the x which I dont know what it is...Hail lucifer
p0_m_omega = constant(1); % inital value for the x which I dont know what it is...Hail lucifer
Fu_m_omega = constant(2);
plot(x2,omega_m(2:end,:),'*')
hold on
plot(x2,yfit_m_omega,'r')
hold off

val
