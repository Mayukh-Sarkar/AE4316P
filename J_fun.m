%%%% Specify function for optimization

load('ae4316P_2021_data_group1.mat')

% N_runs = 5 ;
N_participants = 8 ;
N_parameters = 5 ;
N_runs = 60 ;
human_par_M_Kp = zeros(N_runs,N_participants) ;
human_par_M_TL = zeros(N_runs,N_participants) ;
human_par_M_tau_p = zeros(N_runs,N_participants) ;
human_par_M_zeta_nm = zeros(N_runs,N_participants) ;
human_par_M_omega_nm = zeros(N_runs,N_participants) ;
omega = w ;
omega_test = logspace(-1, 1.5, 200) ;
options = struct('MaxFunEvals', 15000,'MaxIter', 15000);

data_mag1 = zeros(N_runs, length(omega)) ;
data_phase1 = zeros(N_runs, length(omega)) ;
data_mag2 = zeros(N_runs, length(omega_test)) ;
data_phase2 = zeros(N_runs, length(omega_test)) ;

for i = 1 : N_runs
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

    for m = 1 : N_participants
        H_pe = zeros(length(w),1) ; 
        fun = @(x) 0 ;

        Phase = phase_M(:,m);
        Mag = mag_M(:,m);
        for k = 1 : length(w)
            H_pe(k) = Mag(k) * (cos(Phase(k)*pi/180) + 1j * sin(Phase(k)*pi/180)) ;
            g = @(x) (abs(H_pe(k) - x(1)*(1 + x(2)*(1j*w(k))) * exp(-1j*w(k)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*w(k) + (1j*w(k))^2))))^2;
            fun = @(x) fun(x) + g(x) ;
            data_mag1(i,k) = Mag(k) ;
            data_phase1(i,k) = Phase(k) ;
        end

        x0 = [3, 1, 0.35, 0.5, 15] ;
        % [x,fval,exitflag,output] = fminsearch(f, x0, options);
%         x = fminsearch(f,x0,options);
        lb = [0, 0, 0, 0, 5];
        ub = [100, 10, 10, 1, 30];
        x = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

        human_par_M_Kp(i,m) = x(1) ;
        human_par_M_TL(i,m) = x(2) ;
        human_par_M_tau_p(i,m) = x(3) ;
        human_par_M_zeta_nm(i,m) = x(4) ;
        human_par_M_omega_nm(i,m) = x(5) ;

        phase_out = zeros(length(omega_test),1) ;
        mag_out = zeros(length(omega_test),1) ;
        
        get_i = 100000 ;
        get_ii = 100000 ;
        get_iii = 100000 ;
        count = 0 ;
        for l = 1:length(omega_test)
            Hpe_model = x(1)*(1 + x(2)*(1j*omega_test(l))) * exp(-1j*omega_test(l)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*omega_test(l) + (1j*omega_test(l))^2));
            mag_out(l) = abs(Hpe_model) ;
            phase_out(l) = angle(Hpe_model)*180/pi ;
            if (phase_out(l) < -160) && (count == 0)
                get_i = l ;
                count = 1 ;
            end
            if (l >= get_i) && (phase_out(l) > -160)
                phase_out(l) = phase_out(l) - 360 ;
            end
            if (phase_out(l) < -480) && (count == 1)
                get_ii = l ;
                count = 2 ;
            end
            if (l >= get_ii) && (phase_out(l) > -480)
                phase_out(l) = phase_out(l) - 360 ;
            end
%             if (phase_out(l) < -880) && (count == 2)
%                 get_iii = l ;
%                 count = 3 ;
%             end
%             if (l >= get_iii) && (phase_out(l) > -880)
%                 phase_out(l) = phase_out(l) - 360 ;
%             end
            data_mag2(i,l) = mag_out(l) ;
            data_phase2(i,l) = phase_out(l) ;
        end
        
%         title = "Bode plot for run #" + num2str(i) + " from participant #" + num2str(m) ;
% %         figure(N_participants*(i-1)+m)
%         figure(i)
%         subplot(2,1,1)
%         loglog(omega, mag_M(:,m), 'o')
%         hold on
%         loglog(omega_test, mag_out)
%         hold off
%         axis(10.^[-1 1.5 -1 2])
%         legend('Experiment','Parameter estimation','Location','southwest')
% 
%         subplot(2,1,2)
%         semilogx(omega, phase_M(:,m), 'o')
%         hold on
%         semilogx(omega_test, phase_out)
%         hold off
%         axis([10.^-1 10.^1.5 -360 180])
%         legend('Experiment','Parameter estimation','Location','southwest')
%         sgtitle(title)
    end
end


title = "Bode plot for run 1 and 60 from participant 1 in Moving simulator" ;
figure(1)
set(gcf, 'Position', [100 100 700 650])
subplot(2,1,1)
loglog(omega, data_mag1(1,:),'ok')
hold on
loglog(omega, data_mag1(60,:),'*k')
hold on
loglog(omega_test, data_mag2(1,:),'-k')
hold on
loglog(omega_test, data_mag2(60,:),'--k')
hold off
axis(10.^[-1 1.5 -1 2])
xlabel("\omega [rad/s]")
ylabel("|H_{p} (j \omega)| [-]")
% legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southwest')

subplot(2,1,2)
semilogx(omega, data_phase1(1,:),'ok')
hold on
semilogx(omega, data_phase1(60,:),'*k')
hold on
loglog(omega_test, data_phase2(1,:),'-k')
hold on
loglog(omega_test, data_phase2(60,:),'--k')
hold off
axis([10.^-1 10.^1.5 -360 180])
xlabel("\omega [rad/s]")
ylabel("\angle H_{p} (j \omega) [deg]")
legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southwest')
sgtitle(title)

x1 = [1:1:60]';
x2 = [2:1:60]';

disp(human_par_M_Kp(9,:))
disp(human_par_M_TL(45,:))
disp(human_par_M_TL(48,:))
disp(human_par_M_TL(52,:))
disp(human_par_M_TL(55,:))


kp_m = mean(human_par_M_Kp,2) ;
Tl_m = mean(human_par_M_TL,2) ;
Tp_m = mean(human_par_M_tau_p,2) ;
zeta_m = mean(human_par_M_zeta_nm,2) ;
omega_m = mean(human_par_M_omega_nm,2) ;

constant = lsqcurvefit(@f,[0;0;0],x1,kp_m);
yfit_m_kp =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_m_kp = constant(3); % asymptotic value for the x which I dont know what it is...
p0_m_kp = constant(1); % inital value for the x which I dont know what it is...
Fu_m_kp = constant(2);

constant = lsqcurvefit(@f,[0;0;0],x1,Tl_m);
yfit_m_Tl =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_m_Tl = constant(3); % asymptotic value for the x which I dont know what it is...
p0_m_Tl = constant(1); % inital value for the x which I dont know what it is...
Fu_m_Tl = constant(2);

constant = lsqcurvefit(@f,[0;0;0],x1,Tp_m);
yfit_m_Tp =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_m_Tp = constant(3); % asymptotic value for the x which I dont know what it is...
p0_m_Tp = constant(1); % inital value for the x which I dont know what it is...
Fu_m_Tp = constant(2);

constant = lsqcurvefit(@f,[0;0;0],x2,zeta_m(2:end,:));
yfit_m_zeta =  constant(3) +(1-constant(2)).^x2*(constant(1)-constant(3));
pa_m_zeta = constant(3); % asymptotic value for the x which I dont know what it is...
p0_m_zeta = constant(1); % inital value for the x which I dont know what it is...
Fu_m_zeta = constant(2);

constant = lsqcurvefit(@f,[0;0;0],x2,omega_nm(2:end,:));
yfit_m_omega =  constant(3) +(1-constant(2)).^x2*(constant(1)-constant(3));
pa_m_omega = constant(3); % asymptotic value for the x which I dont know what it is...
p0_m_omega = constant(1); % inital value for the x which I dont know what it is...
Fu_m_omega = constant(2);

%% Non moving

% N_runs = 5 ;
N_participants = 1 ;
N_parameters = 5 ;
N_runs = 60 ;
human_par_NM_Kp = zeros(N_runs,N_participants) ;
human_par_NM_TL = zeros(N_runs,N_participants) ;
human_par_NM_tau_p = zeros(N_runs,N_participants) ;
human_par_NM_zeta_nm = zeros(N_runs,N_participants) ;
human_par_NM_omega_nm = zeros(N_runs,N_participants) ;
omega = w ;
omega_test = logspace(-1, 1.5, 200) ;
options = struct('MaxFunEvals', 5000,'MaxIter', 2000);

data_mag1_NM = zeros(N_runs, length(omega)) ;
data_phase1_NM = zeros(N_runs, length(omega)) ;
data_mag2_NM = zeros(N_runs, length(omega_test)) ;
data_phase2_NM = zeros(N_runs, length(omega_test)) ;

for i = 1 : N_runs
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

    for m = 1 : N_participants
        H_pe = zeros(length(w),1) ; 
        fun = @(x) 0 ;

        Phase = phase_NM(:,m);
        Mag = mag_NM(:,m);
        for k = 1 : length(w)
            H_pe(k) = Mag(k) * (cos(Phase(k)*pi/180) + 1j * sin(Phase(k)*pi/180)) ;
            g = @(x) (abs(H_pe(k) - x(1)*(1 + x(2)*(1j*w(k))) * exp(-1j*w(k)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*w(k) + (1j*w(k))^2))))^2;
            fun = @(x) fun(x) + g(x) ;
            data_mag1_NM(i,k) = Mag(k) ;
            data_phase1_NM(i,k) = Phase(k) ;
        end

        x0 = [3, 1, 0.35, 0.5, 15] ;
        % [x,fval,exitflag,output] = fminsearch(f, x0, options);
%         x = fminsearch(f,x0,options);
        lb = [0, 0, 0, 0, 5];
        ub = [100, 10, 10, 1, 30];
        x = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

        human_par_NM_Kp(i,m) = x(1) ;
        human_par_NM_TL(i,m) = x(2) ;
        human_par_NM_tau_p(i,m) = x(3) ;
        human_par_NM_zeta_nm(i,m) = x(4) ;
        human_par_NM_omega_nm(i,m) = x(5) ;

        phase_out = zeros(length(omega_test),1) ;
        mag_out = zeros(length(omega_test),1) ;
        
        get_i = 100000 ;
        get_ii = 100000 ;
        get_iii = 100000 ;
        count = 0 ;
        for l = 1:length(omega_test)
            Hpe_model = x(1)*(1 + x(2)*(1j*omega_test(l))) * exp(-1j*omega_test(l)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*omega_test(l) + (1j*omega_test(l))^2));
            mag_out(l) = abs(Hpe_model) ;
            phase_out(l) = angle(Hpe_model)*180/pi ;
            if (phase_out(l) < -160) && (count == 0)
                get_i = l ;
                count = 1 ;
            end
            if (l >= get_i) && (phase_out(l) > -160)
                phase_out(l) = phase_out(l) - 360 ;
            end
            if (phase_out(l) < -480) && (count == 1)
                get_ii = l ;
                count = 2 ;
            end
            if (l >= get_ii) && (phase_out(l) > -480)
                phase_out(l) = phase_out(l) - 360 ;
            end
%             if (phase_out(l) < -880) && (count == 2)
%                 get_iii = l ;
%                 count = 3 ;
%             end
%             if (l >= get_iii) && (phase_out(l) > -880)
%                 phase_out(l) = phase_out(l) - 360 ;
%             end
            data_mag2_NM(i,l) = mag_out(l) ;
            data_phase2_NM(i,l) = phase_out(l) ;
        end
        
%         title = "Bode plot for run #" + num2str(i) + " from participant #" + num2str(m) ;
% %         figure(N_participants*(i-1)+m)
%         figure(i)
%         subplot(2,1,1)
%         loglog(omega, mag_M(:,m), 'o')
%         hold on
%         loglog(omega_test, mag_out)
%         hold off
%         axis(10.^[-1 1.5 -1 2])
%         legend('Experiment','Parameter estimation','Location','southwest')
% 
%         subplot(2,1,2)
%         semilogx(omega, phase_M(:,m), 'o')
%         hold on
%         semilogx(omega_test, phase_out)
%         hold off
%         axis([10.^-1 10.^1.5 -360 180])
%         legend('Experiment','Parameter estimation','Location','southwest')
%         sgtitle(title)
    end
end


% title = "Bode plot for run 1 and 60 from participant 1 in Non-Moving simulator" ;
figure(2)
set(gcf, 'Position', [100 100 700 650])
subplot(2,1,1)
loglog(omega, data_mag1_NM(1,:),'ok')
hold on
loglog(omega, data_mag1_NM(60,:),'*k')
hold on
loglog(omega_test, data_mag2_NM(1,:),'-k')
hold on
loglog(omega_test, data_mag2_NM(60,:),'--k')
hold off
axis(10.^[-1 1.5 -1 2])
xlabel("\omega [rad/s]")
ylabel("|H_{p} (j \omega)| [-]")
% legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southwest')

subplot(2,1,2)
semilogx(omega, data_phase1_NM(1,:),'ok')
hold on
semilogx(omega, data_phase1_NM(60,:),'*k')
hold on
loglog(omega_test, data_phase2_NM(1,:),'-k')
hold on
loglog(omega_test, data_phase2_NM(60,:),'--k')
hold off
axis([10.^-1 10.^1.5 -360 180])
xlabel("\omega [rad/s]")
ylabel("\angle H_{p} (j \omega) [deg]")
legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southwest')
% sgtitle(title)

x1 = [1:1:60]';
x2 = [2:1:60]';

kp_nm = mean(human_par_NM_Kp,2) ;
Tl_nm = mean(human_par_NM_TL,2) ;
Tp_nm = mean(human_par_NM_tau_p,2) ;
zeta_nm = mean(human_par_NM_zeta_nm,2) ;
omega_nm = mean(human_par_NM_omega_nm,2) ;

constant = lsqcurvefit(@f,[0;0;0],x1,kp_nm);
yfit_nm_kp =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_nm_kp = constant(3); % asymptotic value for the x which I dont know what it is...
p0_nm_kp = constant(1); % inital value for the x which I dont know what it is...
Fu_nm_kp = constant(2);

constant = lsqcurvefit(@f,[0;0;0],x1,Tl_nm);
yfit_nm_Tl =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_nm_Tl = constant(3); % asymptotic value for the x which I dont know what it is...
p0_nm_Tl = constant(1); % inital value for the x which I dont know what it is...
Fu_nm_Tl = constant(2);

constant = lsqcurvefit(@f,[0;0;0],x1,Tp_nm);
yfit_nm_Tp =  constant(3) +(1-constant(2)).^x1*(constant(1)-constant(3));
pa_nm_Tp = constant(3); % asymptotic value for the x which I dont know what it is...
p0_nm_Tp = constant(1); % inital value for the x which I dont know what it is...
Fu_nm_Tp = constant(2);

constant = lsqcurvefit(@f,[0;0;0],x2,zeta_nm(2:end,:));
yfit_nm_zeta =  constant(3) +(1-constant(2)).^x2*(constant(1)-constant(3));
pa_nm_zeta = constant(3); % asymptotic value for the x which I dont know what it is...
p0_nm_zeta = constant(1); % inital value for the x which I dont know what it is...
Fu_nm_zeta = constant(2);

constant = lsqcurvefit(@f,[0;0;0],x2,omega_nm(2:end,:));
yfit_nm_omega =  constant(3) +(1-constant(2)).^x2*(constant(1)-constant(3));
pa_nm_omega = constant(3); % asymptotic value for the x which I dont know what it is...
p0_nm_omega = constant(1); % inital value for the x which I dont know what it is...
Fu_nm_omega = constant(2);


%% Plotting 

x_axis = [1:1:60];
%plotting pilot gain
figure(4)
set(gcf, 'Position', [100 100 700 650])
subplot(3,2,1)
scatter(x_axis,kp_nm,'k','o')
hold on
scatter(x_axis,kp_m,'k','*')
hold on
plot(x_axis,yfit_nm_kp,'k','LineStyle','-','LineWidth',1.2);
hold on
plot(x_axis,yfit_m_kp,'k','LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('K_p')
% ylim([-10,5])
% legend({'K_p','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plotting time lead
subplot(3,2,2)
scatter(x_axis,Tl_nm,'k','o')
hold on
scatter(x_axis,Tl_m,'k','*')
hold on
plot(x_axis,yfit_nm_Tl,'k','LineStyle','-','LineWidth',1.2);
hold on
plot(x_axis,yfit_m_Tl,'k','LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('T_L')
% legend({'\tau_L','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plotting pilot time
subplot(3,2,3)
scatter(x_axis,Tp_nm,'k','o')
hold on
scatter(x_axis,Tp_m,'k','*')
hold on
plot(x_axis,yfit_nm_Tp,'k','LineStyle','-','LineWidth',1.2);
hold on
plot(x_axis,yfit_m_Tp,'k','LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\tau_p')
% legend({'\tau_p','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plotting damping
x_new = [2:1:60];
subplot(3,2,5)
scatter(x_new,zeta_nm(2:end,:),'k','o')
hold on
scatter(x_new,zeta_m(2:end,:),'k','*')
hold on
plot(x_new,yfit_nm_zeta,'k','LineStyle','-','LineWidth',1.2);
hold on
plot(x_new,yfit_m_zeta,'k','LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\zeta_{nm}')
% legend({'\zeta_{nm}','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
%plotting omega
subplot(3,2,6)
scatter(x_new,omega_nm(2:end,:),'k','o')
hold on
scatter(x_new,omega_m(2:end,:),'k','*')
hold on
plot(x_new,yfit_nm_omega,'k','LineStyle','-','LineWidth',1.2);
hold on
plot(x_new,yfit_m_omega,'k','LineStyle','--','LineWidth',1.2);
hold off
xlabel('Trial runs')
ylabel('\omega_{nm}')
% legend({'\omega_{nm}','Learning curve','Upper bound','Lower Bound'},'Location','northeast','Orientation','vertical')
