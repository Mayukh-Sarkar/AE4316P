%%%% Specify function for optimization

load('ae4316P_2021_data_group1.mat')

N_runs = 5 ;
N_participants = 8 ;
N_parameters = 5 ;
% N_runs = 60 ;
human_par_M_Kp = zeros(N_runs,N_participants) ;
human_par_M_TL = zeros(N_runs,N_participants) ;
human_par_M_tau_p = zeros(N_runs,N_participants) ;
human_par_M_zeta_nm = zeros(N_runs,N_participants) ;
human_par_M_omega_nm = zeros(N_runs,N_participants) ;
omega = w ;
options = struct('MaxFunEvals', 2000,'MaxIter', 1000);

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
        f = @(x) 0 ;

        Phase = phase_M(:,m);
        Mag = mag_M(:,m);
        for k = 1 : length(w)
            H_pe(k) = Mag(k) * (cos(Phase(k)*pi/180) + 1j * sin(Phase(k)*pi/180)) ;
            g = @(x) (abs(H_pe(k) - x(1)*(1 + x(2)*(1j*w(k))) * exp(-1j*w(k)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*w(k) + (1j*w(k))^2))))^2;
            f = @(x) f(x) + g(x) ;
        end

        x0 = [3, 1, 0.35, 0.5, 15] ;
        % [x,fval,exitflag,output] = fminsearch(f, x0, options);
        x = fminsearch(f,x0,options);

        human_par_M_Kp(i,m) = x(1) ;
        human_par_M_TL(i,m) = x(2) ;
        human_par_M_tau_p(i,m) = x(3) ;
        human_par_M_zeta_nm(i,m) = x(4) ;
        human_par_M_omega_nm(i,m) = x(5) ;


        omega_test = logspace(-1, 1.5, 200) ;
        phase_out = zeros(length(omega_test),1) ;
        mag_out = zeros(length(omega_test),1) ;

        for l = 1:length(omega_test)
            Hpe_model = x(1)*(1 + x(2)*(1j*omega_test(l))) * exp(-1j*omega_test(l)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*omega_test(l) + (1j*omega_test(l))^2));
            mag_out(l) = abs(Hpe_model) ;
            phase_out(l) = angle(Hpe_model)*180/pi ;
        end
        
        title = "Bode plot for run #" + num2str(i) + " from participant #" + num2str(m) ;
        figure(7+i)
        subplot(2,1,1)
        loglog(omega, mag_M(:,m), 'o')
        hold on
        loglog(omega_test, mag_out)
        hold off
        axis(10.^[-1 1.5 -1 2])
        legend('Experiment','Parameter estimation','Location','southwest')

        subplot(2,1,2)
        semilogx(omega, phase_M(:,m), 'o')
        hold on
        semilogx(omega_test, phase_out)
        hold off
        axis([10.^-1 10.^1.5 -360 180])
        legend('Experiment','Parameter estimation','Location','southwest')
        sgtitle(title)
    end
end