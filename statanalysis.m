testRMS_e_nm = RMS_NM_e(:,56:60);
testRMS_e_m = RMS_M_e(:,56:60);
testRMS_u_nm = RMS_NM_u(:,56:60);
testRMS_u_m = RMS_M_u(:,56:60);
alpha = 0.05;
H_0 = 0;
H_1 = 1;
[Henm , penm, Wenm] = swtest(mean(testRMS_e_nm),0.05);
[Hem , pem, Wem] = swtest(mean(testRMS_e_m),0.05);
[Hunm , punm, Wunmm] = swtest(mean(testRMS_u_nm),0.05);
[Hum , pum, Wum] = swtest(mean(testRMS_u_m),0.05);

if Henm == 0 && Hem == 0 && Hunm == 0 && Hum == 0
    disp('H_0 is accepted , samples are normally distributed ')
    [Htenm , ptenm] = ttest(mean(testRMS_e_nm));
    [Htem , ptem] = ttest(mean(testRMS_e_m));
    [Htunm , ptunm] = ttest(mean(testRMS_u_nm));
    [Htum , ptum] = ttest(mean(testRMS_u_m));
    if Htenm == 0 && Htem == 0 && Htunm == 0 && Htum == 0
        disp('H_0 is true, two distributions are identitcal accept hypothesis 2')
    else 
        disp('H_1 is true, two distributions are non- reject hypothesis 2')
    end
end