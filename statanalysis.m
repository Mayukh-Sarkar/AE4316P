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
        disp("Null hypothesis cannot be rejected, accept H_0,accept hypothesis 2")
    else 
        disp("Null hypothesis cannot be tested , accept H_1,reject hypothesis 2")
    end
elseif Henm == 1 && Hem == 1 && Hunm == 1 && Hum == 1
    [pUe,hUe,statsUe] = ranksum(mean(testRMS_e_nm),mean(testRMS_e_m));
    [pUu,hUu,statsUu] = ranksum(mean(testRMS_u_nm),mean(testRMS_u_m));
    if hUe == 0 && hUu == 0
        disp("Null hypothesis cannot be rejected, accept H_0,accept hypothesis 2")
    else
        disp("Null hypothesis cannot be tested , accept H_1,reject hypothesis 2")
    end
end
