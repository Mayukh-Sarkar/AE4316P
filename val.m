parameters = {'error' ; 'input' ; 'x'};
p0_M = [p0_e;p0_u;p0_x];
pa_M = [pa_e;pa_u;pa_x];
F_M = [Fe;F_u;Fu_x];
p0_NM = [p0_nm_e;p0_nm_u;p0_nm_x];
pa_NM = [pa_nm_e;pa_nm_u;pa_nm_x];
F_NM = [Fu_nm_e;Fu_nm_u;Fu_nm_x];
T = table(p0_M,pa_M,F_M,p0_NM,pa_NM,F_NM,'RowNames',parameters);
disp(T)
