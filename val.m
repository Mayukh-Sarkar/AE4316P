parameters = {'error' ; 'input' ; 'x';'kp';'Tl';'Tp';'zeta';'omega'};
p0_M = [p0_e;p0_u;p0_x;p0_m_kp;p0_m_Tl;p0_m_Tp;p0_m_zeta;p0_m_omega];
pa_M = [pa_e;pa_u;pa_x;pa_m_kp;pa_m_Tl;pa_m_Tp;pa_m_zeta;pa_m_omega];
F_M = [Fe;F_u;Fu_x;Fu_m_kp;Fu_m_Tp;Fu_m_Tl;Fu_m_zeta;Fu_m_omega];
p0_NM = [p0_nm_e;p0_nm_u;p0_nm_x;p0_nm_kp;p0_nm_Tl;p0_nm_Tp;p0_nm_zeta;p0_nm_omega];
pa_NM = [pa_nm_e;pa_nm_u;pa_nm_x;pa_nm_kp;pa_nm_Tl;pa_nm_Tp;pa_nm_zeta;pa_nm_omega];
F_NM = [Fu_nm_e;Fu_nm_u;Fu_nm_x;Fu_nm_kp;Fu_nm_Tp;Fu_nm_Tl;Fu_nm_zeta;Fu_nm_omega];
T = table(p0_M,pa_M,F_M,p0_NM,pa_NM,F_NM,'RowNames',parameters);
disp(T)
