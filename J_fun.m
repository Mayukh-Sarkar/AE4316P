%%%% Specify function for optimization

function x = J_fun(w,Mag,Phase)
H_pe = zeros(length(w),1) ; 
f = @(x) 0 ;
for i = 1 : length(w)
    H_pe(i) = Mag(i) * (cos(Phase(i)*pi/180) + 1j * sin(Phase(i)*pi/180)) ;
    g = @(x)(H_pe(i) - x(1)*(1 + x(2)*(1j*w(i))) * exp(-1j*w(i)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*w(i) + (1j*w(i))^2)))^2;
    f = @(x) f(x) + g(x) ;
end

x0 = [3, 1, 0.35, 0.5, 15] ;
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
% x0 = [4, 2, 0.35, 0.25, 10];
x = fminsearch(f, x0);

disp(f)
end

%f = @f + @(K_v, T_l, tau_p, zeta_nm, omega_nm)(H_pe(i) - K_v*(1 + T_l*(1j*w(i))) * exp(-1j*w(i)*tau_p) * (omega_nm^2/(omega_nm^2 + 2*zeta_nm*omega_nm*1j*w(i) + (1j*w(i))^2)));


