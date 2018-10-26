function[global_error,t_vals,approx_sol]=aprox_ode(f,a,b,initial_val,N)

%h = (b-a)./N; %1
h = .2; %2
t_vals = a:h:b;
approx_sol = (zeros(1, N+1));
approx_sol(1) = initial_val; 

%s = [.1 .2 .3 .4]; %2a
s = [.2 .22 .24 .26]; %2b


for j = 1:4
    for i = 2:N+1 
        k1=h*feval('f',t_vals(i-1),approx_sol(i-1), s(j));
        k2=h*feval('f',t_vals(i-1)+(h./2),approx_sol(i-1)+(k1./2), s(j));
        k3=h*feval('f',t_vals(i-1)+(h./2),approx_sol(i-1)+(k2./2), s(j));
        k4=h*feval('f',t_vals(i-1)+h,approx_sol(i-1)+k3, s(j)); 
        approx_sol(i) =  approx_sol(i-1) + (1/6)*(k1+(2.*k2)+(2*k3)+k4);
    end
    hold on
    plot(approx_sol, 'LineWidth', 2)
    approx_sol(end)
end
legend('S1', 'S2', 'S3', 'S4')

global_error=abs(-b-(1/(exp(b)+1))-approx_sol(end))
t_vals;
approx_sol(end)


end