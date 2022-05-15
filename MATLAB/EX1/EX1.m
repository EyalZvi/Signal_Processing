clc;
clear all;
close all;
%_________________________________________________%
% Q1 - A
% Generate |x(t)| and plot graph to screen

t = 0.2: 0.01 : 3;
wm = 3*pi;
x1= x(t,wm);
plot(t,abs(x1),Color=[1,0,0]);
title("Signal x(t) - Time Domain:");
xlabel("t [sec]");
ylabel("x(t) [V]");
figure;
%_________________________________________________%
% Q1 - B
% Generate |X(ω)| and plot graph to screen

w = -17*pi: 0.1 : 17*pi;
X1 = X(wm,w,0,0);
plot(w,abs(X1),Color=[1,0,0]);
title("Signal X(ω) - Frequency Domain:");
xlabel("ω [rad/sec]");
ylabel("| X(ω) |");
figure;
%_________________________________________________%
% Q1 - C
% Generate x_(ZOH)(t) and plot aside x(t)

Ts = 1/15;
tn = 0.2:Ts:3;
xn = x(tn,wm);
stairs(tn,xn, Color=[0,0,0]);
hold on
plot(t,x1, Color=[1,0,0]);
title("x(t) and x_Z_O_H(t) - Time Domain");
xlabel("t=T_s[sec]");
ylabel("x(t) [V]");
legend('x_Z_O_H(t)','x(t)');
hold off
figure;
%_________________________________________________%
% Q1 - D
% Generate |X_(ZOH)(ω)| and plot graph to screen

ws = 10*wm;
sinc = sin(pi.*w/ws)./(pi.*w/ws);
H=Ts.*sinc.*exp(-1i.*w.*Ts/2);
Xs=(1/Ts).*(X1+X(wm,w,ws,1)+X(wm,w,ws,-1));
X_zoh = Xs.*H;
plot(w,abs(X_zoh),Color=[1,0,0]);
xlabel("ω [rad/sec]");
ylabel("| X_Z_O_H(ω) |");
title("Signal X_Z_O_H(ω) - Frequency Domain");
figure;
%_________________________________________________%
% Q1 - E
% Reconstruction of x(t) with ZOH

% Generate H(ω) - given filter
H_rec = exp((1i.*pi.*w)./ws)./sinc.*rectangularPulse(-ws/2,ws/2,w);
X_rec=X_zoh.*H_rec;

% Inverse Fourier Transform
vec=zeros(size(t));
i=1;
for ti=0.2:0.01:3
    phi=X_rec.*exp(1i.*w.*ti);
    vec(i)=(1/(2*pi))*trapz(w,phi);
    i=i+1;
end

plot(t,x1, Color=[1,0,0]);
hold on
plot(t,vec,'.', Color=[0,0,1]);
legend('x(t)','x_r_e_c(t)');
title('Signals x_r_e_c(t) and x(t) - Time Domain');
xlabel('t [sec]');
ylabel('x(t) [V]');
hold off
figure;
%_________________________________________________%
% Q1 - F
% Repeat reconstruction process with ωs = 9ωm

ws = 9*wm;
Ts2 = 2*pi/ws;
sinc = sin(pi.*w/ws)./(pi.*w/ws);
H=Ts2.*sinc.*exp(-1i.*w.*Ts2/2);
Xs=(1/Ts2).*(X1+X(wm,w,ws,1)+X(wm,w,ws,-1));
X_zoh = Xs.*H;
plot(w,abs(X_zoh),Color=[1,0,0]);
xlabel("ω [rad/sec]");
ylabel("| X_Z_O_H(ω) |");
title("Signal X_Z_O_H(ω) (sampled with ω_s=9ω_m) - Frequency Domain");
figure;
%_______________________________%

H_rec = exp((1i.*pi.*w)./ws)./sinc.*rectangularPulse(-ws/2,ws/2,w);
X_rec=X_zoh.*H_rec;

vec=zeros(size(t));
% Inverse Fourier Transform
i=1;
for ti=0.2:0.01:3
    phi=X_rec.*exp(1i.*w.*ti);
    vec(i)=(1/(2*pi))*trapz(w,phi);
    i=i+1;
end

plot(t,x1, Color=[1,0,0]);
hold on
plot(t,vec,'-', Color=[0,0,1]);
legend('x(t)','x_r_e_c(t)');
title('Signals x(t) and x_r_e_c(t) (sampled with ω_s=9ω_m) - Time Domain');
xlabel('t [sec]');
ylabel('x(t) [V]');
hold off
figure;
%_________________________________________________%
% Q2 - A
% Even sampling of given periodic signal
% 15 Samples

wA = 7.*pi;
wB = 4.*pi;
T = 2;
w0 = 2*pi/T;
t = 0: T/1000 : T;
t_s = 0: T/15: (1-1/15)*T;
x_2 = x2(t,wA,wB);
x_2_s = x2(t_s,wA,wB);
plot(t,x_2,Color=[0,0,0]);
xlabel('t [sec]');
ylabel('x(t) [V]');
title('Signals x(t) and x_s(t) - Time Domain');
hold on
plot(t_s,x_2_s,'*',Color=[1,0,0]);
legend('x(t)','x_s(t)');
hold off
figure;
%_________________________________________________%
% Q2 - B
% Compute (by sampling) Fourier coeffiecients vector (a = F^(-1) * x)

m=(-7:7);
F_rec=exp(1i.*w0.*t_s'.*m); 
a=inv(F_rec)*x_2_s.';
condF=cond(F_rec);
%_________________________________________________%
% Q2 - C
% Reconstruction of x(t) from Fourier coeffiecients vector
% 15 Samples

F_Matrix=exp(1i.*w0.*t'.*m); 
x_vec = F_Matrix*a;

plot(t,x_2,Color=[0,0,0]);
hold on
plot(t,x_vec,'-.',Color=[1,0,0]);
title("Signals x(t) and x_r_e_c (t) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
legend('x(t)','x_r_e_c(t)');
hold off
figure;
%_________________________________________________%
% Q2 - D
% Compute and reconstruction - Random sampling x_s(t)
% 15 Samples

clear all
clc
wA = 7*pi;
wB = 4*pi;
T = 2;
w0 = 2*pi/T;
t = 0: T/1000 : T;
t_s_rand = (T-0.0000000001)*rand(1,15); 
x_2 = x2(t,wA,wB);
x_2_s_rand = x2(t_s_rand,wA,wB);
plot(t,x_2,Color=[0,0,0]);
hold on
plot(t_s_rand,x_2_s_rand,'*',Color=[1,0,0]);
title("Signals x(t) and random sampled x_s(t) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
legend('x(t)','x_s(t)');
hold off
figure;
%_______________________________%
m=(-7:7);
F_rec=exp(1i.*w0.*t_s_rand'.*m); 
a=inv(F_rec)*x_2_s_rand.';
condF=cond(F_rec);
%_______________________________%
F_Matrix=exp(1i.*w0.*t'.*m); 
x_vec = F_Matrix*a;

plot(t,x_2,Color=[0,0,0]);
hold on
plot(t,x_vec,'-.',Color=[1,0,0]);
title("Signals x(t) and x_r_e_c_2(t) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
legend('x(t)','x_r_e_c_2(t)');
hold off
figure;
%_________________________________________________%
% Q2 - E
% Compute and reconstruction - Even sampling with random error
% 15 Samples

clear all
clc
wA = 7.*pi;
wB = 4.*pi;
T = 2;
w0 = 2*pi/T;
t_s = 0: T/15: (1-1/15)*T;
t = 0: T/1000 : T;
t_s_rand = (T-0.0000000001)*rand(1,15);
tn_even_s=t_s+0.01*rand(1,15);
x_2 = x2(t,wA,wB);
x2_even_s=x2(t_s,wA,wB);
plot(t,x_2,Color=[0,0,0]);
hold on
title("Signals x(t) and even sampled x_s(t) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
plot(t_s,x2_even_s,'*',Color=[1,0,0]);
legend('x(t)','x_s(t)');
hold off
figure;
%_______________________________%
m=(-7:7);
F_even=exp(1i.*w0.*tn_even_s'.*m); 
a_even=inv(F_even)*x2_even_s.';
condF_even=cond(F_even);
%_______________________________%
F_Matrix_even=exp(1i.*w0.*t'.*m); 
x_vec_even = F_Matrix_even*a_even;

plot(t,x_2,Color=[0,0,0]);
hold on
plot(t,x_vec_even,'-.',Color=[1,0,0]);
title("Signals x(t) and x_r_e_c_3(t) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
legend('x(t)','x_r_e_c_3(t)');
hold off
figure;
%_______________________________%
% Compute and reconstruction - Uneven sampling with random error
% 15 Samples

tn_uneven_s=t_s_rand+0.01*rand(1,15);
x2_uneven_s=x2(t_s_rand,wA,wB);
plot(t,x_2,Color=[0,0,0]);
hold on
title("Signals x(t) and uneven sampled x_s(t) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
plot(t_s_rand,x2_uneven_s,'*',Color=[1,0,0]);
legend('x(t)','x_s(t)');
hold off
figure;
%_______________________________%
F_uneven=exp(1i.*w0.*tn_uneven_s'.*m); 
a_uneven=inv(F_uneven)*x2_uneven_s.';
condF_uneven=cond(F_uneven);
%_______________________________%
F_Matrix_uneven=exp(1i.*w0.*t'.*m); 
x_vec_uneven = F_Matrix_uneven*a_uneven;

plot(t,x_2,Color=[0,0,0]);
hold on
plot(t,x_vec_uneven,'-.',Color=[1,0,0]);
title("Signals x(t) and x_r_e_c_4(t) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
legend('x(t)','x_r_e_c_4(t)');
hold off
figure;
%_________________________________________________%
% Q2 - F
% Compute and reconstruction - Uneven sampling with random error
% 40 Samples

clear all
clc
wA = 7.*pi;
wB = 4.*pi;
T = 2;
w0 = 2*pi/T;
t = 0: T/1000 : T;
m=(-7:7);
x_2 = x2(t,wA,wB);

F_Matrix=exp(1i*w0*t'*m);
tn_uneven_s=(T-0.0000000001)*rand(1,40);
t_s_rand=tn_uneven_s+0.01*rand(1,40);
x2_uneven_s=x2(tn_uneven_s,wA,wB);

plot(t,x_2,Color=[0,0,0]);
hold on
plot(tn_uneven_s,x2_uneven_s,'*',Color=[1,0,0]);
title("Signals x(t) and uneven sampled x_s(t) (40 samples) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
legend('x(t)','x_s(t)');
hold off


F_uneven_40=exp(1i*w0*t_s_rand'*m); 
c_uneven_40=cond(F_uneven_40);
a_uneven_40=inv(F_uneven_40'*F_uneven_40)*F_uneven_40'*x2_uneven_s.';

condf_rand_f=cond(F_uneven_40);
x_vec_uneven_40=F_Matrix*a_uneven_40;

figure;
plot(t,x_2,Color=[0,0,0]);
hold on
plot(t,x_vec_uneven_40,'-.',Color=[1,0,0]);
title("Signals x(t) and x_r_e_c_5(t) (40 samples) - Time Domain");
xlabel("t [sec]");
ylabel("x(t) [V]");
legend('x(t)','x_r_e_c_5(t)');
hold off
clear all
clc
%_________________________________________________%
% Q3 - B
% Compute Cn projection coefficients vectors 
% For each function with each kernel function

T=10;
t=0:(1/500):T;
f_t=f(T,t);
g_t=g(T,t);

n_1=-20:20;
PHI = phi_func(n_1,t);
n_2=0:19;
PSI=psi_func(n_2,t',0);

phi_f=projection_coef_vec(f_t,PHI,T);
phi_g=projection_coef_vec(g_t,PHI,T);
psi_f=projection_coef_vec(f_t,PSI,T);
psi_g=projection_coef_vec(g_t,PSI,T);
%_________________________________________________%
% Q3 - C
% Reconstruction of f and g from projection coefficients vectors

f_rec_phi=PHI*phi_f.';
f_rec_psi=PSI*psi_f.';
g_rec_phi=PHI*phi_g.';
g_rec_psi=PSI*psi_g.';

figure;
plot(t,f_t,Color = [0,0,0])
hold on
plot(t,f_rec_phi,'--',Color=[1,0,0])
hold on
plot(t,f_rec_psi,'-.',Color = [0,0,1])
title("f(t) and reconstructions of f(t):")
xlabel("t [sec]")
ylabel("f(t)")
legend('f(t)','f_Φ(t)','f_Ψ(t)')
hold off

figure;
plot(t,g_t,Color = [0,0,0])
hold on
plot(t,g_rec_phi,'--',Color=[1,0,0])
hold on
plot(t,g_rec_psi,'-.',Color = [0,0,1])
title("g(t) and reconstructions of g(t):")
xlabel("t [sec]")
ylabel("g(t)")
legend('g(t)','g_Φ(t)','g_Ψ(t)')
hold off
%_________________________________________________%
%%
%_________________ Q1 Functions __________________%
function [res]=x(t,wm)
% x(t) function Q1
    res=4./(wm*pi*t.^2).*(sin(wm*t)).^2.*(cos(wm*t)).*(sin(2*wm*t)); 
end
%_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _%
function [res]=X(wm,w,ws,k)   
% X(ω) function
    tr1 = triangularPulse(-wm - ws*k,wm - ws*k,3*wm - ws*k,w);
    tr2 = triangularPulse(-3*wm - ws*k,-wm - ws*k,wm - ws*k,w);
    tr3 = triangularPulse(wm - ws*k,3*wm - ws*k,5*wm - ws*k,w);
    tr4 = triangularPulse(-5*wm - ws*k,-3*wm - ws*k,-wm - ws*k,w);
    res = (1/1i)*(tr1-tr2+tr3-tr4);
end
%__________________ Q2 Function __________________%

function [res]=x2(t,wA,wB)
% x(t) function Q2 
    res=5*cos(wA*t)-3*sin(wB*t);
end
%_________________ Q3 Functions __________________%
                    % Q3-A
function [Cn]=projection_coef_vec(vec,mat,T)  
% Compute projection coefficients vectors
    dt = T/(length(vec)-1);
    t = 0:dt:T;
    Cn = (trapz(t',vec.'.*conj(mat),1)./trapz(t',mat.*conj(mat),1));  
end
%_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _%
function [res]=psi_func(n,t,k)
% Ψ_n(t) function
   T=10;
   res = 1*(t<=(T/20*(n+0.5)-20*k+T/40)&t>T/20*(n+0.5)-20*k-T/40);
end
%_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _%
function [res]=phi_func(n,t)
% Φ_n(t) function
   T=10;
   res = exp((1i*2*pi/T)*t.'*n);
end
%_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _%
function [res]=f(T,t)
% f(t) function
    res=4*cos((4*pi/T)*t)+sin((10*pi/T)*t);
end
%_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _%
function [res]=g(T,t)
% g(t) function
    res=2*sign(sin((6*pi/T)*t))-4*sign(sin((4*pi/T)*t));
end
%_________________________________________________%