%%
clc
clear all;
close all;

%-------------------------- Q1 -----------------------------------
%--------------------------- A -----------------------------------
N=30;
angle1=pi/10.25;
angle2=2*pi/5;
n=0:N-1;
s = 2*cos(angle1.*n);
v = 3*sin(angle2.*n);
x = s+v;

figure; 
hold on;
plot(n,s,'b-o');
plot(n,v,'r-o');
plot(n,x,'black-o');

title('x[n] , s[n] , v[n] - Time Domain'); 
xlabel('n'); 
ylabel('x[n] , s[n] , v[n] [V]');
legend('s[n]','v[n]','x[n]');

X_d = fft(x);
S_d = fft(s);
V_d = fft(v);

figure;
hold on;
stem(n,abs(S_d),'b'); 
stem(n,abs(V_d),'r');
stem(n,abs(X_d),'--black');
title('| X^d[k] | , | S^d[k] | , | V^d[k] | - Frequency Domain (DFT)'); 
xlabel('k'); 
ylabel('| X^d[k] | , | S^d[k] | , | V^d[k] |');
legend('|S^d[k]|','|V^d[k]|)','|X^d[k]|');

%--------------------------- B -----------------------------------
Xz_d= fft(x,45); 
N_2=0:44;
figure; 
hold on;
stem(n*(2*pi/30),abs(X_d),'r'); 
stem(N_2*(2*pi/45),abs(Xz_d),'black--');
title('Comparing | X^d[k] | and | X_z^d[k] |'); 
xlabel('freq[\theta]');
ylabel('| X^d[k] | , | X_z^d[k] | ');
legend('|X^d[k]|','|X_z^d[k]|');

%--------------------------- C -----------------------------------
N2=45;

s_N2 = 2*cos(angle1*N_2) ;
v_N2 = 3*sin(angle2*N_2) ;
x_N2 = s_N2 + v_N2 ;
X_d2 = fft(x_N2) ;

figure;                                                                     
stem(n*(2*pi/N),abs(X_d),'r');
hold on ;
title('Comparing | X^d[k] | , | X^d_2[k] |');
xlabel('freq[\theta]');
ylabel('| X^d[k] |,| X^d_2 [k] |');
stem(2*pi*(N_2/N2),abs(X_d2),'black--');
legend('| X^d[k] |','| X^d_2 [k] |');
hold off;
%--------------------------- D -----------------------------------
xz=[x,zeros(1,15)];

X_Parseval = (1/N)*(X_d*X_d') ;
x_Parseval = x*x' ;

Xz_Parseval = (1/N2)* (Xz_d*Xz_d') ;
xz_Parseval = xz*xz' ;

%--------------------------- E -----------------------------------

Teta=0:1/1000:2*pi-0.001;
k=0:2*pi/N:(2*pi/N)*(N-1); 
H_teta=1/3*(1+exp(-1j*Teta)+exp(-2j*Teta));
figure;
plot(Teta,abs(H_teta),'r');
hold on;
stem(k,abs(X_d)/45,'black--o');
title('|X^d[k]|, |H^f(\theta)| - Frequency Domain (DFT)');
xlabel('freq [\theta]')
ylabel('| X^d[k] | , | H^d[k] |')
legend('|H^f(\theta)|','|X^d[k]|');
hold off;

x_1a=[x 0 0];
h_1a=[1/3 1/3 1/3 zeros(1,29)];
H_1a=fft(h_1a);
X_1a=fft(x_1a);

Y=X_1a.*H_1a;
k_1a=0:2*pi/32:(2*pi/32)*31;

figure;
stem(k_1a,abs(Y),'r-o');
hold on;
stem(k,abs(X_d),'black--o');
title('|X^d[k]|, |Y^d[k]| - Frequency Domain (DFT)');
xlabel('freq [\theta]')
ylabel('| X^d[k] | , | Y^d[k] |')
legend('|Y^d[k]|','|X^d[k]|');
hold off;

y1=ifft(Y);
y1=y1(1:30);
figure;
plot(n,s,'b-o');
hold on;
plot(n,v,'r-o');
plot(n,x,'black-o');
plot(n,y1,'magenta-o');
title('s[n] ,v[n], x[n], y_1[n] - Time Domain');
xlabel('n')
ylabel('x[n] , s[n] , v[n] , y_1[n] [V]')
legend('s[n]','v[n]','x[n]','y_1[n]');
hold off;
%--------------------------- F -----------------------------------
H2_teta=(1+exp(-1j*Teta));
figure;
plot(Teta,abs(H2_teta),'r');
hold on;
stem(k,abs(X_d)/15,'black--o');
title('|X^d[k]|, |H^f(\theta)| - Frequency Domain (DFT)');
xlabel('freq [\theta]')
ylabel('| X^d[k] | , | H^d[k] |');
legend('|H^f_2(\theta)|','|X^d[k]|');
hold off;

x_2a=[x 0];
h_2a=[1 1 zeros(1,29)];
H_2a=fft(h_2a);
X_2a=fft(x_2a);
Y2=X_2a.*H_2a;
k_2a=0:2*pi/31:2*pi/31*30;

figure;
stem(k_2a,abs(Y2),'r-o');
hold on;
stem(k,abs(X_d),'black--o');
title('|X^d[k]|, |Y^d[k]| - Frequency Domain (DFT)');
xlabel('freq [\theta]')
ylabel('| X^d[k] | , | Y^d[k] |')
legend('|Y^d[k]|','|X^d[k]|');
hold off

y2=ifft(Y2);
y2=y2(1:30);
figure;
plot(n,s,'b-o');
hold on;
plot(n,v,'r-o');
plot(n,x,'black-o');
plot(n,y2,'magenta-o');
title('s[n] ,v[n], x[n], y_2[n] - Time Domain');
xlabel('n')
ylabel('x[n] , s[n] , v[n] , y_2[n] [V]')
legend('s[n]','v[n]','x[n]','y_2[n]');
hold off;

clc
clear all;
close all;

%%

%-------------------------- Q2 -----------------------------------
load('data2020.mat')
n=0:1:length(y)-1;
x_test=[x_test zeros(1,length(y)-length(x_test))];
X_test=fft(x_test);
Y_test=fft(y_test);
Y_z=fft(y_z);
Y=fft(y);

Y_0=Y_test-Y_z;
H_rec=X_test./Y_0;

h_rec=ifft(H_rec);

figure;
stem(n,h_rec,'LineStyle','-','Marker','none');
title('h_{rec}[n]');
xlabel('n');
ylabel('h_{rec}[n] [v]');
legend('h_{rec}[n]');

x_rec=ifft(H_rec.*(Y-Y_z));

figure;
stem(n,y,'b','-','Marker','none');
hold on;
stem(n,x_rec,'r','-','Marker','none');
title('x_{rec}[n]');
xlabel('n');
ylabel('x_{rec}[n][v]');
legend('y[n]','x_{rec}[n]');
Fs=44100;
soundsc(x_rec,Fs);
hold off;