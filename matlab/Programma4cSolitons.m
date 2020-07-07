%% PROPAGATION METHOD FOR THE SOLUTION OF NLSE
% $i\partial{F(z,t)}{\partial z}-\frac{beta^{\prime\prime}}{2}\frac{\partial^2F(z,t)}{\partial t^2}+\gamma\vert F(z,t)\vert^2F(z,t)$

%% 
%

clc
clear all
close all

%% 
% MATERIAL PROPERTIES
beta2=-1;
gamma=1;

%%
% TEMPORAL COORDINATE
t0=-20;
t1=20;
nt=500;
t=linspace(t0,t1,nt);
deltat=(t1-t0)/(nt-1);

%%
% SPATIAL COORDINATE
z0=0;
z1=10;
nz=501;
z=linspace(z0,z1,nz);
deltaz=(z1-z0)/(nz-1);
h=deltaz;

%%
% FREQUENCY COORDINATE
indfreq=-nt/2:1:nt/2-1;
omega=(pi./t1).*indfreq;

%%
% INPUT ENVELOPE
% tau0=1;
% amp=1;
% FIN=amp*exp(-t.^2/(2*tau0^2)); % unchirped gaussian pulse
% c=-1;
% FIN=amp*exp(-(1+1j*c)/2*t.^2/tau0.^2);  % chirped gaussian pulse

% c=-2;
% tau0=1;
% FIN=sech(t/tau0).*exp(-1j*c*t.^2/(2*tau0^2)); % hyperbolic chirped sech pulses

% c=0;
% m=1;
% tau0=1;
% amp=1;
% FIN=amp*exp(-(1+1j*c)/2*(t/tau0).^(2*m));  % chirped super gaussian pulse
% FING=amp*exp(-(1+1j*c)/2*(t/tau0).^(2));  % compare with gaussian

% three unchirped gaussian pulses
% tau0=1;
% amp=1;
% 
% tshift1=-10;
% FIN1 = amp*exp(-(t+tshift1).^2/(2*tau0^2));
% 
% tshift2=0;
% FIN2 = amp*exp(-(t+tshift2).^2/(2*tau0^2));
% 
% tshift3=10;
% FIN3 = amp*exp(-(t+tshift3).^2/(2*tau0^2));
% 
% FIN = FIN1+FIN2+FIN3;

% FIN1 = sech(t)*exp(1j*z0/2); % FUNDAMENTAL SOLITON PULSE
% FIN2 = sech(t-10)*exp(1j*z0/2)*0;  % FUNDAMENTAL SOLITON PULSE
% FIN3 = sech(t+10)*exp(1j*z0/2)*0; %  FUNDAMENTAL SOLITON PULSE
% FIN = FIN1+FIN2+FIN3;

% FIN = (4*(cosh(3*t)+3*exp(4*1j*z0).*cosh(t))./(cosh(4*t)+4*cosh(2*t)+3*cos(4*z0)))*exp(1j*z0/2); % HIGHER ORDER SOLITON PULSE
% sum(abs(FIN) - abs(2*sech(t)))


% FIN = tanh(t)*exp(z0); % DARK SOLITON

% %solitons interaction
% q0=3.5;r=1.1;theta=0;
% FIN = sech(t-q0).*exp(1j*z0/2)+r*sech(r*(t+q0)).*exp(1j*z0/2)*exp(1j*theta);

%solitons collision
q0=10;dw=-2;
FIN = sech(t-q0).*exp(1j*dw*t)+sech(t+q0).*exp(-1j*dw*t);

% % input modulus (time)
% figure 
% plot(t,abs(FIN))
% xlabel('t')
% ylabel('|F|')
% set(findall(gcf,'type','text'),'FontSize',30)
% set(gca,'FontSize',30)
% grid()

% % input phase (time)
% figure 
% plot(t,unwrap(angle(FIN)))
% xlabel('t')
% ylabel('angle')
% set(findall(gcf,'type','text'),'FontSize',30)
% set(gca,'FontSize',30)
% grid()

FFIN = deltat*fftshift(fft(FIN));
% % FFING = deltat*fftshift(fft(FING));
% % input modulus (frequency)
% figure
% plot(omega,abs(FFIN))
% xlabel('\omega')
% ylabel('|Fourier Transform F|')
% set(findall(gcf,'type','text'),'FontSize',30)
% set(gca,'FontSize',30)
% grid()
% 
% % input phase (frequency)
% figure 
% plot(omega,unwrap(angle(FFIN)))
% xlabel('omega')
% ylabel('angle(Fourier Transform F)')
% set(findall(gcf,'type','text'),'FontSize',30)
% set(gca,'FontSize',30)

%%
%%%%%%%%%%%%%%%%%%%%%%%
% CORE OF THE PROGRAM %
%%%%%%%%%%%%%%%%%%%%%%%

F = zeros(ceil(nt),floor(nz));
FF = zeros(ceil(nt),floor(nz));

F(:,1)=FIN;
FF(:,1)=FFIN;

q = FIN;

for loop_step=2:1:nz
    
    % LINEAR DISPOERSIVE STEP
    qs=deltat*fftshift(fft(q));
    qs_old=qs;
    
    prop=beta2/2*omega.^2;
    fact=1i*prop*h;
    qs=qs_old.*exp(fact); % calculation of the propagation in the frequency domain
    q=(1/deltat)*ifft(ifftshift(qs)); % coming back in the time domain
    
    % NON LINEAR CHI3 EFFECT
    q_old=q;
    q=q_old.*exp(1i*gamma*abs(q_old).^2*h);
    
    % SAVE DATA EVERY nz
    F(:,loop_step) = q;
    FF(:,loop_step) = deltat*fftshift(fft(q));
end

save dati

figure(1)
mesh(t,z,abs(F)')
xlabel('t')
ylabel('z')
zlabel('|F|')
set(findall(gcf,'type','text'),'FontSize',30)
set(gca,'FontSize',30)

% figure(2)
% mesh(t,z,unwrap(angle(F)'))
% xlabel('t')
% ylabel('z')
% zlabel('ange(F)')
% 
% figure
% plot(t,abs(F(:,1)),t,abs(F(:,end)))
% xlabel('t')
% ylabel('|F|')
% legend('IN','OUT')
% set(findall(gcf,'type','text'),'FontSize',30)
% set(gca,'FontSize',30)
% grid()

% figure()
% plot(z,unwrap(angle(F(nt/2,:)))) % follow phase in the peak along z
% xlabel('z')
% ylabel('angle F versus z at t=0')
% legend('phase of soliton')
% set(findall(gcf,'type','text'),'FontSize',30)
% set(gca,'FontSize',30)
% grid()

% 
% figure
% plot(t,unwrap(angle(F(:,1))),t,unwrap(angle(F(:,end))))
% xlabel('t')
% ylabel('angle(F)')
% legend('IN','OUT')
% grid()
% 
% figure
% deomegaIN = -gradient(unwrap(angle(F(:,end))),t);
% deomegaOUT = -gradient(unwrap(angle(F(:,1))),t);
% plot(t,deomegaIN,t,deomegaOUT)
% xlabel('t')
% ylabel('\delta omega')
% legend('IN','OUT')
% grid()

%
figure(7)
mesh(omega,z,abs(FF)')
xlabel('omega')
ylabel('z')
zlabel('|FFT F|')
set(findall(gcf,'type','text'),'FontSize',30)
set(gca,'FontSize',30)
grid
% 
% % 
% figure(8)
% plot(omega,abs(FF(:,1)),omega,abs(FF(:,end)))
% xlabel('omega')
% ylabel('|Fourier transform F|')
% legend('IN','OUT')
% set(findall(gcf,'type','text'),'FontSize',30)
% set(gca,'FontSize',30)
% grid

