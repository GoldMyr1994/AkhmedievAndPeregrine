% clc
% clear all
% close all

%% 
% % MATERIAL PROPERTIES
% beta2=-1;
% gamma=1;

t = -10:.2:11;
z = -10:.1:10;

[T, Z] = meshgrid(t,z);

%%
% FREQUENCY COORDINATE
indfreq=-length(t)/2:1:length(t)/2-1;
omega=(pi./10).*indfreq;
[O,Z1] = meshgrid(omega,z);

gamma=1; % we can change this parameter
a=(1/2)*(1-(gamma*2)/4);
b=sqrt(8*a*(1-2*a));

N=(1-4*a)*cosh(b*Z)+sqrt(2*a)*cos(gamma*T)+1j*b*sinh(b*Z);
D=sqrt(2*a)*cos(gamma*T)-cosh(b*Z);

PER = (1-4*(1+2*1j*Z)./(1+4*T.^2+4*Z.^2)).*exp(1j*Z);
AKM=(N./D).*exp(1j*Z);

AKMF = zeros(size(AKM));
for k = 1:size(AKMF,2)
    AKMF(:,k) = (abs(t(2)-t(1)))*fftshift(fft(AKM(:,k)));
end

figure
mesh(T,Z,abs(AKM))
xlabel('t')
ylabel('z')
zlabel('|F|')

figure
mesh(O,Z,abs(AKMF))
xlabel('Omega')
ylabel('z')
zlabel('|FF|')

% figure
% mesh(T,Z,abs(PER))
% xlabel('t')
% ylabel('z')
% zlabel('|F|')



