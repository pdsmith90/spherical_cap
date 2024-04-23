%% calculate SH expansion for a polar cap at some point above surface
% from Pollack 1973
% polar for convenience because then its just zonal

%% earth constants
G=6.6743e-20; % universal gravitational constant,  km^3 * kg^-1 * s^-2
Re=6378.3; % equatorial radius of earth, km

%% cap constants
alpha_d=0.5; % generating angle of cap, deg
sigma = 1e9; % areal density of surface mass distribution, kg/km^2
% 1e12 = area density of 1km deep water at 1000kg/m^3

cosalpha=cosd(alpha_d); % defined for convenience

Acap = 2.*pi.*Re.^2.*(1-cosalpha); % area of spherical cap
Mcap=sigma.*Acap; % mass of cap, kg
% alternatively, just define total mass of cap.
% for reference, 1Gt = 1e12 kg
Gt=1e12;

%% satellite coordinates

% 1 = vary colatitude of satellite
varytheta=0;

if varytheta==1 
    theta_d = 0:1:180;
    altvec=500;
else
    theta_d = 0; % colatitude of satellite
    altvec=(0:1:1000); % altitude of satellite
end

r = Re+altvec; % distance from center of earth to satellite, km

%% SH nmax & legendre polys

% this nmax needs to be +1 since the final formula is summation over n
% which includes Pn+1 term
nmax=101;
n=(0:1:nmax)';

% legendre polys have to be calculated for generating angle of cap
% and for the colatitude of satellite
[Pncap,~,~]=legendremultitheta(nmax,alpha_d);
[Pnsat,~,~]=legendremultitheta(nmax,theta_d);

% right now do not feed vectors as angles for convenience


%% calculate V per degree
% array with dimensions = n x alt
% reminder that matlab is 1-indexed :)

% preallocate
if varytheta==1
    Vn=zeros(nmax,length(theta_d));
else
    Vn=zeros(nmax,length(altvec));
end
% seperately calculate n=0 point mass
Vn(1,:)=G.*Mcap./r;

% big element-wise
Vn(2:end,:)=-G.*Mcap./r.*(Re./r).^n(2:end-1)./ ...
    ((2.*n(2:end-1)+1).*(1-cosalpha)).* ...
    Pnsat(2:end-1,:).*(Pncap(3:end,:)-Pncap(1:end-2,:));

% put it from km^2/s^2 to m^2/s^2
Vn=Vn.*1e6;

figure(1);clf;
% for ii=1:nmax
%     plot(theta_d,Vn(ii,:)); hold on
% end
% plot(theta_d,Vn(20,:))
% plot(theta_d,Vn(40,:))
% plot(theta_d,Vn(40,:))

plot(altvec,Vn(1,:));hold on
plot(altvec,Vn(11,:));
plot(altvec,Vn(21,:));
plot(altvec,Vn(41,:));
plot(altvec,Vn(61,:));
plot(altvec,Vn(101,:));
legend('n=0','n=10','n=20','n=40','n=60','n=100','location','east')
% for ii=1:nmax
%     plot(theta_d,Vn(ii,:)); hold on
% end
title('Gravitational potential due to 10 Gt spherical cap')
xlabel('altitude above cap, km')
ylabel('gravitational potential, m^2/s^2')
% 
% figure(2);clf;
% plot(theta_d,sum(Vn,1))
% xlabel('Angle away from cap, deg')
% ylabel('Gravitational Potential due to cap')

%figure(3);clf;
%plot(Vn(:,1000))