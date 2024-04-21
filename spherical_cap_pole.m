%% calculate SH expansion for a polar cap at some point above surface
% from Pollack 1973
% polar for convenience because then its just zonal

%% earth constants
G=6.6743e-20; % universal gravitational constant,  km^3 * kg^-1 * s^-2
Re=6378.3; % equatorial radius of earth, km

%% cap constants
alpha_d=1; % generating angle of cap, deg
sigma = 1e10; % areal density of surface mass distribution, kg/km^2
% 1e12 = area density of 1km depp water at 1000kg/m^3

cosalpha=cosd(alpha_d); % defined for convenience

Acap = 2.*pi.*Re.^2.*(1-cosalpha); % area of spherical cap
Mcap=sigma.*Acap; % mass of cap, kg
% alternatively, just define total mass of cap.
% for reference, 1Gt = 1e12 kg
Gt=1e12;

%% satellite coordinates

% 1 = vary colatitude of satellite
varytheta=1;

if varytheta==1 
    theta_d = 0:1:180;
    alt=500;
else
    theta_d = 90; % colatitude of satellite
    alt=(0:1:1000); % altitude of satellite
end

r = Re+alt; % distance from center of earth to satellite, km

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
%plot(altvec,Vn(1,:));hold on
%plot(altvec,Vn(2,:));
%plot(altvec,Vn(3,:));
%plot(altvec,Vn(4,:));
%plot(altvec,Vn(100,:))
%legend('n=0','n=1','n=2','n=3','n=100')
for ii=1:nmax
    plot(theta_d,Vn(ii,:)); hold on
end

%xlabel('altitude above cap, km')
%ylabel('gravitational potential, m^2/s^2')

figure(2);clf;
plot(theta_d,sum(Vn,1))

%figure(3);clf;
%plot(Vn(:,1000))