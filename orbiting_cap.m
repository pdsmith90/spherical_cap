%% use Pollack 1973, potential felt by satellite due to spherical cap
% arbitrary lat,long of cap
% arbitrary radius,lat,long of object (satellite)
%
% objective: calculate r,theta,lambda of satellite as function of time
% for a specific orbit, using 2-body + J2, then calculate the
% gravitational potential experienced by the satellite due to an
% anomalous spherical cap of mass with a specific area density and
% angular size. Show how the potential varies through its orbit.
% 
% *** needs functions spherical_cap.m NALF.m, legendremultitheta.m ***
%

%% earth constants
G=6.6743e-20; % universal gravitational constant,  km^3 * kg^-1 * s^-2
mu=398600.440; % gravitational constant for earth, km^3 * s^-2
Re=6378.3; % equatorial radius of earth, km
sid_day_h=24.*(365.24./366.24); % sidereal day, hr
OMdot_e=2.*pi./(sid_day_h.*3600); % rotation rate of earth, rad/s
J2=0.0010826; % J2 oblateness coeff, dimensionless

%% reference altitudes of constellations
altGRACE=470;
altGFO=506;
altCHAMP=450;
altGOCE=255;
altSWARM=460;
altCOSMIC=520;
altSPIRE=400; %??? or is it 500km?

%% generate orbit
% pick an orbit, eg. GRACE-FO
alt=500; %altitude, km
i=89; % inclination, deg

% orbit properties
% assume circular orbit in 2-body motion
a=Re+alt; % semi-major axis, km
p=a; % semi-latus rectum = radius for circle
v=sqrt(mu./a); % velocity is constant for circle, km/s
T= 2.*pi.*sqrt(a.^3./mu); % orbital period, s
n=2.*pi./T; % mean motion, rad/s

% secular precession of ascending node due to J2
% assume this is an even rate throughout the orbit period
OMdot_s=-1.5.*n.*J2.*(Re./p).^2.*cosd(i); % rad/s

%% calculate time series of lat, long for satellite orbit
% can take some shortcuts here since we're just looking at circular orbits
% and ignoring everything but secular J2 perturbation meaning just 2-body 
% arbitrary start so set it as an ascending equator crossing at greenwich

dt=10; % time step, seconds
numorbits=24*3600/T;

% r= constant
% theta= constant = fraction of orbit by period / (# of dts)
% since keplers equation reduces to M = E
timevec=0:dt:(T*numorbits);
thvec=(timevec./T-floor(timevec./T)).*2.*pi;

thetasat=zeros(length(timevec),1);
lambdasat=zeros(length(timevec),1);

parfor ii=1:length(timevec)
    rpqw=[a.*cos(thvec(ii));a.*sin(thvec(ii));0];

    % orbit plane rotated into ECI
    % is rotation about inclination
    % ECI into ECEF is rotation about x3 by ERA
    % simplified ERA = ( OMdot_s + OMdot_e ) * time 
    % r_ecef=R3(-OM)*R1(-i)*r_pqw
    
    OM=(OMdot_s+OMdot_e).*timevec(ii);

    recef=[cos(-OM),sin(-OM),0;-sin(-OM),cos(-OM),0;0,0,1]*...
        [1,0,0;0,cosd(-i),sind(-i);0,-sind(-i),cosd(-i)]*...
        rpqw;

    thetasat(ii)=atan2(sqrt(recef(1).^2+recef(2).^2),recef(3));

    lambdasat(ii)=atan2(recef(2),recef(1));

end

lambdasat(lambdasat<0)=lambdasat(lambdasat<0)+2.*pi;

thetasat=thetasat.*180./pi;
lambdasat=lambdasat.*180./pi;


%% calculate potential due to cap at every time step
% this part takes a while to run

thetacap=90-25;
lambdacap=67;
alphacap=0.5;
sigmacap=1e9;

Acap = 2.*pi.*Re.^2.*(1-cosd(alphacap)); % area of spherical cap
Mcap=sigmacap.*Acap; % mass of cap, kg
Gt=1e12; % kg tp Gt

nmax=120;
Vn=zeros(nmax+1,length(timevec));
Psi_n=zeros(nmax,length(timevec));
Phi_n=zeros(nmax,length(timevec));
grn=zeros(nmax+1,length(timevec));
gtn=zeros(nmax+1,length(timevec));
gln=zeros(nmax+1,length(timevec));
for ii=1:length(timevec)

    [Vn(:,ii), Phi_n(:,ii), Psi_n(:,ii),...
        grn(:,ii), gtn(:,ii), gln(:,ii)]=...
        spherical_cap(a, thetasat(ii), lambdasat(ii),...
        alphacap,sigmacap,thetacap,lambdacap,nmax);
end

%% calculate acceleration due to cap at every time step
% d/dr (V)
%n=(0:1:nmax)';

% radial direction is easy b/c direct function of potential
%grn=-(n+1)./a.*Vn;

gn=grn+gtn+gln;



%% visualize

figure(1);clf;
subplot(3,1,1)
plot(timevec./T,sum(gn,1),'.');
ylabel('acceleration experienced')
hold on
%[pks,locs]=findpeaks(abs(sum(gn,1)),timevec,'MinPeakHeight',1e-6);
%plot(locs./T,pks,'ro')
subplot(3,1,2)
plot(timevec./T,90-thetasat,'.');
ylabel('latitude')
subplot(3,1,3)
plot(timevec./T,lambdasat,'.');
ylabel('longitude')
xlabel('orbit rev')


figure(2); clf
scatter(lambdasat,90-thetasat,1,sum(Vn,1));
colorbar
title(strcat('potential of ',...
    num2str(alphacap),' degree cap, ',...
    num2str(Mcap./Gt),' Gt'))
hold on
plot(lambdacap,90-thetacap,'*')


figure(3); clf
scatter(lambdasat,90-thetasat,1,sum(gn,1));
colorbar
title(strcat('accel of ',...
    num2str(alphacap),' degree cap, ',...
    num2str(Mcap./Gt),' Gt'))
hold on
plot(lambdacap,90-thetacap,'*')


figure(4);clf

[pxx,f]=plomb(sum(gn,1),timevec,dt,10,'normalized');
plot(f*T,pxx);hold on
xlabel('# periods')
xlim([0,timevec(end)./T])
xticks(0:1:floor(timevec(end)./T))
grid on
title('normalized power spectrum')

figure(5);clf

Y=fft(sum(gn,1));
Fs=1./dt;
L=length(timevec);
ff=Fs/L*(0:(L/2));

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(ff*T,P1);
xlim([0,15])

%% todo

% generate time series of inertial position of satellite due to 2-body

% rotate inertial position of satellite over time into earth fixed
% earth rotation angle is: ( OMdot_s + OMdot_e ) * time

% pick a cap

% plot groundtrack + cap + rectangular projection. (GMT)

% for each time step, get Vn Phi_n & Psi_n

% visualize this data somehow
% first, plot total potential experienced as a function of time
% alongside orbital elements (specifically ecef angles)
%

% convert potential to accel by doing derivative per sh degree


