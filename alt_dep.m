%% earth constants
G=6.6743e-20; % universal gravitational constant,  km^3 * kg^-1 * s^-2
mu=398600.440; % gravitational constant for earth, km^3 * s^-2
Re=6378.3; % equatorial radius of earth, km
sid_day_h=24.*(365.24./366.24); % sidereal day, hr
OMdot_e=2.*pi./(sid_day_h.*3600); % rotation rate of earth, rad/s
J2=0.0010826; % J2 oblateness coeff, dimensionless

%% altitude dependence of gravitational attraction
% estimate mass of anomalies as point masses on Earth's surface

altvec1=(0:1:1000)';

altGRACE=470;
altGFO=506;
altCHAMP=450;
altGOCE=255;
altSWARM=460;
altCOSMIC=520;
altSPIRE=400; %??? or is it 500km?

% mass of different anomalies
anomalymass=[1e14,1e13,1e12,1e11];
anomalynames={'100 Gt','10 Gt','1 Gt = 1e12 kg','0.1 Gt'};

anomalyaccel=G.*anomalymass./(altvec1.^2).*1e3;

figure(1);clf
for ii=1:length(anomalymass)
    semilogy(altvec1,anomalyaccel(:,ii));hold on
end
legend(anomalynames,'AutoUpdate','off');
grid on
xlabel('altitude, km');
ylabel('acceleration due to anomaly mass, m/s^2')

ynameplace=ylim;

xline(altGOCE)
text(altGOCE+10,ynameplace(2)./1000,'GOCE')

xline(altCOSMIC)
text(altCOSMIC+10,ynameplace(2)./1000,'COSMIC-2')

xline(altGRACE)
text(altGRACE+10,ynameplace(2)./100,'GRACE')

title('Acceleration due to mass anomalies at orbit overflight')
yline(1e-11,'r','LineWidth',2)
text(10,5e-12,'GRACE accelerometer precision','Color','r')

%% wavelength of gravity anomaly vs height
% for convenience, pick zonal coeffs, n=something, m=0
% plot side view cross section of gravity anomaly
% due to individual unity Cn,0 coefficients

% max degree to compute
nmax=101;

% pick a range of altitudes
altvec2=(0:1:1000)';

% calculate legendre function
% this is easy non - associated legendre polynomials

% range of colatitude, in degrees
theta_d=0:0.1:180;
% precompute sintheta, costheta, using builtin cosd and sind
costheta=cosd(theta_d);
sintheta=sind(theta_d);

% make matrix of Pn & first derivs for whatever high order we want
% matrix is 1-indexed but Pn is 0-indexed.
Pn=nan(nmax+1,length(theta_d));
d1Pn=nan(nmax+1,length(theta_d));

% low orders for seeding recurrence formulae
Pn(1,:)=1; %P0
Pn(2,:)=costheta; %P1
%Pn(3)=3/2.*cos2theta-1/2; %P2
%Pn(4)=5/2.*cos3theta-3/2.*costheta; %P3
%firstderivs
d1Pn(1,:)=0; %d1P0
d1Pn(2,:)=-sintheta; %d1P1
%d1Pn(3)=3.*sintheta.*costheta; %d1P2

% recurrence relation
for n=1:(nmax-1)
    %have to shift indexing
    Pn((n+1)+1,:)=( (2*n+1).*costheta.*Pn((n+1),:) - n.*Pn((n+1)-1,:) )./(n+1);
    d1Pn((n+1)+1,:)=costheta.*d1Pn((n+1),:)-(n+1).*sintheta.*Pn((n+1),:);
end

% set n vector for element-wise multiplications
n=0:1:nmax;


% calculate everything that hangs out in front of legendre fn
% potential
Varr=mu./(Re+altvec2).*(Re./(Re+altvec2)).^n;
% force in radial
Frarr=-mu.*(Re+altvec2).^(-2).*(n+1).*(Re./(Re+altvec2)).^n;
% force in latitude direction
Ftarr=-mu.*(Re+altvec2).^(-2).*(Re./(Re+altvec2)).^n;

F0=(abs(Frarr(:,1)).^2+abs(Ftarr(:,1)).^2).^0.5;

figure(2);clf
% plot a selection of values for n
colors={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE'};
npickv=[10,20,40,60];
for ii=1:4
    npick=npickv(ii);
    % plot a range of altitudes
    for altpick=100:100:400
        % accel as sum of Fr & Ftheta (Flambda=0)
        plot(0.25./F0(1).*...
            ((abs(Frarr(altpick+1,npick+1).*Pn(npick+1,:)')).^2+...
            (abs(Frarr(altpick+1,npick+1).*d1Pn(npick+1,:)')).^2).^0.5+...
            altpick,...
            90-theta_d,'Color',colors{ii})
        hold on
    end

end

title('Normalized acceleration due to different degree SH')
xlabel('Altitude above mass, km')
ylabel('Latitude, degrees')
grid on
ylim([-90,90])
xlim([100,420])
legend('Degree 10','','','','','','',...
    'Degree 20','','','','','','',...
    'Degree 40','','','','','','',...
    'Degree 60','Location','East')


