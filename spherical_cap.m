function [Vn, Phi_n, Psi_n]=...
    spherical_cap(r, theta, lambda,...
    alpha,sigma,theta_prime,lambda_prime,nmax)
%% calculate SH expansion for a cap anywhere on sphere
% compute gravitational potential at a point due to a uniform density
% spherical cap that is anywhere on surface of reference sphere
% 
% in:   r (km), phi, theta (deg) ECEF coords of field point (aka satellite)
%       alpha_d (deg) generating angle of cap
%       sigma (kg/km^2) area density of cap
%       phi_prime, theta_prime (deg) ECEF location of center of cap
%       nmax = max degree
%
% out:  Vn (m^2/s^2) = gravitational potential per degree (order summed)
%       Phi_n = Kaula's power spectrum per degree
%       Psi_n = Pollack's power spectrum per degree
%
% from Pollack 1973
% eqs. 2, 5, 9, 15, 16
% variables in equations:
%
% r = geocentric radius of field point
% theta = colatitude of field point
% lambda = longitude of field point
%
% a = radius of reference sphere on which cap is located
% theta_prime = colatitude of center of spherical cap
% lambda_prime = longitude of center of spherical cap
% alpha = generating angle of spherical cap
% sigma = area density of spherical cap mass
% M = total mass of spherical cap
%   = 2*pi*a^2*(1-cos*alpha)
% 
% note: Pollack uses some shorthand variables & non-convential lat&long
%       which I don't follow:
% phi = colatitude
% theta = longitude
% lambda = cos(alpha)
% mu = cos(phi)
% 
% *** need NALF.m function ***
% *** need legendremultitheta.m function ***
%

%% constants
G=6.6743e-20; % universal gravitational constant,  km^3 * kg^-1 * s^-2
a=6378.3; % equatorial radius of earth, km

%alpha=1; % generating angle of cap, deg
cosalpha=cosd(alpha); % defined for convenience

%sigma = 1e10; % areal density of surface mass distribution, kg/km^2
% 1e12 = area density of 1km depp water at 1000kg/m^3

Acap = 2.*pi.*a.^2.*(1-cosalpha); % area of spherical cap
Mcap=sigma.*Acap; % mass of cap, kg
% alternatively, just define total mass of cap.

% for reference, 1Gt = 1e12 kg
%Gt=1e12;
%fprintf(strcat('total mass of cap: ',num2str(Mcap./Gt),' Gt\n'))

%% calculate fully normalized associated legendre polynomials
% *** calling NALF.m function ***
% NALF.m saves alpha_nm and beta_k coeffs to speed up subsequent runs
% legendre polys have to be calculated for colatitude of cap
% and for the colatitude of satellite.
% !do not feed vectors as angles for convenience!
Pnm_sat=NALF(theta,nmax,0);
Pnm_cap=NALF(theta_prime,nmax,0);

% also need to calculate degree-only legendre polys for generating angle
% *** calling legendremultitheta.m ***
% this nmax needs to be +1 since the final formula is summation over n
% which includes Pn+1 term
[Pnalpha,~,~]=legendremultitheta(nmax+1,alpha);

% remember matlab 1-indexing where degree starts at 0
% Pnm(degree+1,order+1) as lower triangular matrix

%% calculate Rnm and Snm for both cap and point, and Jnm and Knm for cap
% eqs 2, 9
n=(0:1:(nmax+1))';
m=(0:1:nmax);

Rnm_sat=Pnm_sat.*cos(m.*lambda.*pi./180);
Snm_sat=Pnm_sat.*sin(m.*lambda.*pi./180);

Rnm_cap=Pnm_cap.*cos(m.*lambda_prime.*pi./180);
Snm_cap=Pnm_cap.*sin(m.*lambda_prime.*pi./180);

% hanging out in front of the equation 9
JKcoeff=(Pnalpha(3:end)-Pnalpha(1:(end-2)))./...
    (((2.*n(2:(end-1))).^2).*(1-cosalpha));
% note this part now removes the degree-zero and 
% is thus essentially not 1-indexed
Jnm_cap=JKcoeff.*Rnm_cap(2:end,2:end);
Knm_cap=JKcoeff.*Snm_cap(2:end,2:end);

% clean up some intermediate calculated vars
clear JKcoeff RNM_cap SNM_cap cosalpha;
clear Pnm_sat Pnm_cap;

%% gravitatonal potential per-degree (sum order) (eq 5)

% preallocate so that degree-0 can be put in manually
Vn=zeros(nmax+1,1);
Vn(1)=G.*Mcap./r;

% sum over order
Jn_cap=sum(Jnm_cap,2);
Kn_cap=sum(Knm_cap,2);
Rn_sat=sum(Rnm_sat(2:end),2);
Sn_sat=sum(Snm_sat(2:end),2);

% hanging out in front
Vncoeff=-G.*Mcap./r.*(a./r).^n(2:end-1);

% potential per degree
Vn(2:end)=Vncoeff.*(Jn_cap.*Rn_sat+Kn_cap.*Sn_sat);

%% power spectra (Eq 15 & 16)
sigma_n=sqrt(Jn_cap.^2+Kn_cap.^2);

Psi_n=sigma_n./sqrt(2.*n(2:(end-1))+1);
Phi_n=sigma_n.*sqrt(2.*n(2:(end-1))+1);

end