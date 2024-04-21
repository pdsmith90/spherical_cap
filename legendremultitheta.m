function [Pn,d1Pn,d2Pn] = legendremultitheta(nmax,theta_d)
% legendre polys n=0..nmax, k=0
% and first and second derivatives
% usage [Pn,d1Pn,d2Pn]=legendre_pds(5,50)
% in:   theta_d, colatitude in degrees
%       nmax = max degree (must be greater than 1)
% out:  Pn = Pn(cos(theta)) legendre polys evaluated at cos(theta)
%       d1Pn = first derivative of legendre poly
%       d2Pn = second derivative of legendre poly

%% setup

%make sure max order >=1
if nmax<1
    error 'Max order cannot be less than 1'
end

% given a theta, in degrees
% precompute sintheta, costheta, using builtin cosd and sind
costheta=cosd(theta_d);
sintheta=sind(theta_d);

% make matrix of Pn for whatever high order we want
% matrix is 1-indexed but Pn is 0-indexed.
Pn=nan(nmax+1,length(theta_d));
d1Pn=nan(nmax+1,length(theta_d));
d2Pn=nan(nmax+1,length(theta_d));

%% low orders for seeding recurrence formulae

Pn(1,:)=1; %P0
Pn(2,:)=costheta; %P1
%Pn(3)=3/2.*cos2theta-1/2; %P2
%Pn(4)=5/2.*cos3theta-3/2.*costheta; %P3

%firstderivs
d1Pn(1,:)=0; %d1P0
d1Pn(2,:)=-sintheta; %d1P1
%d1Pn(3)=3.*sintheta.*costheta; %d1P2

%second derivatives
d2Pn(1,:)=0; % d2P0
d2Pn(2,:)=-costheta; %d2P1
%d2Pn(3)=-3.*cos2theta+3.*sin2theta; %d2P2

%% loop
for n=1:(nmax-1)
    %again, have to shift indexing
    
    Pn((n+1)+1,:)=( (2*n+1).*costheta.*Pn((n+1),:) - n.*Pn((n+1)-1,:) )./(n+1);
    
    d1Pn((n+1)+1,:)=costheta.*d1Pn((n+1),:)-(n+1).*sintheta.*Pn((n+1),:);

    d2Pn((n+1)+1,:)=d2Pn((n+1)-1,:)-(2*n+1).*(costheta.*Pn((n+1),:)+sintheta.*d1Pn((n+1),:));

end