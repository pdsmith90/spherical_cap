function [P, vargout] = NALF(theta_d,nmax,wanterrorcheck,vargin)
% Geodetic Normalized Associated Legendre Polys for a given max
% degree and value of colatitude in degrees
%
% inputs:  
%    thetad   =   colatitude [degrees] (can accept vectors)
%    nmax     =   max degree N
%    wanterrorcheck = error checking or no (spits into command window)
%    vargin   =   optional arg to say how many derivs to gen
%
% outputs: 
%    P        =   matrix of 1-indexed poly values for given theta
%                 P(1,1)=P_0^0 => everything shifted by 1
%    vargout  =   optional output cell array of derivs reqested
%                 s.t. vargout{n} is the nth deriv of P {curly bracket}
%
% example use 1
%    for theta= 0.1 degrees, Nmax = degree 1800, with 2nd derivative
%    [P, dP] = NALF(0.1,1800,2)
%    P is a matrix
%    dP{1} is 1st derivative & dP{2} is 2nd derivative (curly brackets!)
%
% example use 2
%    for theta= 0.1 degrees, Nmax = degree 2100
%    P = NALF(0.1,2100)
%
% notes:
% 
% 1. matlab uses 16 decimal precision; for 32 bit or more, 
%    wrap the vpa function (from Symbolic Math) around all calculations.
%    define variable digits to set the number of significant digits in vpa
%
% 2. to run subroutines independantly, copy and save them into seperate file
%
% 3. the alpha, beta, d functions save outputs into files that can be
%    opened again later to reduce computing time.
%
% 4. If nmax>1800-2000 ish, then be sure to be around the equator
%
% updates
% Feb 4 2022: can now accept vector of thetas, represented by the
%    third dimension in output P
% Feb 25 2022: added some parallel processing parfors to speed things up
%    but then that doesnt work so nevermind


%% setup

% check to make sure we are dealing with n>1
if nmax<=1
    error 'nmax must be greater than one'
end

% retrieve / calculate alpha_nk values: they are independant of theta!
alpha_nk=alpha_nk_fun(nmax);

% retrieve / calculate beta_k values: they are independant of theta!
beta_k=beta_k_fun(nmax);

% convert theta to rads 
theta=theta_d.*pi./180;

% precompute cos(theta) and sin(theta) because we use it a lot
% mess with adding dimensions for the 3d
costheta=zeros(1,1,length(theta));
sintheta=zeros(1,1,length(theta));
costheta(1,1,:)=cosd(theta_d);
sintheta(1,1,:)=sind(theta_d);

% preallocate P matrix
P=zeros(nmax+1,nmax+1,length(theta)); %always shifted by one

% seed the first two diagonals, P_0^0, P_1^1
P(1,1,:)=1;
P(2,2,:)=sqrt(3).*sintheta;

%% loop 1: diagonal

for jj = 3:(nmax+1)
    P(jj,jj,:)=beta_k(jj).*sintheta.*P(jj-1,jj-1,:);
end

%% loop 2: diagonal incremented by n+1

for jj=1:nmax
    P(jj+1,jj,:)=costheta.*P(jj,jj,:)./alpha_nk(jj+1,jj);
end

%% loop 3: for each k, compute the increasing values of n

% nested loops
for ii=2:nmax
    for jj=1:ii
        P(ii+1,jj,:)= ...
            (costheta.*P(ii,jj,:)-alpha_nk(ii,jj).*P(ii-1,jj,:)) ...
            ./alpha_nk(ii+1,jj);
    end
end

%% error checking sums
if wanterrorcheck==1
    [sanity,~,error2]=errorcheck(P(:,:,1));
    if sanity==1
        disp('P sum error checking success')
    else
        % something went wrong, print total sum to see how bad it is?
        disp('P sum error error checking failed')
        fprintf('total sum error value of : %.2f\n', error2);
    end
end

%% optional derivative output

% check to see if it has been requested as the optional last input
if nargin>3
    % how many derivs requested?
    nderiv=vargin(1);
    % make sure we have a positive number here
    if nderiv<1
        return
    end
    vargout=cell(1,nderiv);
    % first deriv
    vargout{1}=ALF_deriv(P);
    % higer derivs
    if nderiv>=2
        for nn=2:nderiv
            %recursively call derivative solver
            vargout{nn}=ALF_deriv(vargout{nn-1});
        end
    end
end

end

%% subroutines

function [alpha_nk]=alpha_nk_fun(nmax)
%calculate / retrieve alpha values: they are independant of theta!
%input: maximum order
%output: matrix of values for alpha_nk: rows=n, cols=k
%also saves var alpha_nk to file called alpha_nk.mat

% first check if the variable exists in the path
if exist('./alpha_nk.mat','file') == 0
    alpha_nk=0; % first value is zero
else
    load ./alpha_nk.mat alpha_nk% exists in path, so load it up
end

% see if alpha goes to a high enough order
alphamax=length(alpha_nk); % should always be square so we just use length
% big wrapper condition
if alphamax<nmax+1
    % extend alpha matrix to preallocate to correct size
    alpha_nk(nmax+1,nmax+1)=0; % tricky matlab autofills zeros
    % nested loops
    % (too lazy to convert back and forth from single to double indexing)
    for ii=(alphamax+1):(nmax+1)
        n=ii-1; %index shift because you know
        for jj=1:ii
            k=jj-1; %another index shift
            % the equation
            alpha_nk(ii,jj)=sqrt((n+k)./(2.*n-1).*(n-k)./(2.*n+1));
        end
    end
    % save new finalized workspace alpha var to file
    save ./alpha_nk.mat alpha_nk
end

end

function [beta_k]=beta_k_fun(nmax)
%calculate / retrieve beta values: they are independant of theta!
%essentially the same thing as alpha_nk
%input: nmax order
%output beta_k vector of beta values for every k
% also writes var beta_k to file beta_k.mat

% first check if the variable exists in the path or workspace
if exist('./beta_k.mat','file') == 0
    beta_k=inf; % first value is infinite (beta_0)
else
    load ./beta_k.mat beta_k % exists in path, so load it up
end

% see if beta goes to a high enough order
betamax=length(beta_k); % should always be square so we just use length
% big wrapper condition
if betamax<nmax+1
    % extend beta matrix to preallocate to correct size
    beta_k(nmax+1)=0; % tricky matlab autofills zeros
    % no need for nested loops this time
    for jj=(betamax+1):(nmax+1)
        k=jj-1; %index shift because you know
        % the equation
        beta_k(jj)=sqrt(1+0.5./k);
    end
    % save new finalized workspace alpha var to file
    save ./beta_k.mat beta_k
end

end

function [d_nk]=d_nk_fun(nmax)
%calculate / retrieve d values: they are independant of theta!
%essentially the same thing as alpha_nk
%input: maximum order
%output: matrix of values for alpha_nk: rows=n, cols=k
%also saves var d_nk to file called d_nk.mat

% first check if the variable exists in the path
if exist('./d_nk.mat','file') == 0
    d_nk=0; % first value is zero
else
    load ./d_nk.mat d_nk% exists in path, so load it up
end

% see if alpha goes to a high enough order
dmax=length(d_nk); % should always be square so we just use length
% big wrapper condition
if dmax<nmax+1
    % extend alpha matrix to preallocate to correct size
    d_nk(nmax+1,nmax+1)=0; % tricky matlab autofills zeros
    % nested loops
    % (too lazy to convert back and forth from single to double indexing)
    for ii=(dmax+1):(nmax+1)
        n=ii-1; %index shift because you know
        for jj=1:ii
            k=jj-1; %another index shift
            % the equation
            d_nk(ii,jj)=sqrt((n+k+1).*(n-k))./2;
        end
    end
    % save new finalized workspace alpha var to file
    save ./d_nk.mat d_nk
end

end

function [sanity,errorlevel1,errorlevel2]=errorcheck(P)
%error checking against known summations
%input: normalized ALF values for given theta
%output: sanity = overall check
%        errorlevel1 = per n error value
%        errorlevel2 = total error value

nmax=length(P)-1;

errorlevel1=zeros(1,nmax+1);
errortol=1e-8;


% for any n,
% we have two conditions. & can calculate the amount of error

% for all n together, we can calculate the total summation,
% as well as a total error level

sum2=0;
%overall error check
sanity=1;

for ii=1:(nmax+1)

    % start sum over for every n
    sum1=0;
    % define n as offset by one
    n=ii-1;

    for jj=1:ii

        % the summations of interest
        sum1=sum1+P(ii,jj)^2;
    
    end

    % second sum should add up to 2n+1
    if sum1/(2*n+1)<1-errortol || sum1/(2*n+1)>1+errortol
        sanity=0;
    end

    % calculate error level
    errorlevel1(ii)=(sum1-(2*n+1))/(2*n+1);
    
    %keep adding up overall sum
    sum2=sum2+sum1;

end

% overall sum error
if sum2*(nmax+1)^-2>=1+errortol || sum2*(nmax+1)^-2<=1-errortol
    sanity=0;
end

%calculate error level
errorlevel2=(sum2-(nmax+1).^2).*(nmax+1)^-2;

end

function [dP]= ALF_deriv(P)
%recursive derivative calculator
% input:    P  = normalized ALF 
%                (or any order derivative)
% output:   dP = derivative of normalized ALF wrt co-latitude
%                (or any order derivative +1)

nmax=size(P,1)-1;

% calculate/retrieve coeff d_n^k = independant of theta!
d_nk=d_nk_fun(nmax);

%preallocate, this also satisfies k=0,n=0: dP=0
dP=zeros(size(P));

% calculate dP for k=0, all n=1..nmax
for ii=2:(nmax+1)
    % always offset by 1
    dP(ii,1,:)=-sqrt(2).*d_nk(ii,1).*P(ii,2,:);
end

% hard code dP for k=1,n=1
dP(2,2,:)=P(2,1,:);

% here on needs nmax>=2
if nmax<2
    return
end

% calculate dP for k=2..N, n=k diagonals
for jj=3:(nmax+1)
    dP(jj,jj,:)=d_nk(jj,jj-1).*P(jj,jj-1,:);
end

% calculate dP for k=1, n=2..nmax
for ii=3:(nmax+1)
    dP(ii,2,:)=sqrt(2).*d_nk(ii,1).*P(ii,1,:)-d_nk(ii,2).*P(ii,3,:);
end

% calculate for all k>1, n=(k+1)..nmax
for jj=3:(nmax) % don't have to do last element because it is only the diagonal
    for ii=(jj+1):(nmax+1)
        dP(ii,jj,:)=d_nk(ii,jj-1).*P(ii,jj-1,:)-d_nk(ii,jj).*P(ii,jj+1,:);
    end
end


end
