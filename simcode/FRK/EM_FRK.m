    
% Author: Pulong Ma <mpulong@gmail.com>
% Date:  May 30th, 2015
% Last Modified by: Pulong Ma
% Last Modified Date: May 30th, 2015
% Last Modified time: 17:14:58


% Purpose: 
% EM algorithm to find MLEs of K and sig2xi (fine-scale variance) and 
% this code is at least 5 times faster than the code from Fixed Rank Kring tutorial

function [K_em, sig2xi_em, T, diff]=EM_FRK(S,z,sig2eps,V_eps,V_xi,maxit,avgtol)

%% input:
% V_eps is the diagonal covariance matrix for measurement error 
% V_xi is the diagonal covariance matrix for the fine-scale variation
%     (default: V_xi=eye(size(S,1));  )
%% output:
% K_em is the covariance matrix of eta in small-scale variation
% sig2xi_em is the marginal variance for fine-scale variation 

n=size(S,1);
r = size(S,2);

% default values for the last two parameters
if nargin<7, avgtol=1e-5; end
if nargin<6, maxit=300; end
if nargin<5
  V_xi=sparse(1:n,1:n,1); 
end
if nargin<4
  V_eps=sparse(1:n,1:n,1); 
end
if nargin<3
    sig2eps=0;
end

diagV_eps = diag(V_eps);
diagV_xi=diag(V_xi);

% initial values
varest=var(z,1);
K_old=.9*varest*eye(size(S,2));
sig2xi=.1*varest;   % 
t=1;
done=0;
    
while done==0,
   
    % update help terms
    diagDinv=(sig2xi(t)*diagV_xi+ sig2eps*diagV_eps).^(-1);
    DInv = spdiags(diagDinv, 0, n, n);
    SDS = S'*DInv*S;
    temp=K_old\speye(r) + SDS;
    SigInvZ = temp\(S'*(DInv*z));

    % update K
    muEta = K_old*(S'*(DInv*z)) - K_old*(SDS*SigInvZ);
    SigEta = K_old - K_old*SDS*K_old' + K_old*SDS*(temp\SDS)*K_old';
    K_new = SigEta + muEta*muEta';

    % update sigma_xi 
    muXi = sig2xi(t)*(DInv*z - DInv*(S*SigInvZ));
    trSigInv = trace(((S'*DInv)*DInv*S)/temp) - trace(DInv);
    sig2xi(t+1) = sig2xi(t) + (sig2xi(t))^2*trSigInv/n + muXi'*muXi/n;

    % check for convergence
    diff=sum(sum((K_new-K_old).^2,1),2)+(sig2xi(t+1)-sig2xi(t))^2;
    if diff<min(avgtol*r^2, 1), done=1; end
    if t>maxit, 
        done=1; 
        %disp(strcat('Algorithm did not converge after ', num2str(maxit),' iterations')); 
    end
    
    %disp(strcat('t=', ' ', num2str(t)))
    t=t+1;
    K_old=K_new;

end

K_em=K_new;
sig2xi_em=sig2xi(t);
T=t;