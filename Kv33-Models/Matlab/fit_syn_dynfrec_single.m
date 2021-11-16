% Fit synapse model to depression/recovery data.
% Fit over a single frequency
% Least squares fit to data is calculated and returned.
% BPG 31-8-11 & 20-5-21

function f = fit_syn_dynfrec_single(x)

global Pv0 P1 tauf trB trH trR D tauD spt ndat ndatrec


Pv0=x(1);
trH=x(2);
trR=x(3);
trB=x(4);

% Calculate model response
[n,Pv,frD,psr,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt);
% Error
%f = sqrt(sum((ndat-(psr(1:length(ndat))./psr(1))').^2));
ndatall = [ndat' ndatrec']; % stim + recov
f = sqrt(sum((ndatall-(psr./psr(1))).^2));