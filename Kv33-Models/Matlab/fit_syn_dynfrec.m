% Fit synapse model to depression/recovery data.
% Fit over multiple frequencies
% Least squares fit to data is calculated and returned.
% BPG 31-8-11 & 5-3-21

function f = fit_syn_dynfrec(x)

global Pv0 P1 tauf trB trH trR D tauD spt100 spt200 spt600 ndat100 ndat200 ndat600


Pv0=x(1);
trH=x(2);
trR=x(3);
trB=x(4);

% Calculate model response (100 Hz)
[n,Pv,frD,psr,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt100);
% Error
f = sqrt(sum((ndat100-(psr(1:length(ndat100))./psr(1))').^2));

% Calculate model response (200 Hz)
[n,Pv,frD,psr,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt200);
% Error
f = f + sqrt(sum((ndat200-(psr(1:length(ndat200))./psr(1))').^2));

% Calculate model response (600 Hz)
[n,Pv,frD,psr,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt600);
% Error
f = f + sqrt(sum((ndat600-(psr(1:length(ndat600))./psr(1))').^2));