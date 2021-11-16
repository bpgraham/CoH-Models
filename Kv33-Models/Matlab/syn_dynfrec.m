function [n,Pv,frD,psr,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt)
% function [n,Pv,frD,psr,Pr] = syn_fpfad2(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt) - synapse model
% Synapse with facilitation, dynamic (rate-changing) activity-dependent vesicle mobilisation
% and desensitisation of the postsynaptic response
% Mobilisation is from an infinite sized reserve pool
% Pv0 - base level use fraction per AP
% P1 - increment in Use per AP
% tauf - relaxation time constant of facilitation (msecs) 
% trB - time constant of background recovery (msecs)
% trH - recovery time constant instantly following a spike
% trR - time constant of relaxation of tr back to background rate
% D - fraction of desensitised receptors
% tauD - time constant of recovery from desensitisation (msecs)
% spt - vector of spike times (msecs)
% Returns size of RRVP (n), prob. of vesicle release (pv),
% fraction of desensitised receptors (frD),
% postsynaptic response (psr) and vesicle release (Pr)
% (model adapted from Fuhrmann et al, J. Neurophys. 87:140-148, 2002)
% updated to include dynamic changes in recovery rate (BPG 5-3-21)
% BPG 24-6-02 & 5-3-21

tstep = 100;    % number of time steps to integrate between spikes

%Generate spikes
spcnt = length(spt);	% number of spikes
isi = zeros(1,spcnt-1);
for i=2:spcnt
   isi(i-1) = spt(i) - spt(i-1);
end;

% Generate PRs (probability of release)

ffac = 1-exp(-isi./tauf);	% facilitation
%frec = 1-exp(-isi./taur);	% vesicle depletion
frecD = 1-exp(-isi./tauD);	% desensitisation

n = zeros(1,spcnt);	% avg. no. of available vesicles (RRVP)
n(1) = 1;
taur = trB;
Pv = zeros(1,spcnt);
Pv(1) = Pv0;
frD = zeros(1,spcnt);	% fraction desensitised receptors

for i=2:spcnt
   % release probability with facilitation
   Pvp = Pv(i-1)+P1*(1-Pv(i-1));
   Pv(i) = Pvp + ffac(i-1)*(Pv0-Pvp);
   % RRVP following release at previous spike
   np = (1-Pv(i-1))*n(i-1);
   % integrate vesicle recovery between spikes
   dt = isi(i-1) / tstep; % time step for integration
   n(i) = np;
   taur = trH;  % high recovery rate
   for k=1:tstep-1
       n(i) = (dt/taur)*(1-n(i)) + n(i);
       taur = (dt/trR)*(trB-taur) + taur;
   end;
%   frDp = frD(i-1) + (D*(1-frD(i-1)));
   frDp = frD(i-1) + (D*Pv(i-1)*n(i-1)*(1-frD(i-1)));
   frD(i) = frDp - frecD(i-1)*frDp;
end;

Pr = n.*Pv;			% prob. release
psr = Pr.*(1-frD);	% prob. release x fraction not desensitised

