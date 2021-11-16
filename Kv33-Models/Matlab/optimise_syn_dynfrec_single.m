% Optimise synapse model against experimental time series
% Optimise against a single frequency, possibly with recovery data too.
% BPG 20-5-21

% Variables to set:
% smod - synapse model
% freq - stimulation frequency (constant rate spike train)
% snum - number of stimuli
% rect - recovery times to test
% fstem - file stem of experimental data for fitting
% fstemr - file stem of experimental recovery data for fitting

% Synapse parameters
% Pv0 - initial release probability
% P1 - increment in Pv0 per AP (facilitation)
% tauf - relaxation time constant of facilitation (msecs) 
% trB - time constant of background recovery (msecs)
% trH - recovery time constant instantly following a spike
% trR - time constant of relaxation of tr back to background rate
% D - fraction of desensitised receptors
% tauD - time constant of recovery from desensitisation (msecs)
% spt - vector of spike times (msecs)
% ndat - normalised experimental EPSC amplitudes

global Pv0 P1 tauf trB trH trR D tauD spt ndat ndatrec

fdat = 1;   % stimulation data
frecov = 1; % recovery data
if (fdat == 1)
    %fstem = '../Data/Kv33_100Hz_control.txt';
    %fsout = '../Results/Kv33mod_100Hz_control1_stim.txt';
    %fsout = '../Results/Kv33mod_100Hz_control2_stim.txt';
    fstem = '../Data/Kv33_100Hz_Kv33KO.txt';
    fsout = '../Results/Kv33mod_100Hz_Kv33KO1_stim.txt';
    %fsout = '../Results/Kv33mod_100Hz_Kv33KO2_stim.txt';
    expdat = load(fstem);
    ndat = expdat(:,1);  % responses
    se = expdat(:,2);  % standard errors
end;
if (frecov == 1)
    %fstem = '../Data/Kv33_100Hz_control2_recov.txt';
    %frout = '../Results/Kv33mod_100Hz_control2_recov.txt';
    %frout = '../Results/Kv33mod_100Hz_controlTst_recov.txt';
    fstem = '../Data/Kv33_100Hz_Kv33KO2_recov.txt';
    %frout = '../Results/Kv33mod_100Hz_Kv33KO2_recov.txt';
    frout = '../Results/Kv33mod_100Hz_Kv33KOTst_recov.txt';
    expdat = load(fstem);
    ndatrec = expdat(:,1);  % responses
    serec = expdat(:,2);  % standard errors
end;

% Values roughly from Graham, Wong & Forsythe, Neural Computing, 2004
%Pv0 = 0.13;   % WT
Pv0 = 0.266;    % KO
P1 = 0; % facilitation
tauf = 100;
%trH = 66.9;  % fast control (WT) rate
trH = 52.2;  % fast KO rate
trB = 3000; % background rate
trR = 400; % rate of relaxation to background rate
D = 0;  % desensitization
%D = 1;  % desensitization
tauD = 100;

%Generate spikes at different frequencies
slen = 800; % stimulation time (msecs)
freq = 100; % stimulation frequency (Hertz)
isi = 1000/freq;	% interspike interval (msecs)
last = slen - rem(slen,isi);
spt = [isi:isi:last];	% spike times (msecs)
spcnt = (last/isi)+1;	% no. of spikes
rect = [0.05, 0.1, 0.5, 1, 2, 5, 10, 20, 30]*1000+slen;  % recovery times (s)
if (frecov == 1)
    sptrec = [spt rect];
%    spcnt = length(sptrec);
    spt = sptrec;
%    last = spt(length(spt));
end;


% Optimise model against data
fopt = 0;
if (fopt == 1)
%    [op,fval,exitflag,output] = fminsearch(@fit_syn_dynfrec_single, [Pv0, trH], ...
%        optimset('MaxFunEvals', 1e9, 'MaxIter', 1e9));
    [op,fval,exitflag,output] = fminsearch(@fit_syn_dynfrec_single, [Pv0, trH, trR, trB], ...
        optimset('MaxFunEvals', 1e9, 'MaxIter', 1e9));
    Pv0=op(1)
    trH=op(2)
    trR=op(3)
    trB=op(4)
    fval
    output.funcCount
end;

% Get final model results
[n,Pv,frD,psr,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt);

% RMSE over recovery data
ndatr = [ndat(length(ndat)) ndatrec'];   % include final stim point
sqerecov = sqrt(sum((ndatr-(psr(length(ndat):length(ndat)+length(ndatrec))./psr(1))).^2))
%f = sqrt(sum((ndatall-(psr./psr(1))).^2));

% Plot postsynaptic responses
figure();
%if (frecov == 1)
%    subplot(2,1,1);
%end;
if (fdat == 1)
    mline=errorbar(spt(1:spcnt-1),ndat,se,'c-');
    set(mline,'LineWidth',1.5);
    hold on;
end;
mline=plot(spt(1:spcnt-1),psr(1:spcnt-1)./psr(1), 'k-');
set(mline,'LineWidth',1.5);
axis([0 last 0 1]);
title('Normalised Postsynaptic Response');
xlabel('Time (msecs)');
ylabel('PSR');
% save model data to file
dout = [spt(1:spcnt-1); psr(1:spcnt-1)./psr(1)];
fso=fopen(fsout, 'w');
fprintf(fso, '%f %f\n', dout);
fclose(fso);

if (frecov == 1)
    %subplot(2,1,2);
    figure();
    subplot(2,1,1); % full recovery
    mline=errorbar(spt(spcnt:length(spt)),ndatrec,serec,'c-');
    set(mline,'LineWidth',1.5);
    hold on;
    mline=plot(spt(spcnt-1:length(spt)),psr(spcnt-1:length(spt))./psr(1), 'kx-');
    set(mline,'LineWidth',1.5);
    hold on;
    axis([0 rect(length(rect))*1.1 0 1.1]);
    title('Normalised Recovery Response');
    %xlabel('Time (msecs)');
    ylabel('PSR');
    
    subplot(2,1,2); % early recovery
    nr=6;
    mline=errorbar(spt(spcnt:spcnt+nr-1),ndatrec(1:nr),serec(1:nr),'c-');
    set(mline,'LineWidth',1.5);
    hold on;
    mline=plot(spt(spcnt-1:spcnt+nr-1),psr(spcnt-1:spcnt+nr-1)./psr(1), 'kx-');
    set(mline,'LineWidth',1.5);
    hold on;
    axis([slen (rect(nr)-slen)*1.2 0 1.1]);
    %title('Normalised Recovery Response');
    xlabel('Time (msecs)');
    ylabel('PSR');
    
    % save model data to file
    dout = [spt(spcnt-1:length(spt)); psr(spcnt-1:length(spt))./psr(1)];
    fso=fopen(frout, 'w');
    fprintf(fso, '%f %f\n', dout);
    fclose(fso);

end;



