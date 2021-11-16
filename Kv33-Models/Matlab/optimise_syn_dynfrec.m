% Optimise synapse model against experimental time series
% Optimise against multiple frequencies at once
% BPG 19-2-21

% Variables to set:
% smod - synapse model
% freq - stimulation frequency (constant rate spike train)
% snum - number of stimuli
% rect - recovery times to test
% fstem - file stem of experimental data for fitting

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
% normdat - normalised experimental EPSC amplitudes

global Pv0 P1 tauf trB trH trR D tauD spt100 spt200 spt600 ndat100 ndat200 ndat600

%fstem = '../Data/Kv33_100Hz_control.txt';
%fstem = '../Data/Kv33_100Hz_Kv33KO.txt';
%fstem = '../Data/Kv33_200Hz_control.txt';
%fstem = '../Data/Kv33_200Hz_Kv33KO.txt';
%fstem = '../Data/Kv33_600Hz_control.txt';

% Load experimental data
fcont = 0;
if (fcont == 1)
    fstem = '../Data/Kv33_100Hz_control.txt';
    %fstem = '../Data/Kv33_100Hz_Kv33KO.txt';
    expdat = load(fstem);
    ndat100 = expdat(:,1);  % responses
    se100 = expdat(:,2);  % standard errors
    fstem = '../Data/Kv33_200Hz_control.txt';
    %fstem = '../Data/Kv33_200Hz_Kv33KO.txt';
    expdat = load(fstem);
    ndat200 = expdat(:,1);  % responses
    se200 = expdat(:,2);  % standard errors
    fstem = '../Data/Kv33_600Hz_control.txt';
    %fstem = '../Data/Kv33_600Hz_Kv33KO.txt';
    expdat = load(fstem);
    ndat600 = expdat(:,1);  % responses
    se600 = expdat(:,2);  % standard errors
else
    %fstem = '../Data/Kv33_100Hz_control.txt';
    fstem = '../Data/Kv33_100Hz_Kv33KO.txt';
    expdat = load(fstem);
    ndat100 = expdat(:,1);  % responses
    se100 = expdat(:,2);  % standard errors
    %fstem = '../Data/Kv33_200Hz_control.txt';
    fstem = '../Data/Kv33_200Hz_Kv33KO.txt';
    expdat = load(fstem);
    ndat200 = expdat(:,1);  % responses
    se200 = expdat(:,2);  % standard errors
    %fstem = '../Data/Kv33_600Hz_control.txt';
    fstem = '../Data/Kv33_600Hz_Kv33KO.txt';
    expdat = load(fstem);
    ndat600 = expdat(:,1);  % responses
    se600 = expdat(:,2);  % standard errors
end;

% Values roughly from Graham, Wong & Forsythe, Neural Computing, 2004
%Pv0 = 0.13;   % control (WT) all opt
Pv0 = 0.266;    % KO all opt
%Pv0 = 0.1293;   % WT 100Hz+recov opt
%Pv0 = 0.2624;    % KO 100Hz+recov opt
P1 = 0; % facilitation
tauf = 100;
%trH = 66.9;  % fast WT rate all opt
trH = 52.2;  % fast KO rate all opt
%trH = 42.14;  % fast WT rate 100Hz+recov opt
%trH = 44.1145;  % fast KO rate 100Hz+recov opt
trR = 400; % rate of relaxation to background rate (both all opt)
%trR = 131.0415; % rate of relaxation to background rate (WT 100Hz+recov opt)
%trR = 365.5246; % rate of relaxation to background rate (KO 100Hz+recov opt)
trB = 3000; % background rate (both all opt)
%trB = 2130; % background rate (WT 100Hz+recov opt)
%trB = 3409; % background rate (KO 100Hz+recov opt)
D = 0;  % desensitization
%D = 1;  % desensitization
tauD = 100;

%Generate spikes at different frequencies
slen = 800; % stimulation time (msecs)
freq = 100; % stimulation frequency (Hertz)
isi = 1000/freq;	% interspike interval (msecs)
last = slen - rem(slen,isi);
spt100 = [isi:isi:last];	% spike times (msecs)
spcnt100 = (last/isi)+1;	% no. of spikes
freq = 200; % stimulation frequency (Hertz)
isi = 1000/freq;	% interspike interval (msecs)
last = slen - rem(slen,isi);
spt200 = [isi:isi:last];	% spike times (msecs)
spcnt200 = (last/isi)+1;	% no. of spikes
freq = 600; % stimulation frequency (Hertz)
isi = 1000/freq;	% interspike interval (msecs)
last = slen - rem(slen,isi);
spt600 = [isi:isi:last];	% spike times (msecs)
spcnt600 = (last/isi)+1;	% no. of spikes

% Optimise model against data at all frequencies
fopt = 0;
if (fopt == 1)
%    [op,fval,exitflag,output] = fminsearch(@fit_syn_dynfrec, [Pv0, trH], ...
%        optimset('MaxFunEvals', 1e9, 'MaxIter', 1e9));
    [op,fval,exitflag,output] = fminsearch(@fit_syn_dynfrec, [Pv0, trH, trR, trB], ...
        optimset('MaxFunEvals', 1e9, 'MaxIter', 1e9));
    Pv0=op(1)
    trH=op(2)
    trR=op(3)
    trB=op(4)
    fval
    output.funcCount
end;

% Get final model results
%[n,Pv,frD,psr,Pr] = syn_fpfad2(Pv0,P1,tauf,taur,ns,nrp,D,tauD,spt);
[n,Pv,frD,psr100,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt100);
[n,Pv,frD,psr200,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt200);
[n,Pv,frD,psr600,Pr] = syn_dynfrec(Pv0,P1,tauf,trB,trH,trR,D,tauD,spt600);
% save model data to file
%fsout = '../Results/Kv33mod_control_optall_';
fsout = '../Results/Kv33mod_Kv33KO_optall_';
%fsout = '../Results/Kv33mod_control_opt100Hz_';
%fsout = '../Results/Kv33mod_Kv33KO_opt100Hz_';
dout = [spt100; psr100./psr100(1)];
fso=fopen([fsout '100Hz.txt'], 'w');
fprintf(fso, '%f %f\n', dout);
fclose(fso);
dout = [spt200; psr200./psr200(1)];
fso=fopen([fsout '200Hz.txt'], 'w');
fprintf(fso, '%f %f\n', dout);
fclose(fso);
dout = [spt600; psr600./psr600(1)];
fso=fopen([fsout '600Hz.txt'], 'w');
fprintf(fso, '%f %f\n', dout);
fclose(fso);

% Plot postsynaptic responses
%mline=plot(spt(1:lresp),100*psr(1:lresp)./psr(1),lines{li});
%mline=plot(spt(1:snum),normdat(1:snum),'r-');
mline=errorbar(spt100,ndat100,se100,'c-');
set(mline,'LineWidth',1.5);
hold on;
mline=plot(spt100,psr100./psr100(1), 'k-');
set(mline,'LineWidth',1.5);
mline=errorbar(spt200,ndat200,se200,'c-');
set(mline,'LineWidth',1.5);
mline=plot(spt200,psr200./psr200(1), 'k-');
set(mline,'LineWidth',1.5);
mline=errorbar(spt600,ndat600,se600,'c-');
set(mline,'LineWidth',1.5);
mline=plot(spt600,psr600./psr600(1), 'k-');
set(mline,'LineWidth',1.5);
axis([0 last 0 1]);
title('Normalised Postsynaptic Response');
xlabel('Time (msecs)');
ylabel('PSR');


