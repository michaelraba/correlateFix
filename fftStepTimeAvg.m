 function [qq]=fftStepTimeAvg(stepStr,preStr)
[ntimesteps, rMin, rMax, ss, ncs, plotOn, azimuthalSet ,azimuthalSetSize ,printStatus ,lags, blocLength]=constants();
  [azimuthDoneXcorrDoneAnticipate_cs]=initData2("azimuthDoneXcorrDoneAnticipate_cs");
if stepStr=="readDataAndFindVeloFluctuation"
    [qMinusQbar_noCsYet]=initData2("qMinusQbar_noCsYet"); % initialize avg struct
    [qMinusQbar]=initData2("qMinusQbar"); % initialize avg struct
    [myPreFft_noCsYet]=initData2("myPreFft_noCsYet");
    [avgPreFft_noCsYet]=initData2("avgPreFft_noCsYet");
    [azimuthDoneXcorrDone]=initData2("azimuthDoneXcorrDone");
    [qMinusQbar_noCsYet]=initData2("qMinusQbar_noCsYet");
    [xdirNew]=initData2("xdirNew");
    [xdirPostFft]=initData2("xdirPostFft");
    [avgTimeEnd]=initData2("avgTimeEnd");

    for c = ncs:ncs  % crosssection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x-dir fft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in one of the saved azimuthDoneXcorrDone
for timeBloc=1:blocLength
        sprintf('%s','opening xdirPostFft...')
        saveStr=['/mnt/archLv/mike/podData/apr14/xdirPostFft[Case]C' num2str(ncs) 'T' num2str(ntimesteps) '[crossSec]' num2str(c) '[TimeBloc]' num2str(timeBloc) '.mat'       ];
        %save(saveStr,'xdirPostFft','-v7.3');
        qqq=open(saveStr);
        sprintf('%s%s','Saved xdirpostfft into file ',saveStr);
        xdirPostFft = qqq.xdirPostFft;
% Time Averaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

for t=1:ntimesteps
for c=1:ncs
for m=1:azimuthalSetSize
for r=1:1079
   aa = xdirPostFft(t).RadialCircle(r).azimuth(m).dat(c,1);
if t<ntimesteps || timeBloc ~= blocLength
   avgTimeEnd(c).circle(m).dat(r,1) = avgTimeEnd(c).circle(m).dat(r,1) + aa;
elseif t==ntimesteps && timeBloc == blocLength
   avgTimeEnd(c).circle(m).dat(r,1) = (avgTimeEnd(c).circle(m).dat(r,1) + aa)/(ntimesteps*blocLength);
end  % if
end % r
end % m
end % c (little)
end % t (little)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
end % timeBloc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%

% set back in radial direction and time avergae for all timesteps!




qq = xdirPostFft;
        
end % if
end % f
