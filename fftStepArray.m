 function [qq]=fftStepArray(stepStr,preStr)
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

    for c = 1:ncs  % crosssection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step A) load a chonk into memory and read circles in.
    for timeBloc=1:blocLength
    parfor t = 1:ntimesteps % time % <-- nb, this is the parfor loop.
    myPreFft_noCsNoTimeYet=readCircles2(timeBloc*t,c);
    myPreFft_noCsYet(t).circle=myPreFft_noCsNoTimeYet;
    sprintf('%s','pause')
    end % parfor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % after each block is done, find the average
    for t = 1:ntimesteps % time % <-- nb, this is the parfor loop.
      if timeBloc==blocLength && t==ntimesteps
        lastStr="last";
      else
        lastStr="notLast";
      end
    [avgPreFft_noCsYet]=findQbar(t,c,myPreFft_noCsYet,avgPreFft_noCsYet,lastStr); % find temporal average.
    end % end little t
    end % timeBloc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for timeBloc=1:blocLength
    parfor t = 1:ntimesteps % time
% load in time bloc again
        myPreFft_noCsNoTimeYet=readCircles3(timeBloc*t,c);
        myPreFft_noCsYet(t).circle=myPreFft_noCsNoTimeYet;
        % for each loaded timebloc, find qMinusQbar..
        [ qMinusQbar_noCsYet(t) ]=FindqMinusQbar(t,c,myPreFft_noCsYet(t),avgPreFft_noCsYet,qMinusQbar_noCsYet(t),"efficient");
    end % parfor t
    % this should be saved to disk immediately, for each timeBloc...
        sprintf('%s','saving qMinusQbar...')
        saveStr=['/mnt/archLv/mike/podData/apr14/qMinusQbar[Case]C' num2str(ncs) 'T' num2str(ntimesteps) '[crossSec]' num2str(c) '[TimeBloc]' num2str(timeBloc) '.mat'       ];
        save(saveStr,'qMinusQbar_noCsYet','-v7.3');
        sprintf('%s%s','Saved velocity fluctuations into file ',saveStr);

    qq = qMinusQbar_noCsYet(t);
    end % timeblock

    % nb time looping within..
    %[azimuthDoneXcorrDone]=findAzimuthalModes(c, qMinusQbar_noCsYet,azimuthDoneXcorrDone,"alias")
    [azimuthDoneXcorrDone]=findAzimuthalModes(t,c, qMinusQbar_noCsYet,azimuthDoneXcorrDone,"alias")
    sprintf('%s','start azimuthal')
    qq = azimuthDoneXcorrDone;

        end %c
        %elseif stepStr=="azimuth"
        end % if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x-dir fft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in one of the saved azimuthDoneXcorrDone
for timeBloc=1:blocLength
for currentCrossSec=1:ncs
saveStr=['/mnt/archLv/mike/podData/apr14/azimuthDoneXcorrDone[Case]C' num2str(ncs) 'T' num2str(ntimesteps) '[crossSec]' num2str(currentCrossSec) '[TimeBloc]' num2str(timeBloc) '.mat'       ];
open(saveStr);
sprintf('%s','start azimuthal')
% now re-organize:
for t=1:ntimesteps
for r=1:1079
for m=1:azimuthalSetSize
  aa=azimuthDoneXcorrDone(t).circle(m).dat(r,1);
  %ab= xdirNew(t).RadialCircle(r).azimuth(m).dat(currentCrossSec,1);
  xdirNew(t).RadialCircle(r).azimuth(m).dat(currentCrossSec,1) = aa;
end % m
end % r
end % t (little)
sprintf('%s%d%s%d%s','done filling in a crosssec for timeBloc=', timeBloc, ' and t=',t,'.')
end % c

% begin fft x-dir
for t=1:ntimesteps
for r=1:1079
for m=1:azimuthalSetSize
  aa = xdirNew(t).RadialCircle(r).azimuth(m).dat;
  %ab = fft(aa(end/2:end));
  ab = fft(aa);


  xdirPostFft(t).RadialCircle(r).azimuth(m).dat = ab;
  %hold on;
  %plot(real(ab));
  %pause(0.1)
end % m
end % r
end % t (little)
        sprintf('%s','saving xdirPostFft...')
        saveStr=['/mnt/archLv/mike/podData/apr14/xdirPostFft[Case]C' num2str(ncs) 'T' num2str(ntimesteps) '[crossSec]' num2str(c) '[TimeBloc]' num2str(timeBloc) '.mat'       ];
        save(saveStr,'xdirPostFft','-v7.3');
        sprintf('%s%s','Saved xdirpostfft into file ',saveStr);



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
end % c
end % r
end % m (little)
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
        
end % f
