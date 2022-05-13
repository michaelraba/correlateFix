function [qq]=findAzimuthalModes(currentTime, currentCrossSec, qMinusQbar_noCsYet,azimuthDoneXcorrDone,aliasStr)
% [ntimesteps, rMin, rMax, ss, ncs, plotOn, azimuthalSet ,azimuthalSetSize ,printStatus ,lags]=constants();
  [ntimesteps, rMin, rMax, ss, ncs, plotOn, azimuthalSet ,azimuthalSetSize ,printStatus ,lags, blocLength, saveDir]=constants();

  [postAzimuthFft_noCsYet]=initData2("postAzimuthFft_noCsYet");

if aliasStr=="noAlias"
    parfor t = 1:ntimesteps % time
    for m = 1:540 % time
    postAzimuthFft_noCsYet(t).circle(m).dat=fft(qMinusQbar_noCsYet(t).circle(m).dat); % this can perhaps be truncated to 540, then duplicated for the second half, too prevent aliasing.!
    end % m
    end % for t 
    qq = postAzimuthFft_noCsYet;
elseif aliasStr=="alias"
    % do fft for the first half of the circle, then copy the result ot the other half.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Load in the correct qMinusqBar..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    clear qMinusQbar_noCsYet;
    for timeBloc = 1:blocLength% time
        saveStr=['/mnt/archLv/mike/podData/apr18/qMinusQbar[Case]C' num2str(ncs) 'T' num2str(ntimesteps) '[crossSec]' num2str(currentCrossSec) '[TimeBloc]' num2str(timeBloc) '.mat'       ];
        load(saveStr,'qMinusQbar_noCsYet');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *timebloc* should operate on both azimuthal mode and xcorr.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor t = 1:ntimesteps % time
        for m = 1:540 % time
            % first half of string..:
            aa=qMinusQbar_noCsYet(t).circle(m).dat(1:end/2); % this can perhaps be truncated to 540, then duplicated for the second half, too prevent aliasing.!
            aa=fft(aa);
            bb = flip(aa);
            cc = zeros(1080,1);
            for i=1:540
              cc(i) =aa(i);
              cc(1080 - i + 1 ) = aa(i); % get all 1080
            end % i
            postAzimuthFft_noCsYet(t).circle(m).dat=cc;
        end % m ...
        hold on;
        plot(real(cc))
    end % parfor t
    ordStr="xcorrNow";
    if ordStr=="xcorrNow"
        for t=1:ntimesteps%
            vec = zeros(1,540); % collect radial points..
            %for m=1:azimuthalSetSize%
                for m=1:azimuthalSetSize
                mm = azimuthalSet(m);
                for r=1:540% size xcorr
                  aa = postAzimuthFft_noCsYet(t).circle(r).dat(mm,1);
                  vec(r) = aa;
                end % r
                  [bb, lags] = xcorr(vec,"normalized");

                azimuthDoneXcorrDone(t).circle(m).dat=bb';
                
                % this is not init form:
                %azimuthDoneXcorrDone(t).azimuth(m).radial(r).dat=bb;
                %for r=1:1079
                %    azimuthDoneXcorrDone(t).azimuth(m).radial(r).dat=bb;
                %end %r
                end % m
        end % (little)t
     else % empty else-end

    

     end % (timebloc)timebloc
     saveStr=['/mnt/archLv/mike/podData/apr18/azimuthDoneXcorrDone[Case]C' num2str(ncs) 'T' num2str(ntimesteps) '[crossSec]' num2str(currentCrossSec) '[TimeBloc]' num2str(timeBloc) '.mat'       ];
     save(saveStr,'azimuthDoneXcorrDone','-v7.3');
  end % if ordStr
end % if big
qq = azimuthDoneXcorrDone; % asign qq and exi
end % fc
