function [qq]=findAzimuthalModes2(currentTime, currentCrossSec, qMinusQbar_noCsYet,xcorrDone,aliasStr)
% [ntimesteps, rMin, rMax, ss, ncs, plotOn, azimuthalSet ,azimuthalSetSize ,printStatus ,lags]=constants();
  [ntimesteps, rMin, rMax, ss, ncs, plotOn, azimuthalSet ,azimuthalSetSize ,printStatus ,lags, blocLength, saveDir]=constants();

  [postAzimuthFft_noCsYet]=initData2("postAzimuthFft_noCsYet");

if aliasStr=="noAlias"
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
%    for t = 1:ntimesteps % time % parfor
%        for m = 1:540 % time
%            % first half of string..:
%            aa=qMinusQbar_noCsYet(t).circle(m).dat(1:end/2); % this can perhaps be truncated to 540, then duplicated for the second half, too prevent aliasing.!
%            aa=fft(aa);
%            bb = flip(aa);
%            cc = zeros(1080,1);
%            for i=1:540
%              cc(i) =aa(i);
%              cc(1080 - i + 1 ) = aa(i); % get all 1080
%            end % i
%            postAzimuthFft_noCsYet(t).circle(m).dat=cc;
%        end % m ...
%        hold on;
%        plot(real(cc))
%    end % parfor t
    ordStr="xcorrNow";
    if ordStr=="xcorrNow"
        parfor t=1:ntimesteps%
            vec = zeros(1,540); % collect radial points..
            %for m=1:azimuthalSetSize%
                for m=1:1080
                %mm = azimuthalSet(m);
                for r=1:540% size xcorr
                  %aa = postAzimuthFft_noCsYet(t).circle(r).dat(mm,1);
                  % dont need to do this anymore. 
                  %aa = qMinusQbar_noCsYet(t).circle(mm).dat(r,1);
                  aa = qMinusQbar_noCsYet(t).circle(m).dat(r,1);

                  vec(r) = aa;
                end % r
                  [bb, lags] = rcorr(vec,"normalized"); % bb is 1079 because of xcorr ! <- new annotat.
                xcorrDone(t).circle(m).dat=bb';
                end % m
        end % (little)t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% begin azimuthal
% Note: need to adjust aa= name etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

    parfor t = 1:ntimesteps % time % parfor
        for  r = 1:1079 % 1079 because of xcorr has 2x-1 entries..
            vec = zeros(1080,1);
            vec2 = zeros(1080,1);
   
            for zz=1:1080 % there are currently 1080 azimuthal modes. 
            aa=xcorrDone(t).circle(zz).dat(r); % this can perhaps be truncated to 540, then duplicated for the second half, too prevent aliasing.!
            vec(zz)= aa;
            end % for zz 
            aa=fft(vec);
            bb = flip(aa);
            cc = zeros(1080,1);
            for i=1:540
              cc(i) =aa(i);
              cc(1080 - i + 1 ) = aa(i); % get all 1080
            end % i
            postAzimuthFft_noCsYet(t).circle(1,r).dat=cc;
            %plot(real(cc))
            %pause(0.1);
        end % r...
        %hold on;
        
    end % parfor t




     else % empty else-end
    end % end if ordStr
     saveStr=['/mnt/archLv/mike/podData/apr18/xcorrDone[Case]C' num2str(ncs) 'T' num2str(ntimesteps) '[crossSec]' num2str(currentCrossSec) '[TimeBloc]' num2str(timeBloc) '.mat'       ];
     save(saveStr,'xcorrDone','-v7.3');
    end % for big timeBloc
end % if aliasStr
qq = xcorrDone; % asign qq and exi
end % fc