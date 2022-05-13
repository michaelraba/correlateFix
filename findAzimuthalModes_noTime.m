function [qq]=findAzimuthalModes_noTime(currentCrossSec, qMinusQbar_noCsYet,azimuthDoneXcorrDone,aliasStr)
 [ntimesteps, rMin, rMax, ss, ncs, plotOn, azimuthalSet ,azimuthalSetSize ,printStatus ,lags]=constants();

  %sprintf('%s','pause')
  [postAzimuthFft_noCsYet]=initData2("postAzimuthFft_noCsYet");


if aliasStr=="noAlias"
    parfor t = 1:ntimesteps % time
    for m = 1:540 % time
    postAzimuthFft_noCsYet(t).circle(m).dat=fft(qMinusQbar_noCsYet(t).circle(m).dat); % this can perhaps be truncated to 540, then duplicated for the second half, too prevent aliasing.!
    end % m
    end % for t
    qq = postAzimuthFft_noCsYet;
elseif aliasStr=="alias"
    %for t = 1:ntimesteps % time
    for m = 1:540 % time
    aa=qMinusQbar_noCsYet(1).circle(m).dat(1:end/2); % this can perhaps be truncated to 540, then duplicated for the second half, too prevent aliasing.!
    aa=fft(aa);
    bb = flip(aa);
    cc = zeros(1080,1);
    for i=1:540
      cc(i) =aa(i);
      cc(1080 - i + 1 ) = aa(i); % get all 1080
    end
    postAzimuthFft_noCsYet(1).circle(m).dat=cc;
    end % m
    hold on;
    plot(real(cc))
    end % for t
    ordStr="xcorrNow";
    if ordStr=="xcorrNow"
        %sprintf('%s','pause')
        for t=1:ntimesteps%
        vec = zeros(1,540); % collect radial points..
        for m=1:azimuthalSetSize
        mm = azimuthalSet(m);
        for r=1:540% size xcorr
          aa = postAzimuthFft_noCsYet(1).circle(r).dat(mm,1);
          vec(r) = aa;
end % r
          [bb, lags] = xcorr(vec,"normalized");
        for r=1:1079
            azimuthDoneXcorrDone(1).azimuth(m).radial(r).dat=bb;
        end %r
end % m
%end % t % parfor time loop we are removing
%sprintf('%s','pause')
    
    end % if ordStr
end % if big
qq = azimuthDoneXcorrDone; % asign qq and exi
end % fc
