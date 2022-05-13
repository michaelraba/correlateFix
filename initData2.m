function [qq]=initData2(initStr) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
  [ntimesteps ,~ ,~, ~ ,ncs ,~ ,~, azimuthalSetSize ,printStatus, ~]=constants();
Ncs= [ncs,1];
Nts= [ntimesteps,1];
Naz = [azimuthalSetSize,1]; % azimuthal
Nps=[540,1];
  if printStatus=="on" 
    %sprintf('%s%s', '* Initializing ', initStr)
  end
if initStr=="myPreFftDrawCircle" %redo
qq=struct('t', repmat({struct('circle', repmat({  struct('dat',repmat({zeros(3,1080)}, Nps))}, Nts)) }, Ncs));
elseif initStr=="myPreFft" %redo
qq=struct('t', repmat({struct('circle', repmat({  struct('dat',repmat({zeros(1,1080)}, Nps))}, Nts)) }, Ncs));
elseif initStr=="myPreFft_noTimeYet" %redo
%qq=struct('circle', repmat({  struct('dat',repmat({zeros(1,1080)}, Nps))}, Nts))
%qq=struct('dat', repmat({  struct('dat',repmat({zeros(1,1080)}, Nps))}, Nts))
qq=struct('dat', repmat({zeros(1,1080)}, [1,540]));
elseif initStr=="myPreFft_noCsYet" %redo
qq=struct('circle', repmat({struct('dat',repmat({zeros(1,1080)}, [1,540]))} , [1,ntimesteps]));
elseif initStr=="postAzimuthFft_noCsYet" %redo
qq=struct('circle', repmat({struct('dat',repmat({zeros(1,1080)}, [1,540]))} , [1,ntimesteps]));
elseif initStr=="azimuthDoneXcorrDone" % for findAzimuthalModes with *if ordStr="xcorrNow"*
qq=struct('circle', repmat({struct('dat',repmat({zeros(1080,1)}, [1,1080]))} , [1,ntimesteps]));
elseif initStr=="azimuthDoneXcorrDoneAnticipate_cs" % for findAzimuthalModes with *if ordStr="xcorrNow"*
%qq=struct('radial', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(1,100)}, [1,azimuthalSetSize]))}, [1,1079])) }, Nts));
qq=struct('azimuth', repmat({struct('radial', repmat({  struct('dat',repmat({zeros(1,100)}, [1,1079]))}, [1,azimuthalSetSize])) }, Nts));
%azimuthDoneXcorrDoneAnticipate_cs


elseif initStr=="xdirNew"
qq=struct('RadialCircle', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(ncs,1)}, Naz))}, [1,1079])) }, Nts));
elseif initStr=="xdirPostFft"
qq=struct('RadialCircle', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(ncs,1)}, Naz))}, [1,1079])) }, Nts));

elseif initStr=="avgTimeEnd" %redo
qq=struct('circle', repmat({struct('dat',repmat({zeros(1079,1)}, [1,azimuthalSetSize]))} , [1,ncs]));

elseif initStr=="myPreFft_noCsNoTimeYet" %redo
qq=struct('dat', repmat({zeros(1,1080)}, [540,1]));

elseif initStr=="avgPreFft" %redo
%qq=struct('circle', repmat({struct('dat',repmat({zeros(540,1)}, [1,1080]))} , [1,ncs]));
qq=struct('circle', repmat({struct('dat',repmat({zeros(540,1)}, [1080,1]))} , [1,ncs]));

elseif initStr=="avgPreFft_noCsYet" %redo
  %qq=struct('dat',repmat({zeros(540,1)}, [1,540]));
  qq=struct('dat',repmat({zeros(540,1)}, [1080,1]));
elseif initStr=="qMinusQbar"
qq=struct('t', repmat({struct('circle', repmat({  struct('dat',repmat({zeros(1,1080)}, Nps))}, Nts)) }, Ncs));

elseif initStr=="qMinusQbar_noCsYet"
  % too many nesting here.
qq=struct('circle', repmat({struct('dat',repmat({zeros(540,1)}, [1080,1]))} , [1,ntimesteps]));
  %qq=struct('dat',repmat({zeros(1,1080)}, [1,540]));

elseif initStr=="qMinusQbar_noCsNoTimeYet" % for truncating out the time when we feed eahc itme step to the azimuth fucnciton.
  qq=struct('dat',repmat({zeros(540,1)}, [1080,1]));

elseif initStr=="radialXcorrStruct" % redo
qq=struct('crossSec', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(1,540)}, Naz))}, Ncs)) }, Nts));
elseif initStr=="XcorrData"
%%%%%    XcorrData(t).crossSec(c).azimuth(m).dat = zeros(1,1080);
qq=struct('crossSec', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(1,1080)}, Naz))}, Ncs)) }, Nts));
elseif initStr=="XcorrDonePreCsFft"
qq=struct('radial', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(1,ncs)}, Naz))}, [1,1080])) }, Nts));
elseif initStr=="XcorrDonePostCsFft"
qq=struct('radial', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(1,ncs)}, Naz))}, [1,1080])) }, Nts));
elseif initStr=="XcorrDoneCsFftDonePreAzimFft"
qq=struct('radial', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(1,azimuthalSetSize)}, Ncs))}, [1,1080])) }, Nts));
elseif initStr=="XcorrDoneCsFftDonePostAzimFft"
qq=struct('crossSec', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(1,azimuthalSetSize)}, Ncs))}, [1,1080])) }, Nts));
elseif initStr=="XcorrAzimDonePreTimeAvg"
qq=struct('radial', repmat({struct('crosssection', repmat({  struct('dat',repmat({zeros(1,ntimesteps)}, Ncs))}, [1,1080])) }, Naz));
elseif initStr=="XcorrAzimDonePostTimeAvg"
qq=struct('radial', repmat({struct('crosssection', repmat({  struct('dat',repmat({zeros(1,ntimesteps)}, Ncs))}, [1,1080])) }, Naz));
elseif initStr=="Skmr"
qq=struct('radial', repmat({struct('crosssection', repmat({  struct('dat',repmat({zeros(1,ntimesteps)}, Ncs))}, [1,1])) }, Naz));
elseif initStr=="streamwiseFft"
qq=struct('RadialCircle', repmat({struct('azimuth', repmat({  struct('dat',repmat({zeros(1,540)}, Naz))}, [1,540])) }, Nts));
elseif initStr=="SrrPrForFourierPre"
qq=struct('t', repmat({struct('azimuthal', repmat({  struct('dat',repmat({zeros(1,540)}, Naz))}, Nts)) }, Ncs));
elseif initStr=="SrrPrForFourierPost"
qq=struct('t', repmat({struct('azimuthal', repmat({  struct('dat',repmat({zeros(1,540)}, Naz))}, Nts)) }, Ncs));

end %if

end %fc
