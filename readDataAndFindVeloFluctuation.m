function [qMinusQbar_noCsYet ]= readDataAndFindVeloFluctuation()
[ntimesteps ,~ ,~, ~ ,ncs ,~ ,~, ~ ,~, ~]=constants();
  [qMinusQbar_noCsYet]=initData2("qMinusQbar_noCsYet"); % initialize avg struct

for c = 1:ncs  % crosssection
[myPreFft_noCsYet]=initData2("myPreFft_noCsYet");
parfor t = 1:ntimesteps % time
myPreFft_noCsNoTimeYet=readCircles2(t,c);
myPreFft_noCsYet(t).circle=myPreFft_noCsNoTimeYet;

sprintf('%s','pause')
end %t >-

% compute average for each cross section.
[avgPreFft_noCsYet]=initData2("avgPreFft_noCsYet"); % initialize avg struct

for t = 1:ntimesteps % time
[avgPreFft_noCsYet]=findQbar(t,c,myPreFft_noCsYet,avgPreFft_noCsYet); % find temporal average.
end %t
sprintf('%s','initialize QminusQbar for lot of timestep.')

%sprintf('%s','findqminusqbar')
[qMinusQbar_noCsYet]=initData2("qMinusQbar_noCsYet");
%for c = 1:ncs  % crosssection
parfor t = 1:ntimesteps % time
    [ qMinusQbar_noCsYet(t) ]=FindqMinusQbar(t,c,myPreFft_noCsYet(t),avgPreFft_noCsYet,qMinusQbar_noCsYet(t),"efficient");
end %t
%end %c

end % f
