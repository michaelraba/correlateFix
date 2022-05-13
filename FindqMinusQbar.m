function  [ qMinusQbar_noCsYet]=FindqMinusQbar(~,~,myPreFft_noCsYet,avgPreFft_noCsYet,~,performanceString)
  % find difference AND take fft.
[~ ,~ ,~, ~ ,~ ,~ ,~, ~ ,~, ~]=constants();

%fftMatQminusQbar=[];


if performanceString=="efficient"

 % Take FFT of Difference q - \bar{q}
for m=1:1080 % m=1:540
    % set equal to 1 because this function is given 1 timestep only.
    qMinusQbar_noCsYet(1).circle(m).dat = myPreFft_noCsYet(1).circle(m).dat - avgPreFft_noCsYet(m).dat  ; % 3 corresponds ot the streamwise data; 1 2 and are the x coordinates!
end %m

end % if
end %fcn
