function xcorrGpu()
clear all;
maxJ = 7;
maxI = 3;
%qq=struct('t', repmat({struct('circle', repmat({  struct('dat',repmat({zeros(3,1079)}, Nps))}, Nts)) }, Ncs));
qq=struct('t', repmat({zeros(1,7)}, [1,maxJ]));

for i=1:7
  a= zeros(maxJ,1) ;
for j=1:maxJ
  sprintf('%s','hi');
  a(j) = j + maxJ*i;
end
qq(i).t= a';
end
[ntimesteps , rMin, rMax ,ss ,ncs ,plotOn ,azimuthalSet ,azimuthalSetSize ,printStatus ,lags, blocLength, saveDir]=constants();

%%
%genStr=input("(L)oad or (G)enerate? If load eg c2t4\n> ","s");
genStr="G";
if genStr=="G" ||  genStr=="g"
  sprintf('%s%s%s%s%s%s%s%s','**************',newline,'Generating to file C',num2str(ncs),'t',num2str(ntimesteps),newline,'**************'   )
emptyStr=[];
%[qMinusQbar]= fftStep("readDataAndFindVeloFluctuation",emptyStr);
[qMinusQbar]= fftStep("readDataAndFindVeloFluctuation",emptyStr);

sprintf('%s','Saving velocity fluctuations into file...')
saveStr=['/mnt/archLv/mike/podData/structSave/qMinusQbar_CaseC' num2str(ncs) 't' num2str(ntimesteps) '.mat'       ];
%save(saveStr,'qMinusQbar','-v7.3');
sprintf('%s%s','Saved velocity fluctuations into file ',saveStr);
%% 
else % eg load qMinusQbar_C2T400
saveStr=['/mnt/archLv/mike/podData/structSave/qMinusQbar_' genStr '.mat'];
load(saveStr)
sprintf('%s%s','Loaded velocity fluctuations from file ',saveStr);
end % if
%% 
sprintf('%s','done qMinusQbar.')
% Next step is take azimuthal fft...
% [post]=fftStep("azimuth",qMinusQbar)
sprintf('%s','done qMinusQbar.')