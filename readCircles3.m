% write matlb script that imports into multi-d array..
function [myPreFft_noCsYet]=readCircles3(t,c) % redo for parfor loop.
                                                     % % does not take such a big struct.
[~, ~ ,rMax, ~ ,~ ,~ ,~, ~ ,~, ~]=constants();
  tic
  [myPreFft_noCsYet]=initData2("myPreFft_noCsYet");

      sprintf('%s%d%s%d','*t = ',t,'**c=',c)
      M_mat = zeros(540,540); % streamwiseData
      cc = sprintf( '%03d', c  ) ;
      tt = sprintf( '%04d', t  ) ;
      daPath = '/home/mi/podData/snaps2/';
      fileName = [ daPath 'snap_cs_' num2str(cc) '_ts_' num2str(tt) '.dat']; % t starts at 0.
      formatSpec = '%f';
      a=fopen(fileName, 'r');
      importedData = fscanf(a,formatSpec);
      %for r=1:rMax  % r is azimuthal mode
      %for p=1:540 % p is point
      %  pointOnLine_r = p + 540*r;
      %  M_mat(r,p)=importedData(pointOnLine_r,1); % 5 for streamwise entry of the data matrix % 1 for 1 col
      %end %p
      %end %r

      for r=0:rMax-1  % r is azimuthal mode
      for p=1:540 % p is point
        pointOnLine_r = p + 540*r;
        M_mat(r+1,p)=importedData(pointOnLine_r,1); % 5 for streamwise entry of the data matrix % 1 for 1 col
      end %p
      end %r

%&---------------------------------------------------------------------
      for p=1:540 % p is point
        myPreFft_noCsYet(p).dat =    M_mat(:,p);
      end %p
end % fc
