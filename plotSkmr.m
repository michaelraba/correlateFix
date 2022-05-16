function plotSkmr(plotObject,isGraph)
[ntimesteps rMin ,rMax, ss ,ncs ,plotOn ,azimuthalSet, azimuthalSetSize ,printStatus, lags]=constants();
if isGraph=="graph"
  f=figure('Renderer', 'painters', 'Position', [10 10 1900 900])
  
  A=linspace(0,1,541)
  %A=linspace(0,1,1079)
  xlabel('radius $\frac{r}{R}$','interpreter','latex')
  ylabel("$S_{ii}(k,m;r,r')$",'interpreter','latex')  
  %xlabel('radius r');
  %ylabel('S');
  %figure;
  hold on;
  cou = 1
  for c=1:5
      subplot(3,1,cou);

      %for t=1:ntimesteps
      for m=1:18
      
      labelStr = ['(m,k)=(', num2str(azimuthalSet(m)),',',num2str(c),')'];
      pp=plot(A,real(plotObject(c).circle(m).dat((end-1)/2:end) )/ntimesteps,"DisplayName",labelStr);
      tiSt=['Spanwise mode: ' num2str(c)  ];
      title(tiSt, 'FontName','capitana','FontSize',12,'interpreter','latex')
      legend();
      %avgTimeEnd(1).circle(1).dat  
      %pp=plot(A,real(plotObject(c).circle(m).dat /ntimesteps),"DisplayName",labelStr)

      
      
      %pp=plot(A,real(plotObject(c).t(t).azimuthal(m).dat((end-1)/2:end) )/ntimesteps,"DisplayName",labelStr)
                    %SrrPrForFourierPost(c).t(t).azimuthal(m).dat
      %if mod(c,4)==0
      %  pp.LineStyle='-';
      %elseif mod(c,4)==1
      %  pp.LineStyle='--';
      %elseif mod(c,4)==2
      %  pp.LineStyle=':';
      %elseif mod(c,4)==3
      %  pp.LineStyle='-.';
      %end
      hold on;
      %pause(1);
      end % m
      cou = cou + 1;
  end %c
     % end %t
      %legend( );
      %saveas(gcf,'Sij.png','res')
      %sprintf('%s','pause');
elseif isGraph=="graphPause"
  hold on;
    for m=2:38
    for c=1:ncs
    plot(real(plotObject(m).crosssection(c).dat(end/2:end) )/ntimesteps)
    hold on;
    pause(1);
    end
    end
    %sprintf('%s','pause');
      pause(1)
else
end % if
  sgtitle("Time-averaged Correlation $S_{ii}(k,m;r,r')$",'FontName','capitana','FontSize',12,'interpreter','latex')

end % f
