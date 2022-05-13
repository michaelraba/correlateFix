function FinalP2()
    % This is Melborne Data
    sprintf('%s','hi');
    path='/home/mi/Downloads/turbulence_finalProjectData/Data/Boundary\ Layer\ Data\ 1/';
    pathM='/home/mi/Downloads/turbulence_finalProjectData/Data/Boundary Layer Data 1/';
    lVec = [ 64, 64, 94, 58 , 60 ]; % length of file minus 28, which is the length # config lines
    for jj=1:5
    fileN=[path 'Case' int2str(jj) '.txt'];
    fileN_M=[pathM 'Case' int2str(jj) '.txt'];
    %fileName='/home/mi/Downloads/turbulence_finalProjectData/Data/Boundary Layer Data 1/c1';

% define constants
    kappa  = 0.41 % von karman cte
    u_tau = 0.416152 ;  % friction velocity
    nu = 1.51E-05 ;     % kinematic viscosity
    u_tauVec= [0.444139 , 0.416152 , 0.357878 , 0.783982 , 1.22388 ] ;
    nuVec = [ 1.51644e-005 , 1.51E-05 , 1.51E-05 , 1.54E-05 ,  1.58E-05];
    % y^+ = u_tau/nu*y.
%
     
    nFile = [path 'e1'];
    nFileM = [pathM 'e1'];
    %cmd=['tail -n65 "' fileN_M '" > ' nFile]; % 92-28
    cmd=['tail -n' num2str(lVec(jj)) ' "' fileN_M '" > ' nFile]; % 92-28


    system(cmd);
    formatSpec = '%lf%lf';
    a=fopen(nFileM, 'r');
    v=[];
    c = fscanf(a,formatSpec);
    v =[c,v];
    sprintf('%s','hi')
    sprintf('%s','hi')
    %y_mm =  zeros(length(v)/2,1);
    y_m =  zeros(length(v)/2,1);
    y_mLog =  zeros(length(v)/2,1);
    U_ms = zeros(length(v)/2,1);
    b = zeros(length(v)/2,1);
    eCount=1;
    oCount=1;
    for ii = 1:length(v)
      if mod(ii,2) == 0
        U_ms(eCount)=v(ii,1)/u_tauVec(jj); % velocity
        %U_ms(eCount)=v(ii,1); % velocity

        eCount=eCount+1;
      elseif mod(ii,2) ~= 0 % distance from wall
        %y_m(eCount)=1e-2*v(ii,1);
        % define y^+ , and take the log for the semilog plot.
        y_mLog(eCount)=  u_tauVec(jj)/nuVec(jj)*(1e-2)*v(ii,1); % for semi-log plot

    % y^+ = u_tau/nu*y.
        oCount=oCount+1;
      end %if
    end % for
    sprintf('%s','hi')
    % law_of_the_wall_approximation
    for ii = 1:length(y_mLog)
      b(ii) = 1/kappa*log(y_mLog(ii))+0;
    end % for
% plot
    labelStr=['Case ' num2str(jj) ];
    semilogx(y_mLog,U_ms,"DisplayName",labelStr);
    hold on;
 
% end plot
    end % jj case loop      
    plot(y_mLog,b,'b.',"DisplayName","law of the log");
    xline(250, 'r.' ,"DisplayName", "lower bound of log region")
    xline(10500, 'r.', "DisplayName", "upper bound of log region")

    legend();
end % function finalP1.m

