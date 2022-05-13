function [loadedBloc]=loadBloc(t,c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step A) load a chonk into memory and read circles in.
    for timeBloc=1:blocLength
    parfor t = 1:ntimesteps % time % <-- nb, this is the parfor loop.
    myPreFft_noCsNoTimeYet=readCircles2(timeBloc*t,c);
    myPreFft_noCsYet(t).circle=myPreFft_noCsNoTimeYet;
    sprintf('%s','pause')
    end % parfor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % fc