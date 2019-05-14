%     -----------------------------------------------------------------
%
%                              Ex1_1.m
%
%  this file demonstrates example 1-1.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (h)               email dvallado@msn.com
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            30 mar 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************


    fprintf(1,'\n-------- ex 1-1 --------- \n' );

    constmath;
    constastro;
    periodsid = 86164.090518;
%    periodsid = 86400/1.002737909350795;
    
    a = (mu*(periodsid/twopi)^2)^(1.0/3.0);

    fprintf(1,'a %16.8f %18.10f km \n',a,a/re );

