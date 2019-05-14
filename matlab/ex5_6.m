%     -----------------------------------------------------------------
%
%                              Ex5_6.m
%
%  this file demonstrates example 5-6.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            26 mar 11  david vallado
%                         original
%  changes :
%            26 mar 11  david vallado
%                         original baseline
%
%     *****************************************************************

       constmath;

       [jd,jdfrac] = jday( 1995, 2, 15, 12, 0, 0.00 );
       r1 = [0.0 -4464.696 -5102.509 ]; % km
       r2 = [0.0 5740.323 3189.068];

       [los] = sight ( r1, r2, 's' );
       los

       [jd,jdfrac] = jday( 1995, 2, 15, 0, 0, 0.00 );
       [rsun,rtasc,decl] = sun ( jd+jdfrac );
       fprintf(1,'sun MOD %11.9f%11.9f%11.9f au\n',rsun );
       fprintf(1,'sun MOD %14.4f%14.4f%14.4f km\n',rsun*149597870.0 );

       [los] = sight ( r1, rsun*149597870.0, 's' );
       los
       

