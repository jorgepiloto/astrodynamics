%     -----------------------------------------------------------------
%
%                              Ex3_13.m
%
%  this file demonstrates example 3-13.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2013
%                            by david vallado
%
%     (h)               email davallado@gmail.com
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            13 feb 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************

       % -------- page 203
       % -------- 1995 June 8, 20:18:3.703691
       jd = 2449877.3458762;
       
       fprintf(1,'JD %8.7f  \n\n', jd );
       
       jdfrac = jd-floor(jd);
       jd = floor(jd);
       
       fprintf(1,'integer %6i  \n', jd );
       fprintf(1,'fractional  %3.7f  \n\n', jdfrac);       
       
       [year,mon,day,hr,min,sec] = invjday ( jd, jdfrac );

       

       fprintf(1,'year %6i  \n',year);
       fprintf(1,'mon  %3i  \n',mon);
       fprintf(1,'day  %3i  \n',day);
       fprintf(1,'hr   %3i  \n',hr);
       fprintf(1,'min  %3i  \n',min);
       fprintf(1,'sec  %3.6f  \n',sec);

