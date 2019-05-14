%     -----------------------------------------------------------------
%
%                              Ex3_4.m
%
%  this file demonstrates example 3-4.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2007
%                            by david vallado
%
%     (h)               email dvallado@msn.com
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

       % -------- jday         - find julian date
       year = 1996;
       mon = 10;
       day = 26;
       hr  = 14;
       minute = 20;
       secs = 0.00;
       fprintf(1,'\n--------jday test \n' );
       fprintf(1,'year %4i ',year);
       fprintf(1,'mon %4i ',mon);
       fprintf(1,'day %3i ',day);
       fprintf(1,'hr %3i:%2i:%8.6f\n ',hr,minute,secs );

       [jd, jdfrac]= jday(year,mon,day,hr,minute,secs);
       fprintf(1,'jd %18.10f  %18.10f   %18.10f \n',jd, jdfrac, jd + jdfrac);

       [year,mon,day,hr,minute,secs] = invjday ( 2450382.5, jdfrac );
       fprintf(1,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);
     
       [year,mon,day,hr,minute,secs] = invjday ( 2450382.5, jdfrac - 0.2 );
       fprintf(1,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);

       [year,mon,day,hr,minute,secs] = invjday ( 2450382.5 + 1, jdfrac + 1.5 );
       fprintf(1,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);
     
       [year,mon,day,hr,minute,secs] = invjday ( 2450382.5, -0.5 );
       fprintf(1,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);

       [year,mon,day,hr,minute,secs] = invjday ( 2450382.5, +0.5 );
       fprintf(1,'year %5i   mon %4i day %3i %3i:%2i:%8.6f\n',year, mon, day, hr, minute, secs);
