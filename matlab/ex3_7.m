%     -----------------------------------------------------------------
%
%                              Ex3_7.m
%
%  this file demonstrates example 3-7.
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
%            13 feb 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************


        year = 2004;
        mon  =   5;
        day  =  14;
        hr   =  10;
        min  =  43;
        sec  =   0.0;
        dut1 = -0.463326;
        dat  = 32;
        xp   =  0.0;
        yp   =  0.0;
        lod  =  0.0;
        timezone= 6;

        % -------- convtime    - convert time from utc to all the others
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
           = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

        fprintf(1,'ut1 %8.6f tut1 %16.12f jdut1 %18.11f\n',ut1,tut1,jdut1 );
        fprintf(1,'utc %8.6f\n',utc );
        fprintf(1,'tai %8.6f\n',tai );
        fprintf(1,'tt  %8.6f ttt  %16.12f jdtt  %18.11f\n',tt,ttt,jdtt );
        fprintf(1,'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n',tdb,ttdb,jdtdb );




