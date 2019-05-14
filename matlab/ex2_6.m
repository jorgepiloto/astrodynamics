%     -----------------------------------------------------------------
%
%                              Ex2_6.m
%
%  this file demonstrates example 2-6.
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

            constmath;
            constastro;

            fprintf(1,'coe test ----------------------------\n' );
            p     = 11067.790; % km
            ecc   = 0.83285;
            incl  = 87.87/rad;
            omega = 227.89/rad;
            argp  = 53.38/rad;
            nu    = 92.335/rad;
            arglat = 0.0;
            truelon= 0.0;
            lonper = 0.0;

            a = p/(1-ecc*ecc);

            fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
            fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n',...
                    p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad, ...
                    arglat*rad,truelon*rad,lonper*rad );

            % --------  coe2rv       - classical elements to posisiotn and velocity
            [r,v] = coe2rvS(p,ecc,incl,omega,argp,nu,arglat,truelon,lonper);
            fprintf(1,'r    %15.9f%15.9f%15.9f',r );
            fprintf(1,' v %15.10f%15.10f%15.10f\n',v );

