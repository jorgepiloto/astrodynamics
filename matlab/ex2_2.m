%     -----------------------------------------------------------------
%
%                              Ex2_2.m
%
%  this file demonstrates example 2-2.
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
        constmath;

        % --------  newtone      - find true and mean anomaly given ecc and eccentric anomaly
        ecc= 0.4;
        e0 = 334.566986/rad;
        fprintf(1,'              m             nu           e           e0  \n');
        [m,nu]= newtone(ecc,e0);
        fprintf(1,'newe  %14.4f%14.4f%14.4f%14.7f\n',m*rad,nu*rad,ecc,e0*rad );

        % --------  newtonm      - find eccentric and true anomaly given ecc and mean anomaly
        ecc= 0.34;
        m = 235.4/rad;
        fprintf(1,'             e0             nu           e           m   \n');
        [e0,nu]= newtonm(ecc,m);
        fprintf(1,'newm  %14.4f%14.4f%14.4f%14.7f\n',e0*rad,nu*rad,ecc,m*rad );

        % --------  newtonnu     - find eccentric and mean anomaly given ecc and true anomaly
        ecc= 0.34;
        nu = 134.567001/rad;
        fprintf(1,'             e0             m            e           nu  \n');
        [e0,m]= newtonnu(ecc,nu);
        fprintf(1,'newnu %14.4f%14.4f%14.4f%14.7f\n',e0*rad,m*rad,ecc,nu*rad );


        
        p = 25512;
        dt = 53.7874 * 60; % in sec
        ecc = 1.0;
        np = 2.0 * sqrt(398600.4418/(p^3));
        cot2s = 1.5*np*dt;
        s = 0.5 * acot(cot2s);
        w= atan( (tan(s))^(1/3) );
        b = 2.0 * cot(2.0*w);
        
        fprintf(1,'np %11.7f cot2s %11.7f s %11.6f w %11.6f b %11.6f  \n',np, cot2s, s*rad, w*rad, b );
        
        
        
        
