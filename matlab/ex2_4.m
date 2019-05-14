%     -----------------------------------------------------------------
%
%                              Ex2_4.m
%
%  this file demonstrates example 2-4.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2004
%                            by david vallado
%
%     (h)               email dvallado@msn.com
%     (w) 719-573-2600, email dvallado@stk.com
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

    constastro;

    fprintf(1,'\n-------- kepler  ex 2-4, pg 102 --------- \n' );

    % initial vectors in km and km/s
    ro = [ 1131.340  -2282.343  6672.423];
    vo = [ -5.64305  4.30333  2.42879 ];
    fprintf(1,'input: \n' );
    fprintf(1,'ro %16.8f %16.8f %16.8f km \n',ro );
    fprintf(1,'vo %16.8f %16.8f %16.8f km/s \n',vo );

    % convert 40 minutes to seconds
    dtsec = 40.0*60.0;
    fprintf(1,'dt %16.8f sec \n',dtsec );
    fprintf(1,'intermediate values: \n' );

    [r1,v1] =  kepler ( ro,vo, dtsec );

    % answer in km and km/s
    fprintf(1,'output: \n' );
    fprintf(1,'r1 %16.8f %16.8f %16.8f er \n',r1/re );
    fprintf(1,'r1 %16.8f %16.8f %16.8f km \n',r1 );
    fprintf(1,'v1 %16.8f %16.8f %16.8f er/tu \n',v1/velkmps );
    fprintf(1,'v1 %16.8f %16.8f %16.8f km/s \n',v1 );


%     % initial coes
%     rad = 180.0/pi;
%     [ro, vo] = coe2rv (7358.39, 0.0, 28.5/rad, 0.0/rad, 30.0/rad, 0.0/rad, 0.0, 0.0, 0.0);
%     fprintf(1,'input: \n' );
%     fprintf(1,'ro %16.8f %16.8f %16.8f km \n',ro );
%     fprintf(1,'vo %16.8f %16.8f %16.8f km/s \n',vo );
% 
%     % convert 40 minutes to seconds
%     dtsec = 40.0*60.0;
%     fprintf(1,'dt %16.8f sec \n',dtsec );
%     fprintf(1,'intermediate values: \n' );
% 
%     [r1,v1] =  kepler ( ro,vo, dtsec );
% 
%     % answer in km and km/s
%     fprintf(1,'output: \n' );
%     fprintf(1,'r1 %16.8f %16.8f %16.8f er \n',r1/re );
%     fprintf(1,'r1 %16.8f %16.8f %16.8f km \n',r1 );
%     fprintf(1,'v1 %16.8f %16.8f %16.8f er/tu \n',v1/velkmps );
%     fprintf(1,'v1 %16.8f %16.8f %16.8f km/s \n',v1 );
% 

    % initial coes with more than one period = 6281.815597 sec
    rad = 180.0/pi;
    [ro, vo] = coe2rv (7358.39, 0.0, 28.5/rad, 0.0/rad, 30.0/rad, 0.0/rad, 0.0, 0.0, 0.0);
    fprintf(1,'input: \n' );
    fprintf(1,'ro %16.8f %16.8f %16.8f km \n',ro );
    fprintf(1,'vo %16.8f %16.8f %16.8f km/s \n',vo );

    % convert 40 minutes to seconds
    dtsec = 4000.0*60.0;
    dtsec = 1.291007302335531e+03;
    dtsec = 6281.815597;
    fprintf(1,'dt %16.8f sec \n',dtsec );
    fprintf(1,'intermediate values: \n' );

    [r1,v1] =  kepler ( ro,vo, dtsec );

    % answer in km and km/s
    fprintf(1,'output: \n' );
    fprintf(1,'r1 %16.8f %16.8f %16.8f er \n',r1/re );
    fprintf(1,'r1 %16.8f %16.8f %16.8f km \n',r1 );
    fprintf(1,'v1 %16.8f %16.8f %16.8f er/tu \n',v1/velkmps );
    fprintf(1,'v1 %16.8f %16.8f %16.8f km/s \n',v1 );

