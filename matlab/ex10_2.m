%     -----------------------------------------------------------------
%
%                              Ex10_2.m
%
%  this file demonstrates example 10-2.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2007
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%             7 jun 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************

    % problem 2
    fprintf(1,'problem 2 -------------------------\n');
    xo=[1 2 3 4 5 6 7 8];
    yo=[1 1 2 3 3 4 7 6];
    fprintf(1,'xo %3i %3i %3i %3i %3i %3i %3i %3i \n',xo);
    fprintf(1,'yo %3i %3i %3i %3i %3i %3i %3i %3i \n',yo);

    ata=[8 sum(xo); sum(xo) sum(xo*xo')]
    atb=[sum(yo); sum(xo*yo')]
    atai = inv(ata)
    ans= atai*atb

    yc=ans(1)+ans(2)*xo;
    fprintf(1,'yc   %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',yc);
    res=yo-yc;
    fprintf(1,'res  %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',res);

    rms = sqrt(1/8*sum(res*res'));
    conv = sqrt(sum(res*res')/7);

    fprintf(1,'sum res 2 %11.7f  rms %11.7f  sum res /7 %11.7f \n',sum(res*res'),rms,conv);
    fprintf(1,'adjusted cov %11.7f  %11.7f \n',conv*sqrt(atai(1,1)),conv*sqrt(atai(2,2)) );


    fprintf(1,'problem 2 with bad meas so only 7 -------------------------\n');
    clear all;
    xo=[1 2 3 4 5 6 8];
    yo=[1 1 2 3 3 4 6];
    fprintf(1,'xo %3i %3i %3i %3i %3i %3i %3i  \n',xo);
    fprintf(1,'yo %3i %3i %3i %3i %3i %3i %3i  \n',yo);

    ata=[7 sum(xo); sum(xo) sum(xo*xo')]
    atb=[sum(yo); sum(xo*yo')]
    atai = inv(ata)
    ans= atai*atb

    yc=ans(1)+ans(2)*xo;
    res=yo-yc;
    rms = sqrt(1/7*sum(res*res'))
    conv = sqrt(sum(res*res')/6);
    fprintf(1,'sum res 2 %11.7f  rms %11.7f  sum res /7 %11.7f \n',sum(res*res'),rms,conv);
        fprintf(1,'adjusted cov %11.7f  %11.7f \n',conv*sqrt(atai(1,1)),conv*sqrt(atai(2,2)) );
 
        
    % problem 3 example with weighting
    fprintf(1,'problem 3 with weight added in, original 8 obs -------------------------\n');
    xo=[1 2 3 4 5 6 7 8]';
    yo=[1 1 2 3 3 4 4 6]';
    fprintf(1,'xo %3i %3i %3i %3i %3i %3i %3i %3i \n',xo);
    fprintf(1,'yo %3i %3i %3i %3i %3i %3i %3i %3i \n',yo);

    w1 = 1.0/ 0.1; % 0.1
    w2 = 1.0/ 0.02; % 0.02
    w(1, 1) = w1;  % assume 10 5 and 2 % accruacy for the 2 sensors
    w(2,2) = w2;
    w(3,3) = w1;
    w(4,4) = w2;
    w(5,5) = w1;
    w(6,6) = w2;
    w(7,7) = w1;
    w(8,8) = w2;
    
    i = [1;1;1;1;1;1;1;1];
    a = [i xo];
    b = [yo];
    atw = a' * w
    atwa = a' * w * a
    atwb = a' * w * b
    atwai = inv(atwa)
    ans= atwai*atwb
    
%    ata=[8 sum(xo); sum(xo) sum(xo*xo')]
%    atb=[sum(yo); sum(xo*yo')]
%    atai = inv(ata)
%    ans= atai*atb

    yc=ans(1)+ans(2)*xo;
    fprintf(1,'yc   %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',yc);
    res=yo-yc;
    fprintf(1,'res  %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n',res);

    rms = sqrt(1/8*sum(res'*res));
    conv = sqrt(sum(res'*res)/7);

    fprintf(1,'sum res 2 %11.7f  rms %11.7f  sum res /7 %11.7f \n',sum(res'*res),rms,conv);
    fprintf(1,'adjusted cov %11.7f  %11.7f \n',conv*sqrt(atai(1,1)),conv*sqrt(atai(2,2)) );

        
