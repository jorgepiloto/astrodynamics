%     -----------------------------------------------------------------
%
%                              Ex10_3.m
%
%  this file demonstrates example 10-3.
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

    % problem 3
    fprintf(1,'problem 3 -------------------------\n');
    xo=[1 2 3 4];
    yo=[2.5 8.0 19.0 50.0];
    fprintf(1,'xo %8i %8i %8i %8i \n',xo);
    fprintf(1,'yo %8.6f %8.6f %8.6f %8.6f \n',yo);

    beta = log(yo(4)/yo(3))/log(4/3)
    alpha=yo(3)/3^beta

    % do first time manually to get intermdeiate values
    for i = 1:4
        yn(i) = alpha*xo(i)^beta;
        parynalp(i) = xo(i)^beta;
        parynbet(i) = alpha*log(xo(i))*xo(i)^beta;
    end
    yn
    parynalp
    parynbet

    ata=[parynalp*parynalp' parynalp*parynbet'; parynalp*parynbet' parynbet*parynbet']
    atb=[parynalp*(yo-yn)'; parynbet*(yo-yn)']
    atai = inv(ata)
    ans= atai*atb;

    alpha = alpha + ans(1);
    beta = beta + ans(2);

    for j = 1 :4
        b(j) = yo(j) - alpha*(xo(j)^beta);
    end
    b
    rms = sqrt(b*b'/4);
    rmsold = rms;    
     fprintf(1,'0 dx %11.7f  %11.7f ans %11.7f  %11.7f rms %11.7f \n',ans(1), ans(2), alpha, beta, rms);


    for loop = 1 :3
        for i = 1:4
             yn(i) = alpha*xo(i)^beta;
            parynalp(i) = xo(i)^beta;
            parynbet(i) = alpha*log(xo(i))*xo(i)^beta;
        end

        ata=[parynalp*parynalp' parynalp*parynbet'; parynalp*parynbet' parynbet*parynbet'];
        atb=[parynalp*(yo-yn)'; parynbet*(yo-yn)'];
        atai = inv(ata);
        ans= atai*atb;

        alpha = alpha + ans(1);
        beta = beta + ans(2);

        for j = 1 :4
            b(j) = yo(j) - alpha*(xo(j)^beta);
        end    
        b
        rms = sqrt(b*b'/4);
        rmsdel = (rmsold-rms)/rmsold;
         fprintf(1,'%2i dx %11.7f  %11.7f ans %11.7f  %11.7f rms %11.7f %11.7f \n',loop, ans(1), ans(2), alpha, beta, rms,rmsdel);
         rmsold=rms;
    end  








