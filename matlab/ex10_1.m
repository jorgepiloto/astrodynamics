%     -----------------------------------------------------------------
%
%                              Ex10_1.m
%
%  this file demonstrates example 10-1.
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

    % problem 1
    fprintf(1,'problem 1 --------------------------\n');
    xo=[1 2 3 4 5 6 7 8];
    yo=[1 1 2 3 3 4 4 6];
    fprintf(1,'xo %3i %3i %3i %3i %3i %3i %3i %3i \n',xo);
    fprintf(1,'yo %3i %3i %3i %3i %3i %3i %3i %3i \n',yo);

    ata=[8 sum(xo); sum(xo) sum(xo*xo')]
    atb=[sum(yo); sum(xo*yo')]
    atai = inv(ata)
    ans= atai*atb

