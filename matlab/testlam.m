
fid = 1;
constastro;
constmath;
st = 'y';

[fid,outfile] = fopen('fig712.out', 'wt');
for kt= 5 : 5  % 1 : 4, 7 is to test against Ochoa paper
    % ---- use for fixed r and r locations, and fig 7-12, 16, and 20 ---- }
    %      readln( infile,ro(1),ro(2),ro(3),rtgt(1),rtgt(2),rtgt(3),dt,dmy,direc );
    %      fprintf( kt:3,'fig 7-12 ro rtgt ',ro(4):11:7,rtgt(4):11:7 );
    kepmov = 'n';
    if kt == 1
        rinto = [ 1.040000   0.000000  0.000000 ] * re;
        rtgto = [-1.040000  0.2000000  0.000000 ] * re;
        dt = 20*tumin;
        direc = 's';
    end;
    if kt == 2
        rinto = [ 2.000000   0.000000  0.000000 ] * re;
        rtgto = [-2.000000  0.2000000  0.000000 ] * re;
        dt = 20*tumin;
        direc = 's';
    end;
    if kt == 3
        rinto = [ 4.000000   0.000000  0.000000 ] * re;
        rtgto = [-4.000000  0.2000000  0.000000 ] * re;
        dt = 20*tumin;
        direc = 's';
    end;
    if kt == 4
        rinto = [ 6.610700   0.000000  0.000000 ] * re;
        rtgto = [-6.610700  0.2000000  0.000000 ] * re;
        dt = 20*tumin;
        direc = 's';
    end;
    % make up a velocity vector (circular orbit) for dv calcs...
    if kt > 0 & kt < 5
        vinto = [0.0  sqrt(mu/mag(rinto))  0.0];
        vtgto = [0.0  sqrt(mu/mag(rtgto))  0.0]; % close, but not exact
    end;
    if kt == 5  % fixed target test case fig 7-16
        rinto = [ -6518.1083  -2403.8479  -22.1722 ];
        vinto = [  2.604057  -7.105717 -0.263218];
        rtgto = [  6697.4756  1794.5831  0.00000 ];
        vtgto = [  -1.962372  7.323674  0.0000000 ];
        dt = 20*tumin;
        direc = 's';
        kepmov = 'n';
    end;
    if kt == 6   % moving target test case fig 7-20
        rinto = [ 5328.7862   4436.1273   101.4720 ];
        vinto = [  -4.864779  5.816486  0.240163];
        rtgto = [  6697.4756  1794.5831  0.00000 ];
        vtgto = [  -1.962372  7.323674  0.0000000 ];
        dt = 20*tumin;
        direc = 's';
        kepmov = 'y';
    end;
    if kt == 7  % fixed target test case papers...
        trangle = 90.0/rad;
        rinto = [ 1.0  0.0   0.0 ];
        vinto = [  2.604057  -7.105717 -0.263218];
        rtgto = [  1.000*cos(trangle) 1.000*sin(trangle)  0.0 ];
        vtgto = [  -1.962372  7.323674  0.0000000 ];
        dt = 20*tumin;
        direc = 's';
        kepmov = 'n';
        %    mu = 4.0*pi*pi;
    end;
    
    direc = 's';
    nrev = 0;
    dolam;
    
    direc = 'l';
    nrev = 0;
    dolam;
    
    direc = 's';
    nrev = 1;
    dolam;
    
    direc = 'l';
    nrev = 1;
    dolam;
    
    direc = 's';
    nrev = 2;
    dolam;
    
    direc = 'l';
    nrev = 2;
    dolam;
    
end; % for through fig 7-13 case }

fprintf( 'done with fig 7-13 data ' );
fclose (fid);







