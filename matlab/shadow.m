    %
    % calculate algorithm 34 quantities
    %

     rs = 696000.0;
     re = 6378.1363;
     au = 149597870.0;

     angumb = atan((rs-re)/au);
     angpen = atan((rs+re)/au);

     reci1 = [-41221.79149309;      8864.59854079;     0.00000000];
     veci1 = [-0.646416796;    -3.005940793;    -0.000000000];

     % +50 and -80 seem to work here to get the proper angles
     dtsec = input('input dtsec to propagate ');      

     [reci,veci,error] =  kepler  ( reci1,veci1, dtsec );

     year = 2008;
     mon =  3;
     day = 16;
     hr =  6;
     min = 13;
     sec = 0.00;
     [jd,jdfrac] = jday( year,mon,day,hr,min,sec);

     [rsun,rtasc,decl] = sun ( jd+jdfrac );

     dot(reci,rsun);
     umbvert = 0.0;
     penvert = 0.0;
     umb = 'n';
     pen = 'n';

     % now for algorihtm 34
     if dot(reci,rsun) < 0.0
         ang1 = angl(-rsun, reci);

         sathoriz = mag(reci )*cos(ang1);
         satvert  = mag(reci )*sin(ang1);
         x = re/sin(angpen);
         penvert = tan(angpen)*(x + sathoriz);
         if satvert <= penvert
             pen = 'y';
             y = re/sin(angumb);
             umbvert = tan(angumb)*(y-sathoriz);
             if satvert <= umbvert
                 umb = 'y';
             end;
         end;
     end;

     fprintf(1,' %11.7f  %11.4f  %11.4f  %11.4f  %11.4f U %c  P %c \n',ang1*180.0/pi, sathoriz, satvert, penvert, umbvert,umb,pen );
     s = 2.0 * mag(reci) * ang1;
     fprintf(1,' %11.7f  %11.7f \n',s, (s/mag(veci))/60.0);