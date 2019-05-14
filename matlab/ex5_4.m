%     -----------------------------------------------------------------
%
%                              Ex5_4.m
%
%  this file demonstrates example 5-4.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            22 jan 11  david vallado
%                         original
%  changes :
%            22 jan 11  david vallado
%                         original baseline
%
%     *****************************************************************

constmath;

% --------  moon         - moon rise set
[jd,jdfrac] = jday( 1998, 8, 21, 0, 0, 0.00 );
latgd = 40.0/rad;
lon = 0.00 / rad;

[utmoonrise,utmoonset,moonphaseang,error] = moonrise( jd+jdfrac,latgd,lon,'y' )
%fprintf(1,'moon moonrise %14.4f    moonset %14.4f hrs \n',utmoonrise,utmoonset );
fprintf(1,'moon moonrise %14.4f  %14.4f   moonset %14.4f %14.4f  \n',utmoonrise, (utmoonrise-floor(utmoonrise))*60, utmoonset, (utmoonset-floor(utmoonset))*60 );
fprintf(1,'moon phase angle %14.4f   \n',moonphaseang );


[utmoonrise,utmoonset,moonphaseang,error] = moonrise3( jd+jdfrac,latgd,lon,'y' )
fprintf(1,'2moon moonrise %14.4f  %14.4f   moonset %14.4f %14.4f  \n',utmoonrise, (utmoonrise-floor(utmoonrise))*60, utmoonset, (utmoonset-floor(utmoonset))*60 );
fprintf(1,'2moon phase angle %14.4f   \n',moonphaseang );


[jd,jdfrac] = jday( 1990, 3, 5, 0, 0, 0.00 );
latgd =  40.94 /rad;
lon   = -73.97 / rad;

%        [utmoonrise,utmoonset,moonphaseang,error] = moonrise1( jd,latgd,lon,'y' )
%        fprintf(1,'moon moonrise %14.4f  %14.4f   moonset %14.4f %14.4f  \n',utmoonrise, (utmoonrise-floor(utmoonrise))*60, utmoonset, (utmoonset-floor(utmoonset))*60 );
%        fprintf(1,'moon phase angle %14.4f   \n',moonphaseang );

[utmoonrise,utmoonset,moonphaseang,error] = moonrise2( jd+jdfrac,latgd,lon,'y' )
fprintf(1,'moon moonrise %14.4f  %14.4f   moonset %14.4f %14.4f  \n',utmoonrise, (utmoonrise-floor(utmoonrise))*60, utmoonset, (utmoonset-floor(utmoonset))*60 );
fprintf(1,'moon phase angle %14.4f   \n',moonphaseang );



[jd,jdfrac] = jday( 2006, 6, 28, 0, 0, 0.00 );
latgd = 40.0/rad;
lon = 0.00 / rad;

[utmoonrise,utmoonset,moonphaseang,error] = moonrise2( jd+jdfrac,latgd,lon,'y' )
fprintf(1,'moon moonrise %14.4f    moonset %14.4f hrs \n',utmoonrise,utmoonset );
fprintf(1,'moon phase angle %14.4f   \n',moonphaseang );

pause;

fprintf(1,'     40    42    44    46    48    50    52    54    56    58    60    62    64    66  \n' );
for i = 8:30  %8:30
    [jd,jdfrac] = jday( 2006, 6, i, 0, 0, 0.00 );
    for j = 0:13  %0:13
        latgd = (40.0 + j * 2.0)/rad;
        lon = 0.00 / rad;
        
        [utmoonrise,utmoonset,moonphaseang,error] = moonrise3( jd+jdfrac,latgd,lon,'n' );
        %                if strcmp(error,'ok') == 0 % 1 if true, 0 if false
        %                    fprintf(1,'error');
        %                end;
        
        [hr,min,sec] = rad2hms( utmoonrise*15.0 * pi/180.0 ); % be sure to take hrs and convert to rad befor separating
        [hr1,min1,sec1] = rad2hms( utmoonset*15.0 * pi/180.0 ); % be sure to take hrs and convert to rad befor separating
        %  print out header date for each section of results
        if j == 0
            fprintf(1,'%2i ',i );
        end;
        
%         if utmoonrise > 9998.0 utmoonset > 9998.0 
%              fprintf(1,' none ');
%         else
%         if utmoonrise > 9998.0 && utmoonset < 24.0 
%              fprintf(1,' nors ');
%         else 
% %        if utmoonset - utmoonrise > 14.0 
% %             fprintf(1,' none ');
% %        else
%             fprintf(1,'%2i:%2i ',hr, min);            
% %        end;
%         end;
%         end;
%            
%           if hr >= 24
                [jd,jdfrac] = jday( 2006, 6, i, 0, 0, 0.00 );
                [el1] = moonel( jd+jdfrac,latgd,lon );
                [jd,jdfrac] = jday( 2006, 6, i+1, 0, 0, 0.00 );
                [el2] = moonel( jd+jdfrac,latgd,lon );

            if hr >= 24
                fprintf(1,'| nost  ');
            else  
                fprintf(1,'| %2i:%2i ',hr, min);
            end    
         end; % if j = 0
        fprintf(1,'\n');
    end;
        
    fprintf(1,'     40    42    44    46    48    50    52    54    56    58    60    62    64    66  \n' );
for i = 8:30  %8:30
    [jd,jdfrac] = jday( 2006, 6, i, 0, 0, 0.00 );
    for j = 0:13  %0:13
        latgd = (40.0 + j * 2.0)/rad;
        lon = 0.00 / rad;
        
        [utmoonrise,utmoonset,moonphaseang,error] = moonrise3( jd+jdfrac,latgd,lon,'n' );
        %                if strcmp(error,'ok') == 0 % 1 if true, 0 if false
        %                    fprintf(1,'error');
        %                end;
        
        [hr,min,sec] = rad2hms( utmoonrise*15.0 * pi/180.0 ); % be sure to take hrs and convert to rad befor separating
        [hr1,min1,sec1] = rad2hms( utmoonset*15.0 * pi/180.0 ); % be sure to take hrs and convert to rad befor separating
        %  print out header date for each section of results
        if j == 0
            fprintf(1,'%2i ',i );
        end;
        
                [jd,jdfrac] = jday( 2006, 6, i, 0, 0, 0.00 );
                [el1] = moonel( jd+jdfrac,latgd,lon );
                [jd,jdfrac] = jday( 2006, 6, i+1, 0, 0, 0.00 );
                [el2] = moonel( jd+jdfrac,latgd,lon );

            if hr1 >= 24
                fprintf(1,'| nost  ');
            else  
                fprintf(1,'| %2i:%2i ',hr1, min1);
            end    
         end; % if j = 0
        fprintf(1,'\n');
    end;
    
    
    