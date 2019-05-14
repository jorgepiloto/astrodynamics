%     -----------------------------------------------------------------
%
%                              ex6_1.m
%
%  this file demonstrates example 6-1.
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
%            25 nov 08  david vallado
%                         original
%  changes :
%            25 nov 08  david vallado
%                         original baseline
%
%     *****************************************************************

      rad = 180.0 / pi;
      re = 6378.137;  
      fprintf(1,'-------------------- problem ex 6-1 \n');
      rinit  = (re + 191.3411)/re;
      rfinal = (re + 35781.34857)/re;
      einit  = 0.0;
      efinal = 0.0;
      nuinit = 0.0/rad;
      nufinal= 180.0/rad;

      fprintf(1,' rinit  %11.7f %11.7f km \n',rinit, rinit*re);
      fprintf(1,' rfinal %11.7f %11.7f km \n',rfinal, rfinal*re);
      fprintf(1,' einit   %11.7f \n',einit);
      fprintf(1,' efinal  %11.7f \n',efinal);
      fprintf(1,' nuinit  %11.7f deg \n',nuinit * rad);
      fprintf(1,' nufinal %11.7f deg \n',nufinal * rad);

      [deltava,deltavb,dttu ] = hohmann (rinit,rfinal,einit,efinal,nuinit,nufinal);

      constastro;
      fprintf(1,'hohmann answers \n');
      fprintf(1,' deltava  %11.7f  %11.7f km/s \n',deltava, deltava*velkmps );
      fprintf(1,' deltavb  %11.7f  %11.7f km/s \n',deltavb, deltavb*velkmps );
      fprintf(1,' dttu  %11.7f tu %11.7f min \n',dttu,dttu*tumin);


      
      
      
      
