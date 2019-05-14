% ------------------------------------------------------------------------------
%
%                           function eci2hill
%
%  this function finds the relative pos/vel vectors given the ECI target and 
%  interceptor vectors.
%
%  Routine not dependent on km or m for distance unit
%  all vectors are column vectors
%
%  author        : sal alfano              719-573-2600   26 may 2010
%

function [rhill,vhill] = eci2hill(rtgteci, vtgteci, rinteci, vinteci);

%  find rotation matrix from ECI to rsw frame
%  convert target and interceptor, compute vecotr magnitudes
        magrtgt = norm(rtgteci);
        magrint = norm(rinteci);
        [rtgtrsw,vtgtrsw,rotECI2RSW] = rv2rsw(rtgteci, vtgteci);        
        rintrsw  = matvecmult( rotECI2RSW, rinteci, 3);
        vintrsw  = matvecmult( rotECI2RSW, vinteci, 3);
    
%  find rotation angles (radians) to go from target to interecptor
        sinphiint    = rintrsw(3) / magrint;
        phiint       = asin(sinphiint);
        cosphiint    = cos(phiint);
        lambdaint    = atan2(rintrsw(2), rintrsw(1));
        sinlambdaint = sin(lambdaint);
        coslambdaint = cos(lambdaint);
        lambdadottgt = norm(vtgteci)/ magrtgt;  % if circular ==>> norm(vtgteci)/magrtgt; % if not ++>> vtgtrsw(2) / magrtgt
        
%  find position component positions
        rhill(1) = magrint - magrtgt;
        rhill(2) = lambdaint * magrtgt;
        rhill(3) = phiint * magrtgt;

%  find rotation matrix to go from rsw to SEZ of inerceptor
        rotrswtoSEZ(1,1) = sinphiint * coslambdaint;
        rotrswtoSEZ(1,2) = sinphiint * sinlambdaint;
        rotrswtoSEZ(1,3) = -cosphiint;
        rotrswtoSEZ(2,1) = -sinlambdaint;
        rotrswtoSEZ(2,2) = coslambdaint;
        rotrswtoSEZ(2,3) = 0.0;
        rotrswtoSEZ(3,1) = cosphiint * coslambdaint;
        rotrswtoSEZ(3,2) = cosphiint * sinlambdaint;
        rotrswtoSEZ(3,3) = sinphiint;
        
%  find velocity component positions by using angular rates in SEZ frame 
        vintSEZ      = matvecmult( rotrswtoSEZ, vintrsw, 3);
        phidotint    = -vintSEZ(1)/magrint;
        lambdadotint = vintSEZ(2)/(magrint * cosphiint);
        
        vhill(1) = vintSEZ(3);  % if circular ==>> vintSEZ(3); % if not ++>> vintSEZ(3) - vtgtrsw(1)
        vhill(2) = magrtgt * (lambdadotint - lambdadottgt);
        vhill(3) = magrtgt * phidotint;
        
        

      
      
      