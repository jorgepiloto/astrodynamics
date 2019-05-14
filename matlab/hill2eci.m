% ------------------------------------------------------------------------------
%
%                           function hill2eci
%
%  this function finds the interceptor's ECI pos/vel vectors 
%  given the ECI target and hill (relative) interceptor vectors.
%
%  Routine not dependent on km or m for distance unit, but needs to be consistent!
%  all vectors are column vectors
%
%  author        : sal alfano         719-573-2600   13 aug 2010
%

function [rinteci, vinteci] = hill2eci(rtgteci, vtgteci, rinthill, vinthill);

%  find rotation matrix from ECI to rsw frame
%  convert target and interceptor, compute vecotr magnitudes
        magrtgt = norm(rtgteci);
        magrint = magrtgt + rinthill(1);
        [rtgtrsw,vtgtrsw,rotECI2RSW] = rv2rsw(rtgteci, vtgteci);        
 
%  find rotation angles (radians) to go from target to interecptor
        lambdadottgt = norm(vtgteci)/ magrtgt;  % if circular ==>> norm(vtgteci)/magrtgt; % if not ++>> vtgtrsw(2) / magrtgt
        lambdaint    = rinthill(2)/magrtgt;
        phiint       = rinthill(3)/magrtgt;
        sinphiint    = sin(phiint);
        cosphiint    = cos(phiint);
        sinlambdaint = sin(lambdaint);
        coslambdaint = cos(lambdaint);
        
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
        rdotint      = vinthill(1); % if circular ==>> vinthill(1);  % if not ++>> vinthill(1) + vtgtrsw(1)       
        lambdadotint = vinthill(2)/magrtgt + lambdadottgt;
        phidotint    = vinthill(3)/magrtgt;
        vintSEZ(1) = -magrint * phidotint;
        vintSEZ(2) = magrint * lambdadotint * cosphiint;
        vintSEZ(3) = rdotint;                
        vintrsw    = rotrswtoSEZ' * vintSEZ';
        vinteci    = rotECI2RSW' * vintrsw;
        
%  find position component positions
        rintrsw(1) = cosphiint * magrint * coslambdaint;
        rintrsw(2) = cosphiint * magrint * sinlambdaint;
        rintrsw(3) = sinphiint * magrint;
        rinteci    = rotECI2RSW' * rintrsw';

        
        

      
      
      