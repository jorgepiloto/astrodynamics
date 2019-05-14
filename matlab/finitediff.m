% -----------------------------------------------------------------------------
%
%                           procedure finitediff
%
% this procedure perturbs the components of the state vector for processing
% with the finite differencing for the a matrix.
%
%  author        : david vallado                  719-573-2600   15 jan 2008
%
%  inputs          description                    range / units
%    whichconst  - parameter for sgp4 constants   wgs72, wgs721, wgs84
%    pertelem    - which element to perturb
%    percentchg  - amount to modify the vectors   0.001
%                  by in finite differencing
%    deltaamtchg - tolerance for small value in
%                  finite differencing            0.0000001
%    statetype   - type of elements (equinoctial, etc)  'e', 't'
%    xnom        - state vector                   varied
%    scalef      - scale factor for state         all set to 1.0 now
%
%  outputs       :
%    deltaamt    - amount each elemnt is perturbed
%    satrec      - satellite record
%
%  locals        :
%    jj          - index
%
%  coupling      :
%    getgravconst- get the constants for a gravity field for sgp4
%    state2satrec- conversion between state and satellite structure
%    sgp4init    - intiialize sgp4 constants
%
%  references    :
%    vallado       2007, 753-765
% --------------------------------------------------------------------------- */

function [deltaamt, xnomp] = finitediff(pertelem, percentchg, deltaamtchg, xnom);
      deg2rad  =   pi / 180.0;         %   0.0174532925199433
    %  getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

      % chk if perturbing amt is too small. if so, up the percentchg and try again
      % this will execute 5 times, but leaves percentchg the same after each run
      jj = 1;
      deltaamt = 0.0;
      xnomp = xnom; 
    
      while ( (abs(deltaamt) < deltaamtchg) && (jj < 5) );
         if (jj > 1) 
             fprintf(1,'too large\n');  
         end   
          
          deltaamt = xnom(pertelem) * percentchg;
          xnomp(pertelem, 1) = xnom(pertelem, 1) + deltaamt;

%          state2satrec( xnom, scalef, statetype, statesize, eTo, satrec );

          if (abs(deltaamt) < deltaamtchg)       % 0.00001
              percentchg = percentchg * 1.4;  % increase by 40% and try again
              fprintf(1,' %i percentchg chgd %11.8f \n',jj,percentchg);
          end
          jj = jj + 1;
      end

      % printf(" \n");
end  % procedure finitediff

        
