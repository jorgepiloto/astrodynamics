% -----
%   vallado       2007, 370, alg 46, ex 6-10
% function [ttrans,tphase,dvphase,dvtrans1,dvtrans2,aphase ] = noncoplr(phasei,aint,atgt,ktgt,kint,arglatint,nodeint,truelon,deltai);
%------ }

function [ttrans,tphase,dvphase,dvtrans1,dvtrans2,aphase ] = noncoplr(phasei,aint,atgt,ktgt,kint,arglatint,nodeint,truelon,deltai);
     twopi =  6.28318530717959;
     rad   = 57.29577951308230;
     mu = 1.0;  % cannonical

     angvelint = sqrt( mu / (aint * aint * aint) );
     angveltgt = sqrt( mu / (atgt * atgt * atgt) );
     atrans   = (aint + atgt) * 0.5;
     ttrans = pi * sqrt( (atrans * atrans * atrans) / mu );

     deltatnode = phasei / angvelint;

     lead = angveltgt * ttrans;

     omeganode = angveltgt * deltatnode;

     phasenew = nodeint + pi - (truelon + omeganode);
 
     leadnew = pi + phasenew;

     tphase= (leadnew - lead + twopi * ktgt) / angveltgt;

     aphase = (mu * (tphase/(twopi * kint)) ^ 2 ) ^(1.0/3.0);

     % -----------------  find deltav's  ----------------- }
     vint= sqrt(mu/aint);
     vphase= sqrt(2.0 * mu/aint - mu/aphase);
     dvphase= vphase - vint;

     vtrans1= sqrt(2.0 * mu/aint - mu/atrans);
     dvtrans1= vtrans1 - vphase;

     vtrans2= sqrt(2.0 * mu/atgt - mu/atrans);
     vtgt= sqrt(mu/atgt);
     dvtrans2= sqrt(vtgt * vtgt + vtrans2 * vtrans2 - 2.0 * vtgt * vtrans2 * cos(deltai));

     dvtotal = dvphase + dvtrans1 + dvtrans2;
     ttotal = deltatnode + ttrans + tphase;
     
     constastro;
     fprintf(1,' angvelint %11.7f %11.7f rad/s  \n',angvelint, angvelint /tusec);
     fprintf(1,' angveltgt %11.7f %11.7f rad/s  \n',angveltgt, angveltgt /tusec);
     fprintf(1,' atrans %11.7f %11.7f km \n',atrans, atrans*6378.137  );
     fprintf(1,' ttrans %11.7f %11.7f min \n',ttrans, ttrans*tumin  );
     fprintf(1,' deltatnode  %11.7f %11.7f min \n',deltatnode, deltatnode*tumin );
     fprintf(1,' lead  %11.7f \n',lead*rad );
     fprintf(1,' omeganode  %11.7f \n',omeganode*rad );
     fprintf(1,' phasenew  %11.7f \n',phasenew*rad );
     fprintf(1,' leadnew  %11.7f \n',leadnew*rad );
     fprintf(1,' tphase %11.7f %11.7f min \n',tphase, tphase*tumin  );
     fprintf(1,' aphase  %11.7f %11.7f km \n',aphase, aphase*re );
     fprintf(1,' vint   %11.7f %11.7f km/s \n',vint,vint*velkmps );
     fprintf(1,' vphase  %11.7f %11.7f km/s \n',vphase, vphase*velkmps );
     fprintf(1,' dvphase  %11.7f  %11.7f \n',dvphase, dvphase*velkmps );
     fprintf(1,' vtrans1   %11.7f %11.7f km/s \n',vtrans1, vtrans1*velkmps );
     fprintf(1,' dvtrans1   %11.7f %11.7f km/s \n',dvtrans1, dvtrans1*velkmps );
     fprintf(1,' vtrans2   %11.7f  %11.7f km/s \n',vtrans2, vtrans2*velkmps );
     fprintf(1,' vtgt  %11.7f %11.7f km/s \n',vtgt,vtgt*velkmps );
     fprintf(1,' dvtrans2   %11.7f %11.7f km/s \n',dvtrans2, dvtrans2*velkmps );
     fprintf(1,' dvtotal   %11.7f %11.7f km/s \n',dvtotal, dvtotal*velkmps );
     fprintf(1,' ttotal   %11.7f %11.7f min \n',ttotal, ttotal*tumin );
     