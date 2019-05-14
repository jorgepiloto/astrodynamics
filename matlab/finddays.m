% ------------------------------------------------------------------------------
%
%                           function finddays
%
%  this function finds the fractional days through a year given the year,
%    month, day, hour, minute and second.
%
%  author        : david vallado                  719-573-2600   22 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    year        - year                           1900 .. 2100
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - hour                           0 .. 23
%    min         - minute                         0 .. 59
%    sec         - second                         0.0 .. 59.999
%
%  outputs       :
%    days        - day of year plus fraction of a
%                    day                          days
%
%  locals        :
%    lmonth      - length of months of year
%    i           - index
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2007, 207, ex 3-12
%
% [days] = finddays ( year,month,day,hr,min,sec);
% -----------------------------------------------------------------------------

function [days] = finddays ( year,month,day,hr,min,sec);

        for i= 1 : 12
            lmonth(i) = 31;
            if i == 2
                lmonth(i)= 28;
              end;
            if i == 4 | i == 6 | i == 9 | i == 11
                lmonth(i)= 30;
              end;
        end

        if (rem (year,4) == 0)
	    lmonth(2)= 29;
	    if (rem (year,100) == 0) & (rem (year,400) ~= 0)
		lmonth(2)= 28;
            end
	end

        i   = 1;
        days= 0.0;
        while (i < month) & ( i < 12 )
            days= days + lmonth(i);
            i= i + 1;
          end

        days= days + day + hr/24.0 + min/1440.0 + sec/86400.0;

