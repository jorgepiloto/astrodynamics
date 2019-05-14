% -----------------------------------------------------------------------------
%
%                           function inittime.m
%
%  this function initializes some of the time variables.
%
%  revisions
%                -
%
%  inputs          description                    range / units
%
%
%  outputs       :
%
%
%  locals        :
%    i           - index
%
%  coupling      :
%    none
%
% -----------------------------------------------------------------------------

        % ------------------------  implementation   ------------------

        % inittime.m

        for i= 1 : 12
            switch i
             case{1,3,5,7,8,10,12}
                 lmonth(i)= 31;
             case{4,6,9,11}
                 lmonth(i)= 30;
             case 2
                 lmonth(i)= 28;
            end
        end

        monthtitle( 1)= 'jan';
        monthtitle( 2)= 'feb';
        monthtitle( 3)= 'mar';
        monthtitle( 4)= 'apr';
        monthtitle( 5)= 'may';
        monthtitle( 6)= 'jun';
        monthtitle( 7)= 'jul';
        monthtitle( 8)= 'aug';
        monthtitle( 9)= 'sep';
        monthtitle(10)= 'oct';
        monthtitle(11)= 'nov';
        monthtitle(12)= 'dec';

        daytitle(1)= 'sun';
        daytitle(2)= 'mon';
        daytitle(3)= 'tue';
        daytitle(4)= 'wed';
        daytitle(5)= 'thr';
        daytitle(6)= 'fri';
        daytitle(7)= 'sat';

