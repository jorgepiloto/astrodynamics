% ----------------------------------------------------------------------------
%
%                           function coveq2cl
%
%  this function transforms a six by six covariance matrix expressed in
%    equinoctial elements into one expressed in classical orbital elements.
%
%  author        : david vallado                  719-573-2600   18 jun 2002
%
%  revisions
%    vallado     - major update                                  26 aug 2015
%
%  inputs          description                    range / units
%    eqcov       - 6x6 equinoctial covariance matrix
%    eqstate     - 6x1 equinoctial orbit state    (a/n af ag chi psi lm/ln)
%    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
%    fr          - retrograde factor               +1, -1
%
%  outputs       :
%    classcov    - 6x6 classical covariance matrix
%    tm          - transformation matrix
%
%  locals        :
%    n           - mean motion                    rad
%    af          - component of ecc vector
%    ag          - component of ecc vector
%    chi         - component of node vector in eqw
%    psi         - component of node vector in eqw
%    meanlon     - mean longitude                 rad
%    mum         - gravitational paramater        m^3/s^2 NOTE Meters!
%
%  coupling      :
%    constastro
%
%  references    :
%    Vallado and Alfano 2015
%
% [classcov, tm] = coveq2cl eqcov, eqstate, anom, fr);
% ----------------------------------------------------------------------------
 
function [classcov, tm] = coveq2cl(eqcov, eqstate, anom, fr)
 
        % -------- define gravitational constant
        constastro;

        % -------- parse the orbit state
        % --------- determine which set of variables is in use ---------
        if strcmp(anom,'truea') == 1 || strcmp(anom,'meana') == 1  % 1 is true
            a = eqstate(1);  % in m
            n = sqrt(mum/a^3);
        else
            if strcmp(anom,'truen') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                 n      = eqstate(1);
                 a = (mum / n^2)^(1/3);
            end
        end
        af      = eqstate(2);
        ag      = eqstate(3);
        chi     = eqstate(4);
        psi     = eqstate(5);
        if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
            meanlonM = eqstate(6);
          else
            if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
                meanlonNu = eqstate(6);
            end
        end

        
        % ---------------- calculate matrix elements ------------------
        % ---- partials of a wrt (a af ag chi psi l)
        if strcmp(anom,'truea') == 1 || strcmp(anom,'meana') == 1  % 1 is true
            tm(1,1) =  1.0;
        else
            if strcmp(anom,'truen') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                tm(1,1) = -2.0/(3*n) * (mum/n^2)^(1/3);  % if class = a, equin = n
            end
        end    
        tm(1,2) = 0.0;
        tm(1,3) = 0.0;
        tm(1,4) = 0.0;
        tm(1,5) = 0.0;
        tm(1,6) = 0.0;

        % ---- partials of ecc wrt (n af ag chi psi l)
        p0 = 1.0 / sqrt(af^2 + ag^2);
        tm(2,1) = 0.0;
        tm(2,2) = p0*af;
        tm(2,3) = p0*ag;
        tm(2,4) = 0.0;
        tm(2,5) = 0.0;
        tm(2,6) = 0.0;

        % ---- partials of incl wrt (n af ag chi psi l)
        p1 = 2.0*fr / ((1.0 + chi^2 + psi^2) * sqrt(chi^2 + psi^2));
        tm(3,1) = 0.0;
        tm(3,2) = 0.0;
        tm(3,3) = 0.0;
        tm(3,4) = p1*chi;
        tm(3,5) = p1*psi;
        tm(3,6) = 0;

        % ---- partials of omega wrt (n af ag chi psi l)
        p2 = 1.0 / (chi^2 + psi^2);
        tm(4,1) =  0;
        tm(4,2) =  0;
        tm(4,3) =  0;
        tm(4,4) =  p2*psi;
        tm(4,5) = -p2*chi;
        tm(4,6) =  0;

        % ---- partials of argp wrt (n af ag chi psi l)
        p3 = 1.0 / (af^2 + ag^2);
        tm(5,1) = 0.0;
        tm(5,2) = -p3*ag;
        tm(5,3) =  p3*af;
        tm(5,4) = -fr*p2*psi;
        tm(5,5) =  fr*p2*chi;
        tm(5,6) = 0.0;

        % ---- partials of anomaly wrt (n af ag chi psi l)
        p4 = 1.0 / (af^2 + ag^2);
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            tm(6,1) = 0.0;
            tm(6,2) = p4*ag;  % p3 on these. same
            tm(6,3) = -p4*af;
            tm(6,4) = 0.0;
            tm(6,5) = 0.0;
            tm(6,6) = 1.0;
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
               tm(6,1) = 0.0;
               tm(6,2) = p3*ag; 
               tm(6,3) = -p3*af; 
               tm(6,4) = 0.0;   
               tm(6,5) = 0.0;    
               tm(6,6) = 1.0; 
            end
        end

        % ---------- calculate the output covariance matrix -----------
        [classcov] =  tm*eqcov*tm';


