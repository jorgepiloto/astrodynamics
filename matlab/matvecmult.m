% ------------------------------------------------------------------------------
%
%                           function matvecmult
%
%  this function multiplies a matrix by a vector.
%
%  author        : david vallado                  719-573-2600    4 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    mat         - matrix (square)
%    vec         - vector
%    size        - dimension of matrix
%
%  outputs       :
%    outvec      - unit vector
%
%  locals        :
%    i,j         - index
%
%  coupling      :
%    none
%
% [outvec] = matvecmult ( mat, vec, sizeof );
% ------------------------------------------------------------------------------

function [outvec] = matvecmult ( mat, vec, sizeof );
    % -------------------------  implementation   -----------------
    for i = 1 : sizeof
        outvec(i) = 0.0;
        for j = 1 : sizeof
            outvec(i) = outvec(i) + mat(i,j) * vec(j);
        end    
    end
            

