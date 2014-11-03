% Updated on July 5, 2012 by Bala
% Now using glpk to solve LP in place of linprog

function [ norm, x, s, v, w, cons] = msfn( p, K, M, t, lambda, v, w, cons)
%MSFN Calculates multiscale flatnorm on simplicial complexes
%    [ norm, x, s, v, w, bd] = msfn( p, K, M, t, lambda, v, w, bd)
%
% p is a matrix with a column for each corner point in simplex
%   Generally, it will be the coordinates of each point in the
%   ambient space, but this isn't necessary.
%   p is only used to calculate volumes and orientations if not
%   provided
% K is a (d+2) by n matrix encoding all d+1-dimensional simplices 
%   in the complex.  Each column specifies a single simplex by
%   listing the d+2 point indices (where the indices correspond
%   to columns of p)
% M is a (d+1) by m matrix encoding all d-dimensional simplices
%   in the complex.  This may seem redundant since this could be
%   generated from K (see simpbd), but the relevant detail here
%   is that it establishes an ordering of the d-dimensional
%   simplices (and thus how to interpret t).
% t is the d-dimensional simplicial current for which we want 
%   to calculate the flat norm.  Pass [] to do everything but
%   calculate flat norm (for precomputing complex parameters).
% lambda is the flat norm scale constant
% v,w are the volumes of the d+1- and d-dimensional simplices
%   v is a vector of length n, w is of length m
%   The entries in v and w may be signed to indicate orientation.
%   If not provided, these will be calculated using Euclidean
%   volumes.
% cons is the constraint matrix [I -I bd -bd] described in the 
%   paper.  There is generally no need other than performance
%   to specify this directly as it can be calculated from v & w.

% Output
% norm is the flat norm
% x,s are the flat norm decomposition into the d- and 
%   d+1-dimensional components.
% v,w,bd are the same as in the input. These outputs are provided
%   so they can be reused if calling the function multiple times
%   on the same complex.

% n = number of simplices, s(1:n)
n = size(K, 2);
% m = number of codimension 1 simplices
m = size(M, 2);

d = size(K, 1) - 2;

%TODO: Check for argument consistency; make sure that sizes are compatible.

if ~exist('v', 'var')
    v = simpvol2(p, K);
end
if ~exist('w', 'var')
    w = simpvol2(p, M);
end

% Create constraint matrix if necessary
if ~exist('cons', 'var')
    % Goal: build cons = [eye(m) -eye(m) bd -bd]   
    % Initialize with eye(m) -eye(m)
    cons = speye(m,2*n+2*m);
    cons(:,m+1:2*m) = -speye(m);
    % Put vertices in ascending index order for each simplex
    % so they can be compared easily.  However, this could
    % change the orientation so we store the parity of the
    % permutation applied.
    [K, kcolindex] = sort(K, 1);
    kparity = permutationparity(kcolindex,1);
    [M, mcolindex] = sort(M, 1);
    mparity = permutationparity(mcolindex,1);
    for ii = 1:n % For each d+1 dimensional simplex
        % We consider each d dimensional subsimplex by
        % taking all possible sets of d+1 vertices from
        % the d+2 available by excluding them one at a time.
        for jj = 1:d+2 
            subsimplex = K(:,ii);
            subsimplex(jj) = [];
            [~, ssidx] = ismember(subsimplex', M', 'rows');
            if (ssidx == 0)
                error('Unable to find subsimplex! Make sure M contains all boundary simplices.');
            end
            val = (-1)^(jj + 1 + kparity(ii) + mparity(ssidx));
            cons(ssidx,2*m+ii) = val;
            cons(ssidx,2*m+n+ii) = -val;
        end
    end
end

% If no input was specified, return without calculating flat norm.
% This is useful for precomputing the complex-dependent parameters.
if ~exist('t','var') || isempty(t)
    norm = 0;
    x = [];
    s = [];
    return
end

% Variable order: [xplus, xminus, splus, sminus]
% Objection function = w^T(xplus + xminus) + lambda v^T(splus+sminus)
% w and v should be positive for this so we take absolute values.
c = [abs(w);abs(w);lambda*abs(v);lambda*abs(v)];
% cons (defined above) encodes the xplus - xminus + bd(splus - sminus) = t constraints
% Force Simplex solver so we get a corner optimal point
% SWITCHING OFF FOR GLPK!!
% options = optimset('LargeScale', 'off', 'Simplex', 'on');
% Uncomment to restore default (faster) interior point method
%options = optimset(@linprog);
% TODO: Implement algorithm to go from interior point solution to a corner
% point.

% Changing from linprog to glpk now...
% [arg, norm, exitflag] = linprog(c, [], [], cons, t, zeros(2*(n+m),1), [], [], options);

% ctype = repmat("S",1,m);
% vartype = repmat("C",1,2*(m+n));

[arg, norm, exitflag, extras] = glpk(c, cons, t);

% if exitflag ~= 1 
if exitflag ~= 180 % in glpk now
    error('Optimization failed.');
end

x = arg(1:m)-arg(m+1:2*m); % xplus - xminus
s = arg(2*m+1:2*m+n)-arg(2*m+n+1:end); % splus - sminus

end

