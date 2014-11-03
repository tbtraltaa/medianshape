function M = simpbd( K )
%simpbd Given a d+1-dimensional simplicial complex K, 
%   find all d-dimensional subsimplices.
% K is an (d+2) by n matrix where n is the number of simplices
%   in the complex and the top-level simplices are d+1 dimensional

% Put indices in ascending order for each simplex so
% they can be compared easily.
[K, kcolindex] = sort(K, 1);
kparity = permutationparity(kcolindex,1);

% n = number of simplices, s(1:n)
n = size(K, 2);
% Top-level simplices are d+1-dimensional
d = size(K, 1) - 2;

nextm = 1;
% Allocate M with its minimum possible size (and grow later, see TODO)
M = zeros((d+2)*(d+3)/2, d+1);
mparity = ones((d+2)*(d+3)/2, 1);

for ii = 1:n % For each d+1 dimensional simplex
    simplex = K(:,ii);
    % We consider each d dimensional subsimplex by taking all
    % possible sets of d+1 vertices from the d+2 available
    % (i.e., exclude them one at a time).
    for jj = 1:d+2
        subsimplex = simplex;
        subsimplex(jj) = [];
        M(nextm,:) = subsimplex;
        mparity(nextm) = xor(kparity(ii), mod(jj, 2) == 1);
        
        % TODO: Grow M and mparity more intelligently.  
        % Either allocate M at maximum size n*(d+2)
        % or double its size whenever we need to expand it 
        % (up to the maximum n*(d+2)).
        % In either case, trim to size using nextm after loop.
        % Low priority since this function is not a computational
        % hot spot.
        nextm = nextm + 1;
    end
end
[M, s] = unique(M, 'rows');

% Fix up orientations on returned boundary simplices
if d >= 1 % Otherwise, boundary simplices are 0-dimensional and our fix won't work
    for ii = 1:size(M,1)
        if mparity(s(ii))
            M(ii,[1 2]) = M(ii, [2 1]);
        end
    end
end

M = M';