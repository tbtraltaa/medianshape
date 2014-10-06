function v=simpvol2(p,t)
%SIMPVOL2 Simplex volume in Euclidean space.
%   V=SIMPVOL2(P,T)
% p = matrix of point coordinates.  Each column gives the
%   coordinates of a single point in Euclidean space.
% t = matrix of simplices.  Each column corresponds to a
%   single d-dimensional simplex and has d+1 indices corresponding
%   to columns of p.

% V = vector of simplex volumes, one entry per column of t

% Notes:
% simpvol from DistMesh can not calculate volumes of d-dimensional
% simplices where d is less than the dimension of the ambient space.
% This function does so by projecting into the d-dimensional
% subspace containing the simplex and doing the calculation there.
% This will fail if the d+1 points are contained in a subspace of
% dimension d-1 or less.
% Note that p and t are transposed from their use in simpvol.
% This is so that this function can be more easily called from 
% msfn and so that things are processed in a column major order.

% d = simplex dimension
d = size(t,1) - 1;

% Initialize simplex volume array
v=zeros(size(t,2),1);

for ii=1:size(t,2)
    pts = p(:,t(:,ii));
    vec = pts(:,2:end) - pts(:,1)*ones(1,d);
    if (size(vec,1) == size(vec,2))
        v(ii) = det(vec)/factorial(d);
        continue
    end
    
    % Several possible implementations.  Still working on this.
    
    %[Q,R] = qr(vec,0);
    %v(ii) = prod(diag(R))/factorial(d); % actually gives abs(vol)

    %[Q,R] = qr(vec);
    %v0 = det(Q)*prod(R(1:size(R,2)+1:size(R,2)))/factorial(d);
    % prod(diag(R)) is goal since we want the ``diagonal'' of the
    % rectangular matrix R. If R is 2 by 1 (like a vector), diag
    % gives the wrong result (i.e., a matrix instead of a single 
    % entry).  Thus the prod(R(1:size(R,2)+1:size(R,2)).
    % This attempts to provide a signed volume by the det(Q) factor
    % but this isn't rigorous.
    
    % Calculate an orthonormal basis for the simplex subspace
    % so we can use the fact that it is an isometry on points in
    % the subspace and have det(ob'*vec)/factorial(d) give the 
    % correct volume.
    ob = orth(vec);
    v(ii) = det(vec'*ob)/factorial(d);
end