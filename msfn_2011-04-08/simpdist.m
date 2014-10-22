function [dist, proj]=simpdist(p,t,x)
%SIMPDIST Find distance from points to simplices in Euclidean space.
%   V=SIMPDIST(P,T,X)
% p = matrix of point coordinates.  Each column gives the
%   coordinates of a single point in R^n.
% t = matrix of simplices.  Each column corresponds to a
%   single d-dimensional simplex and has d+1 indices corresponding
%   to columns of p.
% x = matrix of points for which to find simplex distances.  Each
%   column contains coordinates for a single point in R^n for which
%   to calculate Euclidean distance to each simplex from t.
% dist = distances to each simplex (one row per simplex,
%   one column per point in x)
% proj = 3D array with closest point for each point/simplex pair. 
%   proj(:,ii,jj) = coordinates for closest point on simplex ii
%                   to point x(:,jj).

m = size(t, 2);
numpts = size(x, 2);
n = size(x, 1);
d = size(t, 1)-1;
dist = zeros(m, numpts);
proj = zeros(n, m, numpts);

% Handle point case separately
if (d == 0)
    proj = p(:,t);
    dist = pdist2(proj', x');
    return
end

options = optimset('Display', 'off', 'LargeScale', 'off');

for ii = 1:m
    pts = p(:,t(:,ii));
    %vec = pts(:,2:end) - pts(:,1)*ones(1,d);
    H = [eye(n) zeros(n, d+1); zeros(d+1, n+d+1)];
    Aeq = [eye(n) pts; zeros(1,n) ones(1, d+1)];
    for jj = 1:numpts
        beq = [x(:,jj); 1];
        lb = [-inf*ones(1, n) zeros(1, d+1)];
        argmin = quadprog(H, zeros(n+d+1,1), [], [], Aeq, beq, lb, [], [], options);
        % Uncomment to minimize l1 distance
        %argmin = linprog([ones(1,n) zeros(1, d+1)], [], [], Aeq, beq, lb, [], [], options);
        dist(ii,jj) = sqrt(argmin(1:n)'*argmin(1:n));
        proj(:,ii,jj) = pts*argmin(n+1:end);
    end
end