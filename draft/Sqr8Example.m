% Examplw on a square domain triangulated by *EIGHT* triangles
% Modified from SqrExample.m
% updated Nov 22, 2014

MofT = 1;
FofT = 0;

% vertices are 0-8
% 
%  4-5-6
%  |\|/|
%  1-2-7
%  |/|\| 
%  0-3-8
%
%  Triangles oriented counter clock-wise; edges oriented lexicographically
%  Edges and triangles are numbered (i.e., ordered) lexicographically

Edges = [ 0 1
          0 2
          0 3
          1 2
          1 4
          2 3
          2 4
          2 5
          2 6
          2 7
          2 8
          3 8
          4 5
          5 6
          6 7
          7 8 ];

m = length(Edges); % # edges = 16

Edges = [(0:m-1)' Edges];

ld = sqrt(2);  % length of diagonal edges
ls = 1;        % length of straight (horizontal/vetical) edges
w = [ls ld ls*ones(1,4) repmat([ld ls],1,3) ls*ones(1,4)]';

Tris = [ 0 2 1
         0 3 2
         1 2 4 
         2 3 8
         2 5 4
	 2 6 5
	 2 7 6
	 2 8 7 ];

n = length(Tris); % # triangles = 8
Tris = [(0:n-1)' Tris];
at = (1/2)*ls^2; % area of each triangle
v = at*ones(n,1);

% Boundary matrix nonzeros (24 nonzeros, 3 per triangle)
nzB = [ 0  0 -1
        1  0  1
        3  0 -1
        1  1 -1
        2  1  1
        5  1 -1
        3  2  1
        4  2 -1
        6  2  1
        5  3  1
       10  3 -1
       11  3  1
        6  4 -1
        7  4  1
       12  4 -1
        7  5 -1
        8  5  1
       13  5 -1
        8  6 -1
        9  6  1
       14  6 -1
        9  7 -1
       10  7  1
       15  7 -1 ];

% Boundary matrix 
B = full(sparse(nzB(:,1)+1,nzB(:,2)+1,nzB(:,3)));

k = 2; % three input currents (or curves)
ti = zeros(m,k);
ti([0 4 12 13]+1,1)  = 1;
ti([2 11 14 15]+1,2) = [1 1 -1 -1];
% ti([1 8]+1,3) = -1;

% Mean Shape Linear Program (MSLP)
% --------------------------------
 
Lambda = 0.01; 

% # rows and # columns in A_MSLP (constraint matrix)
mm = k*m;
nn = 2*m + k*(2*m+2*n);

A_MSLP = repmat([eye(m) -eye(m)],k,1);
for ik=1:k
    A_MSLP = [A_MSLP [zeros((ik-1)*m,2*m+2*n); [-eye(m) eye(m) -B B]; zeros((k-ik)*m,2*m+2*n)] ];
end

% rhs vector
b_MSLP = reshape(ti,k*m,1);

% objective function vector
c_MSLP = (1/k)*[zeros(1,2*m) repmat([w' w' Lambda*[v' v']],1,k)]';

% penalizing mass(T) same as M(Q_i)
cMofT_MSLP = (1/k)*[w' w' repmat([w' w' Lambda*[v' v']],1,k)]';

if MofT==1
   c_MSLP = cMofT_MSLP;
end


% Now adding flat norm of T to the objective function as well
AFofT_MSLP = [              A_MSLP                   zeros(k*m,2*m+2*n)
              eye(m) -eye(m)  zeros(m,k*(2*m+2*n)) [-eye(m) eye(m) -B B] ];

bFofT_MSLP = [b_MSLP; zeros(m,1)];
cFofT_MSLP = [c_MSLP' (1/k)*[w' w' Lambda*[v' v']] ]'; 

if FofT == 1              % solve FofT_MSLP in place of MSLP
   A_MSLP = AFofT_MSLP;
   b_MSLP = bFofT_MSLP;
   c_MSLP = cFofT_MSLP;
end

A_MSLP = round(A_MSLP+0.0000001); % redundant, just to avoid -0 :-)!


[z,MFD,status] = glpk (c_MSLP, A_MSLP, b_MSLP);

t  = z(1:m) - z(m+(1:m));

qi = zeros(m,k);
ri = zeros(n,k);

for ik=1:k
    qi(:,ik) = z(2*m+(ik-1)*(2*m+2*n)+(1:m)) - z(2*m+(ik-1)*(2*m+2*n)+m+(1:m));
    ri(:,ik) = z(2*m+(ik-1)*(2*m+2*n)+2*m+(1:n)) - z(2*m+(ik-1)*(2*m+2*n)+2*m+n+(1:n));
end

disp('t: ');
disp([ find(t)-1 t(find(t))])

if FofT == 1
   x  = z(nn+(1:m)) - z(nn+m+(1:m));
   s  = z(nn+2*m+(1:n)) - z(nn+2*m+n+(1:n));
end

% mySol
