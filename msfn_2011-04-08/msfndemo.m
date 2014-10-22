function msfndemo
%MSFNDEMO Automatic MSFN mesh creation and solving.

rand('state',111); % Always the same results
set(gcf,'rend','opengl');

disp('(10) Cylinder with hole')
tic
[p,K]=distmeshnd(@fd10,@fh10,0.1,[-1,-1,-1;1,1,1],[]);
elapsed = toc;
fprintf('Mesh calculation took %f seconds.\n', elapsed);
tic
p = p';
K = K';
if (size(K,1) >= 2)
    kparity = simpvol2(p,K) < 0;
    for ii = 1:size(K,2)
        if kparity(ii)
            K([1 2],ii) = K([2 1],ii);
        end
    end
end
M = simpbd(K);

% Create t with entries in -1, 0, 1 with equal probability
randvec = randn(1,size(M, 2));
t = double(randvec >= 2/3) - double(randvec <= 1/3);
%t = double(randn(1,size(M, 2)) >= 0.5);
lambda = 1;
elapsed = toc;
fprintf('Boundary calculation took %f seconds.\n', elapsed);
tic
% Calculate flat norm, saving simplicial complex invariants.
% Since v,w and bd are independent of t and lambda, it's a
% good idea to save them whenever we want to solve more than
% one instance of MSFN on the same grid.
[norm, x, s, v, w, cons] = msfn(p, K, M, t, lambda);
elapsed = toc;
fprintf('MSFN (uncached) took %f seconds.\n', elapsed);
%cons = full(cons);
tic
% Calculate flat norm using saved invariants
[norm, x, s] = msfn(p, K, M, t, lambda, v, w, cons);
elapsed = toc;
fprintf('MSFN (cached) took %f seconds.\n', elapsed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d=fd10(p)

r=sqrt(p(:,1).^2+p(:,2).^2);
z=p(:,3);

d1=r-1;
d2=z-1;
d3=-z-1;
d4=sqrt(d1.^2+d2.^2);
d5=sqrt(d1.^2+d3.^2);
d=dintersect(dintersect(d1,d2),d3);
ix=d1>0 & d2>0;
d(ix)=d4(ix);
ix=d1>0 & d3>0;
d(ix)=d5(ix);

d=ddiff(d,dsphere(p,0,0,0,0.5));

function h=fh10(p)

h1=4*sqrt(sum(p.^2,2))-1;
h=min(h1,2);

