function [p,K,M,t,v,w,cons] = msfndemo2(p,K,M,t,v,w,cons)
%MSFNDEMO Automatic MSFN mesh creation and solving.

rand('state',112); % Always the same results
set(gcf,'rend','opengl');

if ~exist('M', 'var')
    tic
end

if ~exist('K', 'var')
    %disp('Adapted from DistMesh example (1a) Unit circle')
    fd=inline('sqrt(sum(p.^2,2))-1','p');
    [p,K]=distmesh2d(fd,@huniform,0.05,[-1,-1;1,1],[]);
    
    % Uncomment for another 2D mesh (DistMesh tends to get stuck on this one)
    %disp('Adapted from DistMesh example (5) Geometric Adaptivity')
    %fix=[-1,0;-.95,0;.1,0;1,0];
    %[p,K]=distmesh2d(@fd5,@fh5,0.025,[-1,0;1,1],fix);
    
    %disp('DistMesh example (10) Cylinder with hole')
    %[p,K]=distmeshnd(@fd10,@fh10,0.1,[-1,-1,-1;1,1,1],[]);
    
    p = p';
    K = K';
end

n = size(p,1);
d = size(K,1)-2;

if ~exist('M', 'var')
    v = simpvol2(p,K);
    % Use signed volumes to fix orientations (works when dimension
    % of complex is same as the ambient space, not harmful in other cases)
    if (size(K,1) >= 2)
        kparity = v < 0;
        for ii = 1:size(K,2)
            if kparity(ii)
                K([1 2],ii) = K([2 1],ii);
                v(ii) = -v(ii);
            end
        end
    end
    M = simpbd(K);
    w = abs(simpvol2(p,M));
    elapsed = toc;
    fprintf('Mesh calculation took %f seconds.\n', elapsed);
    
end

if ~exist('t', 'var')
    tic;
    
    % Tweak this to sample more/less points along our curve
    numpts = 500;
    domain = linspace(0,1,numpts);
    dhalf = linspace(0,1,numpts/2);
    % Piecewise linear curve from left to top to right of circle
    y = [[0; .99]*dhalf+[-.99; 0]*(1-dhalf)  [.99; 0]*dhalf+[0;.99]*(1-dhalf)];
    % Uncomment to get straight line from left to right of circle
    %y = [.99; 0]*domain + [-.99; 0]*(1-domain);
    % Uncomment to get the upper half of the circle in R2
    %x = 2*(domain-.5);
    %y = [1.01*x; 1.01*sqrt(1-x.*x) + .0001];
    % Uncomment for a curve in R3:
    %y = [[0;0]*domain + [1;1]*(1-domain) ;sin(2*pi*domain)];
    
    t = zeros(size(M,2), 1);
    if d == 1
        % Method that follows is for finding 1D boundaries on a 2D mesh (in any
        % dimension).  See below for other versions.
        P =  1:size(p,2);
        [~, closestpoint] = min(simpdist(p, P, y));
        
        [Msorted, mcolindex] = sort(M, 1);
        mparity = permutationparity(mcolindex,1);
        Msorted = Msorted';
        
        lastpt = closestpoint(1);
        for ii = 2:length(closestpoint)
            if closestpoint(ii) == lastpt
                continue
            end
            linesegment = [lastpt closestpoint(ii)];
            line(p(1,linesegment)', p(2,linesegment)','Color','red');
            orientation = 1;
            % Sort segment indices to look it up, noting original orientation.
            if closestpoint(ii) < lastpt
                linesegment = linesegment([2 1]);
                orientation = -1;
            end
            [~, loc] = ismember(linesegment, Msorted, 'rows');
            if loc == 0
                error('Line not found.  Try sampling more points.')
                continue
            end
            if mparity(loc) == 1
                orientation = -orientation;
            end
            t(loc) = orientation;
            lastpt = closestpoint(ii);
        end
    else
        % Variant on previous method for higher dimensions.  Very slow (not
        % particularly optimized and does a quadratic programming problem
        % for every pair of sampled points on curve and simplices in M.
        % Notwithstanding the performance, the results are generally not
        % nice either.
        [~, closestsimplex] = min(simpdist(p, M, y));
        simplices = unique(closestsimplex);
        %%simplices = closestsimplex(simpidx);
        t(simplices) = ones(length(simplices), 1);
    end
    
    % Old method for generating t randomly:
    % Create t with entries in -1, 0, 1 with equal probability
    %randvec = randn(1,size(M, 2));
    %t = double(randvec >= 2/3) - double(randvec <= 1/3);
    %t = double(randn(1,size(M, 2)) >= 0.5);
    
    elapsed = toc;
    fprintf('Coercing points to complex took %f seconds.\n', elapsed);
    
end

if ~exist('v', 'var')
    v = abs(simpvol2(p, K));
end
if ~exist('w', 'var')
    w = abs(simpvol2(p, M));
end


if ~exist('cons', 'var')
    tic
    % Precompute simplicial complex invariants.
    % Since v,w and cons are independent of t and lambda, it's a
    % good idea to save them whenever we want to solve more than
    % one instance of MSFN on the same grid.
    [~, ~, ~, ~, ~, cons] = msfn(p, K, M, [], 0, v, w);
    elapsed = toc;
    fprintf('MSFN (init) took %f seconds.\n', elapsed);
end

tic
inputsimp = M(:,find(t));
% Calculate flat norm using saved invariants for various values of lambda
for lambda = [0 0.1 0.5 1 10 50]
    [norm, x, s] = msfn(p, K, M, t, lambda, v, w, cons);
    % TODO: Quick fix for when using interior point method; may not work in
    % most cases but seems to return corner points for most instances.
    x = round(x);
    if n == 2
        figure('Name',sprintf('Lambda = %f',lambda));
        %axis([-1.5 1.5 -.5 1.5]);
        hold on;
        for ii = 1:size(inputsimp,2)
            line(p(1,inputsimp(:,ii))', p(2,inputsimp(:,ii))','Color','blue');
        end
        outputsimp = M(:,find(x));
        for ii = 1:size(outputsimp,2)
            line(p(1,outputsimp(:,ii))', p(2,outputsimp(:,ii))','Color','red');
        end
        title(sprintf('Lambda = %f',lambda));
        text(.5, 1,sprintf('Norm = %f',norm));
        hold off;
    end
end

elapsed = toc;
fprintf('MSFN took %f seconds.\n', elapsed);


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


function d=fd5(p)

d1=dcircle(p,0,0,1);
d2=dcircle(p,-0.01,0.01,.7);
d=dintersect(-p(:,2),ddiff(d1,d2));

function h=fh5(p)

d1=dcircle(p,0,0,1);
d2=dcircle(p,0,0,.55);

h1=(0.15-0.2*d1);
h2=(0.06+0.2*d2);
h3=(d2-d1)/3;

h=min(min(h1,h2),h3);
