function [zx,zy] = trigradient(x,y,z,t,opt)
%   
if nargin<4; t = delaunay(x,y); opt ='normal';end
p = [x,y];

if nargin<5; opt='normal'; end

nt = size(t,1); np = size(p,1);
% [x1 y1 1; x2 y2 1; x3 y3 1][a; b; c] = [z1; z2; z3] 
C = sparse(repmat(1:3*nt,3,1)',kron(reshape(1:3*nt,3,nt)',[1; 1; 1]),[p(t',:),ones(3*nt,1)])\z(t',:);

if strcmp(opt,'face')
    zx = C(1:3:end,:);
    zy = C(2:3:end,:);
    return
end

areas = repmat(.5*abs((p(t(:,2),1)-p(t(:,1),1)).*(p(t(:,3),2)-p(t(:,1),2))-(p(t(:,3),1)-p(t(:,1),1)).*(p(t(:,2),2)-p(t(:,1),2))),1,size(z,2));

M = sparse(repmat((1:nt)',3,1),t,1,nt,np);

zx = (M'*(C(1:3:end,:).*areas))./(M'*areas);
zy = (M'*(C(2:3:end,:).*areas))./(M'*areas);