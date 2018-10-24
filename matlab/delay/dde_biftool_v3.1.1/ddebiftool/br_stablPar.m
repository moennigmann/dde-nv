function branch=br_stablPar(funcs,branch,skip,recompute)
%% compute stability information along branch
% function st_branch=br_stabl(funcs,branch,skip,recompute)
% INPUT:
%   funcs problem function
%	branch 
%	skip number of points to skip between stability computations
%	recompute if nonzero recompute stability info already present
% OUTPUT:
%	st_branch branch with stability information

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%%
ll=length(branch.point);

if ll<1 
  err=ll;
  error('BR_STABL: branch is empty: %d points!',err);
end

if ~isfield(branch.point(1),'stability')
  branch.point(1).stability=[];
end

indexToCalcStab = 1:skip+1:ll-1;
points=branch.point(indexToCalcStab);

stability=branch.method.stability;
parfor i=1:numel(indexToCalcStab)
  if isempty(points(i).stability) || recompute
    points(i).stability=p_stabil(funcs,points(i),stability);
  end
end
branch.point(indexToCalcStab)=points;

if isempty(branch.point(ll).stability) || recompute
  branch.point(ll).stability=p_stabil(funcs,branch.point(ll),branch.method.stability);
end

return;
