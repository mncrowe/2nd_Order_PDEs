function R = create_BC(type,N,pos,M,params)
% Creates a boundary condition row
% - type:   1 - Dirichlet (value specified) (default)
%           2 - Neumann (derivative specified)
%           3 - Mixed (a*dp/dx + b*p = c)
% - N: number of gridpoints (default: 101)
% - pos: index of boundary row, i.e. 1 - left, N - right, default: 1
% - M: differentiation matrix (not required for Dirichlet)
% - params: vector [a b], required only for type = 3, default: [1 1]

if nargin < 1; type = 1; end
if nargin < 2; N = 101; end
if nargin < 3; pos = 1; end
if nargin < 5; params = [1 1]; end

I = eye(N);

if type == 1
    R = I(pos,:);
end

if type == 2
    R = M(pos,:);
end

if type == 3
    R = params(1)*M(pos,:)+params(2)*I(pos,:);
end

end

