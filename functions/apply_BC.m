function D = apply_BC(D,N,varargin)
% Applies BCs to an operator D
% - D: operator matrix
% - N: dimensions, [Nx1 Nx2 Nx3 ...]
% - BC: enter length(N) comma separated entries;
%       2 x Nxi matrix; here row 1 is the top BC and row 2 the bottom BC
%       2 x 1 vector; value of 'a' for Dirichlet condition a*f = b, enter [0 0] to set BC rows to zero
%       m; any scalar, do not apply BC in this dimension, e.g. periodic condition

n = size(N); if n(1)>1; N=N'; end; n = length(N);
r = prod(N);

if ~isequal(size(D),[r r]); error('Dimensions of D do not match product of [Nx1 Nx2 ...]'); end

for in = 1:n
    
    M_bc{in} = varargin{in};

    if isequal(size(M_bc{in}),[2 1]); M_bc{in} = [M_bc{in}(1) zeros(1,N(in)-1); zeros(1,N(in)-1) M_bc{in}(2)]; end
    if numel(M_bc{in}) == 1; M_bc{in} = 0; end
    
end

for in = 1:n                        % apply BCs
    
    if ~isequal(M_bc{in},0)
    
        b = prod([1 N(1:in-1)]);   % block width
        d = prod([1 N(1:in)]);     % distance between blocks

        rows1 = [];                 % top BC row
        rows2 = [];                 % bottom BC row

        for i2 = 1:(r/prod(N(1:in)))
            rows1 = [rows1 (1:b)+(i2-1)*d];
            rows2 = [rows2 (d-b+1:d)+(i2-1)*d];
        end
        
        ii = [];
        jj = [];
        kk = [];
        
        for r1 = rows1              % set top BC
            ii = [ii; r1*ones(N(in),1)];
            jj = [jj; (r1:b:r1+(N(in)-1)*b)'];
            kk = [kk; (M_bc{in}(1,:))'];
        end
        for r2 = rows2              % set bottom BC
            ii = [ii; r2*ones(N(in),1)];
            jj = [jj; (r2-(N(in)-1)*b:b:r2)'];
            kk = [kk; (M_bc{in}(2,:))'];
        end
        
        D2 = sparse(ii,jj,kk,r,r,max(nnz(D),length(ii)));   % create sparse matrix of boundary condition rows
        
        i1 = setxor(setxor(1:r,rows1),rows2);
        D2(i1,:) = D(i1,:);                 % set non BC rows to the rows of D
        
        D = D2;
    
    end
    
end

end