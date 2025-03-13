function F = create_RHS(f,varargin)
% Creates a reshaped RHS vector with boundary condition values included
% - f: RHS field, n dimensional array
% - BC_val1, BC_val2: pairs of RHS boundary condition values;
%       [a b] - bottom and top values of a and b, scalars
%       array - n dimensional array, with last dimension of size 2 corresponding to bottom and top BCs
%       s - scalar, do not apply BC, useful for perodic directions

s = size(f);
n = length(s(s~=1));
app_BC = ones(1,n);

for in = 1:n
    
    temp = varargin{in};

    if numel(temp) == 2
        t1 = s; t1(in) = 2;
        t2 = zeros([ones(1,n-1) 2]); t2(1) = temp(1); t2(2) = temp(2);
        temp = ones(t1).*t2;
    end
    
    if numel(temp) == 1
        app_BC(in) = 0;
    end

    BC_val{in} = temp;
    
end

if n == 1
    
    if app_BC == 1
        f(1) = BC_val{1}(1);
        f(end) = BC_val{1}(2);
    end
    
else
    
    s1 = s(s~=1);
    
    for in = 1:n
        
        if app_BC(in) == 1
            
            f = reshape(f,[prod(s1(1:in-1)) s1(in) prod(s1(in+1:end))]);
            BC = reshape(BC_val{in},[prod(s1(1:in-1)) 1 prod(s1(in+1:end)) 2]);

            f(:,1,:) = BC(:,:,:,1);
            f(:,end,:) = BC(:,:,:,2);
        
        end
        
    end
    
end

F = reshape(f,[prod(s1) 1]);

end