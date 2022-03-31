function ind(n1,n2,val)
% Adds the stamp of an inductor with a value of "val" (Henrys)
% connected between nodes n1 and n2 to the matrices in
% circuit representation.
%
%                 val
%      n1 O-----@@@@@@-----O n2      where L=val (henry)
%              IL ---->
%---------------------------------------------------------------
global G C b   %define global variables

d = size(C,1); % current size of the MNA
xr = d+1;      % new (extera)row/column

% Using an index bigger than the current size,  Matlab
... automatically increases the size of the matrix:
b(xr) = 0;      % add new row
G(xr,xr) = 0;   % add new row/column

C(xr,xr) = -val;  
if (n1 ~= 0)
    G(n1,xr) = 1;
    G(xr,n1) = 1;
end
if (n2 ~= 0)
    G(n2,xr) = -1;
    G(xr,n2) = -1;
end

end %func



