function cap(n1,n2,val)
% Adds the stamp of a capacitor with a value of "val" 
% (Farads) connected between nodes n1 and n2 to the 
% C matrix in circuit representation.
%
%              val
%       n1 o---||---o n2    where C= val (Farads)
%-----------------------------------------------------------
global C   %define global variable

if (n1 ~= 0)
    C(n1,n1) = C(n1,n1) + val;
end

if (n2 ~= 0)
    C(n2,n2) = C(n2,n2) + val;

if (n1 ~= 0) && (n2 ~= 0)
    C(n1,n2) = C(n1,n2) - val;
    C(n2,n1) = C(n2,n1) - val;
end

end %func

