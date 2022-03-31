function res(n1,n2,val)
% Adds the stamp of a resistor with a value of "val" (Ohms)
% connected between nodes n1 and n2 to the G matrix in 
% circuit representation.
%
%                   val
%       n1 O-----/\/\/\/\---O n2    where R=val (ohms)
%----------------------------------------------------------
global G   %define global variable

if (n1 ~= 0)
    G(n1,n1) = G(n1,n1) + 1/val;
end

if (n2 ~= 0)
    G(n2,n2) = G(n2,n2) + 1/val;
end

if (n1 ~= 0) && (n2 ~= 0)
    G(n1,n2) = G(n1,n2) - 1/val;
    G(n2,n1) = G(n2,n1) - 1/val;
end

end %func
