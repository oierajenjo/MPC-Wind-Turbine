function [Ac,Ae,M,W,B,D,To] = extract_struct(var)
Ac = getfield(var,'Ac');
Ae = getfield(var,'Ae');
M = getfield(var,'M');
W = getfield(var,'W');
B = getfield(var,'B');
D = getfield(var,'D');
To = getfield(var,'To');
end
