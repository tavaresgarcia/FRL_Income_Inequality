function A = varcompanion(Aols)

[nA kA] = size(Aols);
if nA > kA; Aols = Aols'; end;
[nA kA] = size(Aols);
porder = kA/nA;

A = [Aols; [kron(eye(porder-1),eye(nA)) zeros((porder-1)*nA,nA)]];