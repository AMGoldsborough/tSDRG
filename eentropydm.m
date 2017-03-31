function ee = eentropydm(dm)
% Entanglement (Von Neumann) Entropy calculated from the density matrix
% input: density matrix dm
% output: the entanglement entropy

% Andrew Goldsborough 19/02/2013
% function to calculate the entanglement (Von Neumann) entropy. The
% Equation is ee = - tr(dm * log2(dm));

%make sure that the dm's are hermitean
dm = (dm + dm')/2;

[dmv,dmnrg] = eig(dm);

nzloc = find(diag(dmnrg)>=1E-15);

ee = -trace(dm*dmv*((diag(log2(diag(dmnrg)))+diag(log2(diag(dmnrg)))')/2)*dmv');

