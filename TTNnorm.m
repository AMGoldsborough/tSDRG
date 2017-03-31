function norm = TTNnorm(L,w,Jorder)
% norm = TTNnorm(w)
% 
% input: length (L), set of isomteries (w), order in which the TTN is 
% constructed (Jorder)
% output: norm of wavefunction

% Andrew Goldsborough 27/08/2013
% function that calculates the norm of a TTN wavefunction

%create a set of identities as the function
MPOeye = cell(L,1);

for i=1:L
    MPOeye{i} = eye(2);
end

for i=1:L-1
    
    %contract iso with operator
    MPOeye{Jorder(i)} = tcon(w{i,1},MPOeye{Jorder(i)},[-1,1,-3],[1,-2]);
    MPOeye{Jorder(i)} = tcon(MPOeye{Jorder(i)},MPOeye{Jorder(i)+1},[-1,-2,1],[1,-3]);
    MPOeye{Jorder(i)} = tcon(MPOeye{Jorder(i)},conj(w{i,1}),[-1,1,2],[-2,1,2]);
    
    %remove extra MPO tensor
    MPOeye(Jorder(i)+1) = [];
end

% the norm should be the final MPO
norm = MPOeye{1};