function [W,d] = heishamhalfSD(L,J,Jz,h)
% Spin 1/2 Heisenberg model hamiltonian
% input: chain length L, array of interaction strengths J, relative z strength Jz, array of on site strengths h
% output: MPOs for the Heisenberg Hamiltonian W, spin dimension d

% Andrew Goldsborough 09/06/2012
% Function that generates the W tensors for a spin-1/2 Heisenberg Hamiltonian and defines the spin dimension

W = cell(L,1);

%spin-1/2 Heisenberg 
W{1} = zeros(2,2,5,5);
W{1}(:,:,1,1) = [1 0;0 1];
W{1}(:,:,1,2) = [0 0.5*J(1);0 0];
W{1}(:,:,1,3) = [0 0;0.5*J(1) 0];
W{1}(:,:,1,4) = [0.5*J(1)*Jz 0;0 -0.5*J(1)*Jz];
W{1}(:,:,1,5) = [0.5*h(1) 0;0 -0.5*h(1)];

for i=2:(L-1)
    W{i} = zeros(2,2,5,5);
    W{i}(:,:,1,1) = [1 0;0 1];
    W{i}(:,:,1,2) = [0 0.5*J(i);0 0];
    W{i}(:,:,1,3) = [0 0;0.5*J(i) 0];
    W{i}(:,:,1,4) = [0.5*J(i)*Jz 0;0 -0.5*J(i)*Jz];
    W{i}(:,:,1,5) = [0.5*h(i) 0;0 -0.5*h(i)];
    W{i}(:,:,2,5) = [0 0;1 0];
    W{i}(:,:,3,5) = [0 1;0 0];
    W{i}(:,:,4,5) = [0.5 0;0 -0.5];
    W{i}(:,:,5,5) = [1 0;0 1];
end

W{L} = zeros(2,2,5,5);
W{L}(:,:,1,5) = [0.5*h(L) 0;0 -0.5*h(L)];
W{L}(:,:,2,5) = [0 0;1 0];
W{L}(:,:,3,5) = [0 1;0 0];
W{L}(:,:,4,5) = [0.5 0;0 -0.5];
W{L}(:,:,5,5) = [1 0;0 1];

%spin dimension for spin-1/2 is 2 (up, down)
d = 2;
