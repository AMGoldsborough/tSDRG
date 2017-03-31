function szt = TTNSzT(L,w,Jorder)
% Total Sz for TTN
% input: Length (L), Tree tensor network (w), order in which the network is
% built (Jorder)
% output: total Sz (szt)

% Andrew Goldsborough 03/05/2013
% function to calculate the total Sz expectation value

%create SzT MPO
Wz = cell(L,1);

Wz{1} = zeros(2,2,2,2);
Wz{1}(:,:,1,1) = [1 0;0 1];
Wz{1}(:,:,1,2) = [0.5 0;0 -0.5];

for i=2:(L-1)
    Wz{i} = zeros(2,2,2,2);
    Wz{i}(:,:,1,1) = [1 0;0 1];
    Wz{i}(:,:,1,2) = [0.5 0;0 -0.5];
    Wz{i}(:,:,2,2) = [1 0;0 1];
end

Wz{L} = zeros(2,2,2,2);
Wz{L}(:,:,1,2) = [0.5 0;0 -0.5];
Wz{L}(:,:,2,2) = [1 0;0 1];

%calculate correlator
for i=1:L-1
    
    pos = Jorder(i);

    Wz{pos} = tcon(Wz{pos},Wz{pos+1},[-1,-2,-5,1],[-3,-4,1,-6]);
    Wz{pos} = tcon(Wz{pos},w{i},[1,-2,2,-3,-4,-5],[-1,1,2]);
    Wz{pos} = tcon(Wz{pos},conj(w{i}),[-1,1,2,-3,-4],[-2,1,2]);
    
    %remove site
    Wz(pos+1) = [];
end

%in the end corr should just be the value of the 1 element cell
szt = Wz{1}(:,:,1,2);



