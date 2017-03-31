function [corr,tnum] = TTNspincorr(L,w,Jorder,site1,site2)
% Two-Point spin Correlation function for tree tensor network
% input: Length L, Tree tensor network w, order in which the network is
% built Jorder, the sites under investigation site1 < site2
% output: correlator corr, number of tensors in the geodesic tnum

% Andrew Goldsborough 08/04/2013
% function to calculate the two-point Sz correlation function
% C(|i-j|) = < \vec{S}_{i} . \vec{S}_{j} >

%check that site1 is left of site2
if (site1 >= site2)
    fprintf('site 1 has to be to the left of site 2\n');
    
    error('TTNSzcorr:sitepos', 'incompatible site definitions');
end

%create vector showing the active sites
vec = [zeros(site1-1,1);1;zeros(site2-site1-1,1);1;zeros(L-site2,1)];

%create the operators as cells with zeros for the empty sites
coro = cell(L,1);
coro{site1} = zeros(2,2,3,3);
coro{site1}(:,:,1,1) = [0 0.5;0 0];
coro{site1}(:,:,1,2) = [0 0;0.5 0];
coro{site1}(:,:,1,3) = [0.5 0;0 -0.5];

coro{site2} = zeros(2,2,3,3);
coro{site2}(:,:,1,3) = [0 0;1 0];
coro{site2}(:,:,2,3) = [0 1;0 0];
coro{site2}(:,:,3,3) = [0.5 0;0 -0.5];

%count how many tensors have been used tnum outputs, tnumr is a running total
tnumr = 0;
tnum = 0;

%calculate correlator
for i=1:L-1
    
    pos = Jorder(i);
        
    %check to see whether tensor is used
    if vec(pos) == 0 && vec(pos+1) == 0
        %not used remove site
        coro(pos+1) = [];
        vec(pos+1) = [];
        
    elseif vec(pos+1) == 0
        %left leg used
        eyet = zeros(size(w{i},3),size(w{i},3),3,3);
        eyet(:,:,1,1) = eye(size(w{i},3));
        eyet(:,:,2,2) = eye(size(w{i},3));
        eyet(:,:,3,3) = eye(size(w{i},3));
        
        coro{pos} = tcon(coro{pos},eyet,[-1,-2,-5,1],[-3,-4,1,-6]);
        coro{pos} = tcon(coro{pos},w{i},[1,-2,2,-3,-4,-5],[-1,1,2]);
        coro{pos} = tcon(coro{pos},conj(w{i}),[-1,1,2,-3,-4],[-2,1,2]);
        
        %propagate vec
        vec(pos) = 1;
        
        %count tensors
        tnumr = tnumr + 1;

        %remove site
        coro(pos+1) = [];
        vec(pos+1) = [];
        
    elseif vec(pos) == 0
        %right leg used
        eyet = zeros(size(w{i},2),size(w{i},2),3,3);
        eyet(:,:,1,1) = eye(size(w{i},2));
        eyet(:,:,2,2) = eye(size(w{i},2));
        eyet(:,:,3,3) = eye(size(w{i},2));
        
        coro{pos} = tcon(eyet,coro{pos+1},[-1,-2,-5,1],[-3,-4,1,-6]);
        coro{pos} = tcon(coro{pos},w{i},[1,-2,2,-3,-4,-5],[-1,1,2]);
        coro{pos} = tcon(coro{pos},conj(w{i}),[-1,1,2,-3,-4],[-2,1,2]);
        
        %propagate vec
        vec(pos) = 1;
        
        %count tensors
        tnumr = tnumr + 1;
        
        %remove site
        coro(pos+1) = [];
        vec(pos+1) = [];
        
    else
        %both legs used
        coro{pos} = tcon(coro{pos},coro{pos+1},[-1,-2,-5,1],[-3,-4,1,-6]);
        coro{pos} = tcon(coro{pos},w{i},[1,-2,2,-3,-4,-5],[-1,1,2]);
        coro{pos} = tcon(coro{pos},conj(w{i}),[-1,1,2,-3,-4],[-2,1,2]);
    
        %propagate vec
        vec(pos) = 1;
        
        %output the number of tensors before the path is connected
        tnum = tnumr + 1;
            
        %remove site
        coro(pos+1) = [];
        vec(pos+1) = [];
    end 
end

%in the end corr should just be the value of the 1 element cell
corr = coro{1}(1,1,1,3);