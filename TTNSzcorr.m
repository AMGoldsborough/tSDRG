function [corr,tnum] = TTNSzcorr(L,w,Jorder,site1,site2)
% Two-Point Sz Correlation function for tree tensor network
% input: Length L, Tree tensor network w, order in which the network is
% built Jorder, the sites under investigation site1 < site2
% output: correlator corr, number of tensors in the geodesic tnum

% Andrew Goldsborough 08/04/2013
% function to calculate the two-point Sz correlation function
% C(|i-j|) = < S^{z}_{i} S^{z}_{j} >

%check that site1 is left of site2
if (site1 >= site2)
    fprintf('site 1 has to be to the left of site 2\n');
    
    error('TTNSzcorr:sitepos', 'incompatible site definitions');
end

%create vector showing the active sites
vec = [zeros(site1-1,1);1;zeros(site2-site1-1,1);1;zeros(L-site2,1)];

%create the operators (Sz) as cells with zeros for the empty sites
coro = cell(L,1);
coro{site1} = [0.5 0;0 -0.5];
coro{site2} = [0.5 0;0 -0.5];

%count how mant tensors have been used
tleft = 0;
tright = 0;

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
        coro{pos} = tcon(coro{pos},eye(size(w{i},3)),[-1,-2],[-3,-4]);
        coro{pos} = tcon(coro{pos},w{i},[1,-2,2,-3],[-1,1,2]);
        coro{pos} = tcon(coro{pos},conj(w{i}),[-1,1,2],[-2,1,2]);
        
        %propagate vec
        vec(pos) = 1;
        
        %count tensors
        tleft = tleft + 1;
        
        %remove site
        coro(pos+1) = [];
        vec(pos+1) = [];
        
    elseif vec(pos) == 0
        %right leg used
        coro{pos} = tcon(eye(size(w{i},2)),coro{pos+1},[-1,-2],[-3,-4]);
        coro{pos} = tcon(coro{pos},w{i},[1,-2,2,-3],[-1,1,2]);
        coro{pos} = tcon(coro{pos},conj(w{i}),[-1,1,2],[-2,1,2]);
        
        %propagate vec
        vec(pos) = 1;
        
        %count tensors
        tright = tright + 1;
        
        %remove site
        coro(pos+1) = [];
        vec(pos+1) = [];
        
    else
        %both legs used
        coro{pos} = tcon(coro{pos},coro{pos+1},[-1,-2],[-3,-4]);
        coro{pos} = tcon(coro{pos},w{i},[1,-2,2,-3],[-1,1,2]);
        coro{pos} = tcon(coro{pos},conj(w{i}),[-1,1,2],[-2,1,2]);
    
        %propagate vec
        vec(pos) = 1;
        
        %output the number of tensors before the path is connected
        tnum = tleft + tright + 1;
        
        %remove site
        coro(pos+1) = [];
        vec(pos+1) = [];
    end 
end

%in the end corr should just be the value of the 1 element cell
corr = coro{1};