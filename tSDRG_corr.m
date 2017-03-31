function tSDRG_corr(L,Jstr,Jdis,Jz,m,Pdist,Jseed)
% Tensor network SDRG algroithm calculating correlation functions
% Keeps only full SU2 blocks

% Andrew Goldsborough 23/04/2013

%inputs:
%L = length of chain
%Jstr = coupling strength
%Jdis = coupling disorder
%Jz = anisotropy strenghth
%m = number of states kept
%Pdist = disorder distribution
%Jseed = random number generator seed

%when compiled the command line inputs are strings, convert to numbers
if ischar(L)==1
  L = str2double(L);
end
if ischar(Jstr)==1
  Jstr = str2double(Jstr);
end
if ischar(Jdis)==1
  Jdis = str2double(Jdis);
end
if ischar(Jz)==1
  Jz = str2double(Jz);
end
if ischar(m)==1
  m = str2double(m);
end
if ischar(Pdist)==1
  Pdist = str2double(Pdist);
end
if ischar(Jseed)==1
  Jseed = str2double(Jseed);
end

%print parameters to screen
fprintf('L = %d, Jstr = %f, Jdis = %f, Jz = %f, m = %d, Pdist = %d, Jseed = %d\n',L,Jstr,Jdis,Jz,m,Pdist,Jseed);

%turn off sr to sa warning
warning('off','MATLAB:eigs:SigmaChangedToSA');

%Jstr = 1, Jdis = 2
%actual for L=8 is -1.979894900894407
%actual for L=10 is -2.678266900089179
%actual for L=20 is -7.607040369704038
%actual for L=100 is -48.365676667617251 

%create arrays for interaction strengths where rand is between (0,1)
%set seed, set by clock if negative
if Jseed >=0
    rng(Jseed);
else
    rng('shuffle');
end

%set the probability distribution
if Pdist==1
    %P(K) = 2 theta(K-1/2)
    J = zeros((L-1),1) + Jstr*random('unif',0.5,1,[(L-1),1]);
elseif Pdist==2
    %P(K) = 1
    J = zeros((L-1),1) + Jstr*(rand((L-1),1));
elseif Pdist==7
    %my original disorder 1
    J = zeros((L-1),1) + Jstr + Jstr*Jdis*(rand((L-1),1) - 0.5);
elseif Pdist==8
    %my original disorder 2
    J = zeros((L-1),1) + Jstr + Jdis*(rand((L-1),1) - 0.5);
elseif Pdist==9
    %box distribution of Hikihara AF
    J = zeros((L-1),1) + Jdis*random('unif',0,1,[(L-1),1]);
elseif Pdist==10
    %Laflorencie's infinite disorder distribution
    J = rand(L-1,1).^Jdis;
end

tic

h = zeros(L,1);

%import hamiltonian MPOs "W" and spin dimension d
[W,~] = heishamhalfSD(L,J,Jz,h);

%create a cell for the course graining tensors
w = cell(1,1);

%create array for the order
Jorder = zeros((L-1),1);

%save the original interaction strengths for possible use later
% Ji = J;

%want to know which tensors connect to the the ground at each site
groundL = zeros(L,1);
groundR = zeros(L,1);

%to find grouondL/R need vector of the number of sites represented by each MPO tensor
nofsites = ones(L,1);

%create an array for the tensor on the left and right leg of each
tR = zeros((L-1),1);
tL = zeros((L-1),1);

%create an array to find tL and tR (curt = current tensor below)
curt = zeros(L,1);

%start algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(L-2)
    
    %find the maximum J
    [~,Jmaxpos] = max(abs(J));
    
    %store Jmaxpos
    Jorder(i) = Jmaxpos;
    
    %save tL and tR
    tL(i,1) = curt(Jmaxpos);
    tR(i,1) = curt(Jmaxpos+1);
    
    %update curt
    curt(Jmaxpos) = i;
    curt(Jmaxpos+1) = [];
    
    %block Jmaxpos and Jmaxpos+1
    block = tcon(W{Jmaxpos}(:,:,1,:),W{Jmaxpos+1}(:,:,:,end),[-1,-2,-5,1],[-3,-4,1,-6]);
    
    %create a tensor for course graining ansatz
    leg2 = size(block,1);
    leg3 = size(block,3);
    
    %fuse to create matrix
    block = tfuse(block,[1,-2,1,-3]);
    block = tfuse(block,[-1,2,2]);
    
    %force hermiticity
    block = 0.5*(block+block');
    
    sizeblock = size(block,1);
    
    [v,nrg] = eig(block);
    
    if m < (sizeblock)
        %take only full blocks
        nrgvec = sort(real(diag(nrg)),'ascend');
        diffnrg = diff(nrgvec);
        diffloc = find(diffnrg>=1E-12);
        ltm = find(diffloc<=m);
        blockend = diffloc(ltm(end));
        nrg = nrg(1:blockend,1:blockend);
        v = v(:,1:blockend);
    end
    
    %split set of vectors into isometry
    w{i,1} = tsplit((v'),2,[leg2,leg3]);
    
    %find sizes of matrices
    sizenew = size(nrg,1);
    sizeJ = size(J);
    
    %contract isometries with full two site Hamiltonian
    %faster to do each element individually than explicitly contract
    if Jmaxpos ~= 1
        
        if Jmaxpos ~= sizeJ(1)
            %away form edges
            
            blockfull = zeros(sizenew,sizenew,5,5);
            blockfull(:,:,1,1) = eye(sizenew);
            blockfull(:,:,1,2) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{Jmaxpos+1}(:,:,1,2)))*v;
            blockfull(:,:,1,3) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{Jmaxpos+1}(:,:,1,3)))*v;
            blockfull(:,:,1,4) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{Jmaxpos+1}(:,:,1,4)))*v;
            blockfull(:,:,1,5) = nrg;
            blockfull(:,:,2,5) = v'*(kron(W{Jmaxpos}(:,:,2,5),W{Jmaxpos+1}(:,:,5,5)))*v;
            blockfull(:,:,3,5) = v'*(kron(W{Jmaxpos}(:,:,3,5),W{Jmaxpos+1}(:,:,5,5)))*v;
            blockfull(:,:,4,5) = v'*(kron(W{Jmaxpos}(:,:,4,5),W{Jmaxpos+1}(:,:,5,5)))*v;
            blockfull(:,:,5,5) = eye(sizenew);
        else
            %RHS
            
            blockfull = zeros(sizenew,sizenew,5,5);
            blockfull(:,:,1,5) = nrg;
            blockfull(:,:,2,5) = v'*(kron(W{Jmaxpos}(:,:,2,5),W{Jmaxpos+1}(:,:,5,5)))*v;
            blockfull(:,:,3,5) = v'*(kron(W{Jmaxpos}(:,:,3,5),W{Jmaxpos+1}(:,:,5,5)))*v;
            blockfull(:,:,4,5) = v'*(kron(W{Jmaxpos}(:,:,4,5),W{Jmaxpos+1}(:,:,5,5)))*v;
            blockfull(:,:,5,5) = eye(sizenew);
        end
    else
        %LHS
        
        blockfull = zeros(sizenew,sizenew,5,5);
        blockfull(:,:,1,1) = eye(sizenew);
        blockfull(:,:,1,2) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{Jmaxpos+1}(:,:,1,2)))*v;
        blockfull(:,:,1,3) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{Jmaxpos+1}(:,:,1,3)))*v;
        blockfull(:,:,1,4) = v'*(kron(W{Jmaxpos}(:,:,1,1),W{Jmaxpos+1}(:,:,1,4)))*v;
        blockfull(:,:,1,5) = nrg;
    end
    
    %set this new block as ham for site Jmaxpos
    W{Jmaxpos} = blockfull;
    
    %remove combined site
    W(Jmaxpos+1) = [];
    J(Jmaxpos) = [];
    
    %store ground
    if tL(i) == 0
        groundL(sum(nofsites(1:Jmaxpos))) = i;
    end
    
    if tR(i) == 0
        groundR(sum(nofsites(1:Jmaxpos+1))) = i;
    end
    
    %update the number of sites
    nofsites(Jmaxpos) = nofsites(Jmaxpos)+nofsites(Jmaxpos+1);
    nofsites(Jmaxpos+1) = [];
    
    %need to get the new gaps for each side of the block
%     sizeblock = size(block);
    
    if Jmaxpos ~= 1
        %LHS of block
        
        % H' = HB + HC + HB
        Hc = tcon(W{Jmaxpos-1}(:,:,1,:),W{Jmaxpos}(:,:,:,end),[-1,-2,-5,1],[-3,-4,1,-6]);
        
        %fuse to create matrix
        Hc = tfuse(Hc,[1,-2,1,-3]);
        Hc = tfuse(Hc,[-1,2,2]);
        
        %force hermiticity
        Hc = 0.5*(Hc+Hc');
        gap = eig(Hc);
        gap = sort(gap,'ascend');
        diffgap = diff(gap);
        
        %heighest gap in m
        if size(gap,1) <= m
            loc = find(abs(diffgap)>=1E-12);
        else
            loc = find(abs(diffgap(1:m))>=1E-12);
        end
        
        if size(loc,1) == 0
            %no gap found, exit
            error('SDRGgaptestSU2:nogap','no gap detected');
        else
            gaploc = loc(end);
            
            J(Jmaxpos-1) = diffgap(gaploc);
        end
    end
    
    sizeJ = size(J);
    
    if Jmaxpos ~= sizeJ(1)+1
        %RHS of block
        
        % H' = HB + HC + HB
        Hc = tcon(W{Jmaxpos}(:,:,1,:),W{Jmaxpos+1}(:,:,:,end),[-1,-2,-5,1],[-3,-4,1,-6]);
        
        %fuse to create matrix
        Hc = tfuse(Hc,[1,-2,1,-3]);
        Hc = tfuse(Hc,[-1,2,2]);
        
        %force hermiticity
        Hc = 0.5*(Hc+Hc');
        gap = eig(Hc);
        gap = sort(gap,'ascend');
        diffgap = diff(gap);
        
        %heighest gap in m
        if size(gap,1) <= m
            loc = find(abs(diffgap)>=1E-12);
        else
            loc = find(abs(diffgap(1:m))>=1E-12);
        end
        
        if size(loc,1) == 0
            %no gap found, exit
            error('SDRGgaptestSU2:nogap','no gap detected');
        else
            gaploc = loc(end);
            
            J(Jmaxpos) = diffgap(gaploc);
        end        
    end
end

%final step builds the full system block
i = L-1;
Jmaxpos = 1;
Jorder(i) = Jmaxpos;

%save tL and tR
tL(i,1) = curt(Jmaxpos);
tR(i,1) = curt(Jmaxpos+1);

%remove curt
clear curt;

%block Jmaxpos and Jmaxpos+1
block = tcon(W{Jmaxpos}(:,:,1,:),W{Jmaxpos+1}(:,:,:,end),[-1,-2,-5,1],[-3,-4,1,-6]);

leg2 = size(block,1);
leg3 = size(block,3);
% leg1 = m;

block = tfuse(block,[1,-2,1,-3]);
block = tfuse(block,[-1,2,2]);

%force hermiticity
block = 0.5*(block+block');
[v,nrg] = eig(block);

%take ground state as the top tensor
v = v(:,1);

%split into top tensor
w{i,1} = tsplit(v',2,[leg2,leg3]);
fprintf('gs energy is %.15f\n',nrg(1,1));

%store ground
if tL(i) == 0
    groundL(sum(nofsites(1:Jmaxpos))) = i;
end

if tR(i) == 0
    groundR(sum(nofsites(1:Jmaxpos+1))) = i;
end

%update the number of sites
nofsites(Jmaxpos) = nofsites(Jmaxpos)+nofsites(Jmaxpos+1);
nofsites(Jmaxpos+1) = [];

%expectation values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check norm of state
% fprintf('%f\n',TTNnorm(L,w,Jorder))

%open file to write correlations to    
fprintf('printing correlation functions\n');

%spincorr2 finds spin-spin correlations for all pairs of sites
%output format [i j spincorr(i,j) tnum]
%tnum = number of tensors connecting sites i and j
fname = strcat('./spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseed),'_spcorr2.txt');
TTNspincorr2(L,w,tL,tR,groundL,groundR,fname);

%spincorr finds the spin-spin correlation for a given pair of sites (i,j)
% [corr,tnum] = TTNspincorr(L,w,Jorder,i,j);

toc
