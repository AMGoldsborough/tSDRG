function TTNSzcorr2(L,w,tL,tR,groundL,groundR,fname)
% Two-Point Sz Correlation function for tree tensor network
% input: Tree tensor network w, tensors below each leg tL & tR, which
% tensor connects to the ground on left or right groundL & groundR, file
% name fname.
% output: correlator corr, number of tensors in the geodesic tnum

% Andrew Goldsborough 04/12/2013
% function to calculate the two-point Sz correlation function
% C(|i-j|) = < S^{z}_{i} S^{z}_{j} >

%version 2 does all at once and saves blocks for use later to improve speed

%need to store the leg for each path left = 0, right =1

%open file for output
fidspcorr = fopen(fname, 'w');

for site1 = 1:L-1
    
    %create the operators as cells with zeros for the empty sites
    coroL = [0.5 0;0 -0.5];

    %initialise stepprev
    stepprev = 1;

    %find the path to the root node
    if groundL(site1) ~= 0 
        prev = groundL(site1);
        legL = 0;
    else
        prev = groundR(site1);
        legL = 1;
    end    

    tpathL = prev;
    
    while prev ~= L-1
        
        next = find(tL==prev);
        
        if size(next,1) ~= 0
            prev = next;
            legL = cat(1,legL,0);
        else
            prev = find(tR==prev);
            legL = cat(1,legL,1);
        end
        
        tpathL = cat(1,tpathL,prev);
    end
    
    for site2 = site1+1:L
        
        %create the spin operator for site 2
        coroR = [0.5 0;0 -0.5];

        %find the path to the root node
        if groundL(site2) ~= 0 
            prev = groundL(site2);
            legR = 0;
        else
            prev = groundR(site2);
            legR = 1;
        end 
        tpathR = prev;
        
        while prev ~= L-1
            
            next = find(tL==prev);
            
            if size(next,1) ~= 0
                prev = next;
                legR = cat(1,legR,0);
            else
                prev = find(tR==prev);
                legR = cat(1,legR,1);
            end
            
            tpathR = cat(1,tpathR,prev);
        end
        
        %find the common node
        [pathint,idxL,idxR] = intersect(tpathL,tpathR);
        cnode = pathint(1);
        
        idxL = idxL(1);
        idxR = idxR(1);

        pathdist = idxL + idxR - 1;

        %contract the left coro up to the common node
        for pathstep = stepprev:idxL-1
            
            %if groundL(tpathL(pathstep)) ~= 0
            if legL(pathstep) == 0
  
                i = tpathL(pathstep);
                
                %left leg used
                coroL = tcon(coroL,w{i},[1,-2],[-1,1,-3]);
                coroL = tcon(coroL,conj(w{i}),[-1,1,2],[-2,1,2]);
            else
                
                i = tpathL(pathstep);
                
                %right leg used
                coroL = tcon(coroL,w{i},[1,-3],[-1,-2,1]);
                coroL = tcon(coroL,conj(w{i}),[-1,1,2],[-2,1,2]);
            end
        end
        
        stepprev = idxL;

        %contract the right coro up to the common node
        for pathstep = 1:idxR-1
            
            %if groundL(tpathR(pathstep)) ~= 0
            if legR(pathstep) == 0
    
                i = tpathR(pathstep);
                
                %left leg used               
                coroR = tcon(coroR,w{i},[1,-2],[-1,1,-3]);
                coroR = tcon(coroR,conj(w{i}),[-1,1,2],[-2,1,2]);
            else
                
                i = tpathR(pathstep);
                
                %right leg used              
                coroR = tcon(coroR,w{i},[1,-3],[-1,-2,1]);
                coroR = tcon(coroR,conj(w{i}),[-1,1,2],[-2,1,2]);
            end
        end

        %contract the common node
        i = cnode;
        
        %both legs used
        coroR = tcon(w{i},coroR,[-1,-2,1],[1,-3,-4]);
        coroR = tcon(coroR,coroL,[-1,1,-3,2],[1,-2,2]);
        coroR = tcon(coroR,conj(w{i}),[-1,1,2],[-2,1,2]);
              
        %now need to contract up to the root node if they are not the same
        if cnode ~= L-1
            for pathstep = idxR+1:size(tpathR,1)
                if legR(pathstep) == 0
                    
                    i = tpathR(pathstep);
                    
                    %left leg used
                    coroR = tcon(coroR,w{i},[1,-2],[-1,1,-3]);
                    coroR = tcon(coroR,conj(w{i}),[-1,1,2],[-2,1,2]);
                else
                    
                    i = tpathR(pathstep);
                    
                    %right leg used
                    coroR = tcon(coroR,w{i},[1,-3],[-1,-2,1]);
                    coroR = tcon(coroR,conj(w{i}),[-1,1,2],[-2,1,2]);
                end
            end
        end
       
        %print to file
        fprintf(fidspcorr,'%d %d %.15e %d\n',site1,site2,coroR,pathdist);
    end
end
 
fclose(fidspcorr);
