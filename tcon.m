function [T] = tcon(R,S,Rid,Sid)
% Tensor contraction function. Contracts specified indices of tensors R and
% S to give T.
% input: tensors R with indices Rid and S with indices Sid
% Rid & Sid: vectors of tensor indices labelled with positive numbers to be
% contracted and negative numbers giving the desired order of the output
% indices
% output: tensor T with indices contracted as set by Rid and Sid

% tcon - Andrew Goldsborough 04/07/2012
% based on chapter IIB in arXiv:1008.4774v1 by Singh, Pfeifer and Vidal

%note: does not support contraction of two legs of the same tensor. For
%that contract with an identity .

%tic
%check that there are no zero indices
zeroidxR = find(Rid==0);
zeroidxS = find(Sid==0);

if size(zeroidxR,2) ~= 0 || size(zeroidxS,2) ~= 0
    fprintf('index cannot equal zero\n');
        
    error('tcon:zeroidx', 'zero index detected');
end

%sort in increasing order keeping the ordering index array
[Rid,Ridoix] = sort(Rid);
[Sid,Sidoix] = sort(Sid);

%permute indices st the indices to be contracted are on the right
R = permute(R,Ridoix);
S = permute(S,Sidoix);

%get size of arrays
sizeR = size(R);
sizeS = size(S);
 
%need the number of negative indices
numnegR = size(find(Rid<0));
numnegS = size(find(Sid<0));

%build array of T index sizes
Tidxsize = cat(2,sizeR(1:numnegR(2)),sizeS(1:numnegS(2)));

%build Tid array
Tid = cat(2,Rid(1:numnegR(2)),Sid(1:numnegS(2)));
sizeTid = size(Tid);

%check that there are the same number of contracted indices are equal for
%the two tensors. Note that in general this shouldn't be disallowed, but I
%dont think this routine can contract two legs on the same tensor
if (size(Rid,2) - numnegR(2)) ~= (size(Sid,2) - numnegS(2))
    fprintf('the number of contracted indices should be equal for both tensors\n');
    
    error('tcon:numneg', 'unequal number of positive numbers on each tensor detected');
end

%check that the size of the contracted indices match (does not work for
%null indices, need to rewrite)
% for i=1:(sizeRid(2) - numnegR-1)
%     if sizeR(sizeRid(2)-i+1) ~= sizeS(sizeSid(2)-i+1)
%         fprintf('the size of the contracted indices must be equal\n');
%         
%         error('tcon:idxsize', 'size of contracted index does not match');
%     end    
% end

%check that the uncontracted indices are not repeated
if (sizeTid(2) ~= 0) && (isequal(sort(Tid), unique(Tid)) ~= 1)
    fprintf('uncontracted index numbers must be unique\n');
    error('tcon:Tidxunique','repeated uncontracted index detected');
end

%need the size of the joined indices
sizeu = prod(sizeR(1:numnegR(2)));
sizev = prod(sizeR((numnegR(2)+1):end));
sizew = prod(sizeS(1:numnegS(2)));

%reshape: join 'indices to be contracted' and those not
R = reshape(R,sizeu,sizev);
S = reshape(S,sizew,sizev);

%multiply the matrices noting that S is transposed
T = R*S.';
clear R;
clear S;

%if not fully contracted
if sizeTid(2) ~= 0
    
    %split non contracted indices
    T = reshape(T,Tidxsize);
    
    %permute to desired form
    [~,Tidoix] = sort(Tid,'descend');
    T = permute(T,Tidoix);
else
    %T=T;
end
%toc


