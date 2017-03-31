function corranalyse_Jseed_dist_num(L,Jstr,Jdis,Jz,m,Pdist,Jseedmin,Jseedmax)
% function to collect correlation data for multiple seeds
% input: Standard system data L,Jstr,Jdis,Jz,m,Pdist, minimum and maximum
% Jseed values Jseemin, Jseedmax
% output: [distance(j-i), average number of tensors, error]

% Andrew Goldsborough - 01/05/2013

tic

%first sweep
Jseed = Jseedmin;

%open files to read in data
fname = strcat('../spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseed),'_spcorr2.txt');
corr = importdata(fname);

dist = corr(:,2) - corr(:,1);
maxdist = max(dist);

%preallocate
P = zeros(maxdist,1);
Q = zeros(maxdist,1);
G = zeros(maxdist,1);
success = 1;

%import data
for i=1:maxdist
    idx = find(dist==i);
    
    %tnum
    P(i) = sum(corr(idx,4));
    
    %tnum^2
    Q(i) = sum(corr(idx,4).^2);
    
    %degeneracy
    G(i) = size(idx,1);
end

for Jseed=Jseedmin+1:Jseedmax
    
    %open files to read in data
    fname = strcat('../spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseed),'_spcorr2.txt');
    
    try
        corr = importdata(fname);
    catch err
        continue
    end
    
    success = success+1;

    dist = corr(:,2) - corr(:,1);
    
    for i=1:maxdist
        idx = find(dist==i);
        P(i) = P(i) + sum(corr(idx,4));
        Q(i) = Q(i) + sum(corr(idx,4).^2);
        G(i) = G(i) + size(idx,1);
    end
end

%just have the components, need to calculate averages and errors
average = P./G;
error = sqrt((Q./G) - (P./G).^2) ./ sqrt(G);

%print number of successful files
fprintf('%d files successfully imported\n',success);

%open files to write to
fname = strcat('../spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseedmin),'-',num2str(Jseedmax),'_dist_tnum.txt');
fiddisttnum = fopen(fname, 'w');

%print to file
for i=1:maxdist
    fprintf(fiddisttnum,'%d %.15e %.15e\n',i,average(i),error(i));
end

%close file
fclose(fiddisttnum);

toc
