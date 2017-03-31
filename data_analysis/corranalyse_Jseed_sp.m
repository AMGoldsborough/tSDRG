function corranalyse_Jseed_sp(L,Jstr,Jdis,Jz,m,Pdist,Jseedmin,Jseedmax)
% function to collect correlation data for multiple seeds
% input: Standard system data L,Jstr,Jdis,Jz,m,Pdist, minimum and maximum
% Jseed values Jseemin, Jseedmax
% output:
% spdist - [distance(j-i), average spcorr, error, number of correlations]
% sptnum - [number of tensors, average spcorr, error, number of correlations]

% Andrew Goldsborough - 01/05/2013

% staggered correlation

tic

%first sweep does no averaging
Jseed = Jseedmin;

%open files to read in data
fname = strcat('../spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseed),'_spcorr2.txt');
corr = importdata(fname);

maxnum = max(corr(:,4));
maxnummax = maxnum;
dist = corr(:,2) - corr(:,1);
maxdist = max(dist);
corr(:,3) = (-1).^(dist) .* corr(:,3);

%preallocate
Pd = zeros(maxdist,1);
Qd = zeros(maxdist,1);
Gd = zeros(maxdist,1);
Pn = zeros(maxnum,1);
Qn = zeros(maxnum,1);
Gn = zeros(maxnum,1);
success = 0;

%import data
for i=1:maxdist
    idx = find(dist==i);
    Pd(i) = sum(corr(idx,3));
    Qd(i) = sum(corr(idx,3).^2);
    Gd(i) = size(idx,1);
end

for i=1:maxnum
    idx = find(corr(:,4)==i);
    if size(idx,1) ~= 0
        Pn(i) = sum(corr(idx,3));
        Qn(i) = sum(corr(idx,3).^2);
    else
        Pn(i) = 0;
        Qn(i) = 0;
    end
    Gn(i) = size(idx,1);
end

success = success + 1;

for Jseed=Jseedmin+1:Jseedmax
    
    %open files to read in data
    fname = strcat('../spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseed),'_spcorr2.txt');
    try
        corr = importdata(fname);
    catch err
        continue
    end
    
    success = success + 1;
    
    maxnum = max(corr(:,4));
    dist = corr(:,2) - corr(:,1);
    maxdist = max(dist);
    corr(:,3) = (-1).^(dist) .* corr(:,3);
    
     for i=1:maxdist
        idx = find(dist==i);
        Pd(i) = Pd(i) + sum(corr(idx,3));
        Qd(i) = Qd(i) + sum(corr(idx,3).^2);
        Gd(i) = Gd(i) + size(idx,1);
    end
    
    %maxnum will not necessarily have the same size each time
    if maxnum > maxnummax
        %average the ones that are there
        for i=1:maxnummax
            idx = find(corr(:,4)==i);
            if size(idx,1) ~= 0
                Pn(i) = Pn(i) + sum(corr(idx,3));
                Qn(i) = Qn(i) + sum(corr(idx,3).^2);
            else 
                Pn(i) = Pn(i);
                Qn(i) = Qn(i);
            end
            Gn(i) = Gn(i) + size(idx,1);
        end
        
        %add the new ones
        for i=maxnummax+1:maxnum
            idx = find(corr(:,4)==i);
            if size(idx,1) ~= 0
                Pn(i) = sum(corr(idx,3));
                Qn(i) = sum(corr(idx,3).^2);
            else
                Pn(i) = 0;
                Qn(i) = 0;
            end
            Gn(i) = size(idx,1);
        end
        
        %update maxnummax
        maxnummax = maxnum;
    else
        
        for i=1:maxnum
            idx = find(corr(:,4)==i);
            if size(idx,1) ~= 0
                Pn(i) = Pn(i) + sum(corr(idx,3));
                Qn(i) = Qn(i) + sum(corr(idx,3).^2);
            else 
                Pn(i) = Pn(i);
                Qn(i) = Qn(i);
            end
            Gn(i) = Gn(i) + size(idx,1);
        end
    end
end

%just have the components, need to calculate averages and errors
corrdist = zeros(maxdist,3);
corrdist(:,1) = Pd./Gd;
corrdist(:,2) = sqrt((Qd./Gd) - (Pd./Gd).^2) ./ sqrt(Gd);
corrdist(:,3) = Gd;

clear Pd;
clear Qd;
clear Gd;

corrtnum = zeros(maxnummax,3);
corrtnum(:,1) = Pn./Gn;
corrtnum(:,2) = sqrt((Qn./Gn) - (Pn./Gn).^2) ./ sqrt(Gn);
corrtnum(:,3) = Gn;

clear Pn;
clear Qn;
clear Gn;

%open files to write to
fname = strcat('../spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseedmin),'-',num2str(Jseedmax),'_spdist.txt');
fidspdist = fopen(fname, 'w');

%print number of successful files
fprintf('%d files successfully imported\n',success);

%print to file
for i=1:maxdist
    fprintf(fidspdist,'%d %.15e %.15e %d\n',i,corrdist(i,1),corrdist(i,2),corrdist(i,3));
end

%close file
fclose(fidspdist);

%open file
fname = strcat('../spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseedmin),'-',num2str(Jseedmax),'_sptnum.txt');
fidsptnum = fopen(fname, 'w');

%print to file
for i=1:maxnummax
    fprintf(fidsptnum,'%d %.15e %.15e %d\n',i,corrtnum(i,1),corrtnum(i,2),corrtnum(i,3));
end

%close file
fclose(fidsptnum);

toc
