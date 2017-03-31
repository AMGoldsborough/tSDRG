function eeanalyse_nb_average(L,Jstr,Jdis,Jz,m,Pdist,Jseedmin,Jseedmax)
%function eeanalyse_nb_average(L,Jstr,Jdis,Jz,m,Pdist,Jseedmin,Jseedmax)
%
% function to average over the ee and minimal surface data for multiple seeds
% central blocks rather than bipartitions

% Andrew Goldsborough - 01/11/2013

tic

%open files to write to
fname = strcat('../ee/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseedmin),'-',num2str(Jseedmax),'_eeb_av.txt');
fidee = fopen(fname, 'w');

%preallocate
Pee = zeros((L/2)-1,1);
Qee = zeros((L/2)-1,1);
P = zeros((L/2)-1,1);
Q = zeros((L/2)-1,1);
PX = zeros((L/2)-1,1);
QX = zeros((L/2)-1,1);
G = 0;

% get data from files
for Jseed=Jseedmin:Jseedmax
    
    %open files to read in data
    fname = strcat('../ee/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseed),'_eeb.txt');
    try
        n_A = importdata(fname);
    catch err
        continue
    end
    
    if size(n_A,1) ~= (L/2)-1
        continue
    end
    
    G = G + 1;
    
    %get data from ee
    Pee = Pee + n_A(:,2);
    Qee = Qee + n_A(:,2).^2;
    
    P = P + n_A(:,3);
    Q = Q + n_A(:,3).^2;
    
    PX = PX + n_A(:,4);
    QX = QX + n_A(:,4).^2;
    
end

%calculate average and standard error
ee = Pee./G;
ee_error = sqrt((Qee./G) - (Pee./G).^2) ./ sqrt(G);

n_A = P./G;
na_error = sqrt((Q./G) - (P./G).^2) ./ sqrt(G);

chi = PX./G;
chi_error = sqrt((QX./G) - (PX./G).^2) ./ sqrt(G);

%print number of successful files
fprintf('%d files successfully imported\n',G);

%print to file
for i=1:(L/2)-1
    fprintf(fidee,'%d %.15e %.15e %.15e %.15e %.15e %.15e\n',2*i,ee(i),ee_error(i),n_A(i),na_error(i),chi(i),chi_error(i));
end

%close file
fclose(fidee);

toc
