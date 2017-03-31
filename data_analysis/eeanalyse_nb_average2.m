function eeanalyse_nb_average2(L,Jstr,Jdis,Jz,m,Pdist,Jseedmin,Jseedmax)
%function eeanalyse_nb_average(L,Jstr,Jdis,Jz,m,Pdist,Jseedmin,Jseedmax)
%
% function to average over the minimal surface data for multiple seeds for
% central blocks rather than bisections
%
% version 2 doesn't discard non complete files, just extracts what can be used

% Andrew Goldsborough - 15/11/2013

tic

%open files to write to
fname = strcat('../ee/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseedmin),'-',num2str(Jseedmax),'_eeb_av2.txt');
fidee = fopen(fname, 'w');

%preallocate
Pee = zeros((L/2)-1,1);
Qee = zeros((L/2)-1,1);
Pna = zeros((L/2)-1,1);
Qna = zeros((L/2)-1,1);
PX = zeros((L/2)-1,1);
QX = zeros((L/2)-1,1);
G = zeros((L/2)-1,1);
success = 0;

% get data from files
for Jseed=Jseedmin:Jseedmax
    
    %open files to read in data
    fname = strcat('../ee/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(m),'_',num2str(Pdist),'_',num2str(Jseed),'_eeb.txt');
    try
        eefile = importdata(fname);
    catch err
        continue
    end
    
    %get data from ee file
    for i=1:size(eefile,1)
        Pee(eefile(i,1)/2,1) = Pee(eefile(i,1)/2,1) + eefile(i,2);
        Qee(eefile(i,1)/2,1) = Qee(eefile(i,1)/2,1) + eefile(i,2)^2;
        
        Pna(eefile(i,1)/2,1) = Pna(eefile(i,1)/2,1) + eefile(i,3);
        Qna(eefile(i,1)/2,1) = Qna(eefile(i,1)/2,1) + eefile(i,3)^2;
        
        PX(eefile(i,1)/2,1) = PX(eefile(i,1)/2,1) + eefile(i,4);
        QX(eefile(i,1)/2,1) = QX(eefile(i,1)/2,1) + eefile(i,4)^2;
        
        G(eefile(i,1)/2,1) = G(eefile(i,1)/2,1) + 1;
    end
    
    success = success + 1;
end

%calculate average and standard error
ee = Pee./G;
ee_error = sqrt((Qee./G) - (Pee./G).^2) ./ sqrt(G);

n_A = Pna./G;
na_error = sqrt((Qna./G) - (Pna./G).^2) ./ sqrt(G);

chi = PX./G;
chi_error = sqrt((QX./G) - (PX./G).^2) ./ sqrt(G);

%print number of successful files
fprintf('%d files successfully imported\n',success);

%print to file
for i=1:(L/2)-1
    fprintf(fidee,'%d %.15e %.15e %.15e %.15e %.15e %.15e %d\n',2*i,ee(i),ee_error(i),n_A(i),na_error(i),chi(i),chi_error(i),G(i));
end

%close file
fclose(fidee);

toc