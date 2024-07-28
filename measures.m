%Compute observables of interest


%Trace out the auxialliary system
if m > 1
    SYS = 2; %Tracing over the second index. m. 
    DIM = [nn,m]; %[2 is the dim of the desired matrix, m is the dimension we want to cut]
    DensityE = zeros(length(tt),nn,nn,iP,iK,iQ);
    DensityS = zeros(length(tt),nn,nn,iP,iK,iQ);
    DensityEss = zeros(nn,nn,iP,iK,iQ);
    DensitySss = zeros(nn,nn,iP,iK,iQ);
    for i=1:length(tt)
        for j=1:length(Temps)
            for a=1:length(GA)
                for k=1:length(lam)
                    DensityE(i,:,:,j,a,k) = Ptrace(squeeze(PE(i,:,:,j,a,k)),SYS,DIM);
                    DensityS(i,:,:,j,a,k) = Ptrace(squeeze(PS(i,:,:,j,a,k)),SYS,DIM);

                    DensityEss(:,:,j,a,k) = Ptrace(squeeze(rhoEss1(:,:,j,a,k)),SYS,DIM);
                    DensitySss(:,:,j,a,k) = Ptrace(squeeze(rhoSss1(:,:,j,a,k)),SYS,DIM);
                end
            end 
        end
    end
end


%Calculate the BLP non-markovianity measure
if m > 1
    NblpRC = zeros(length(TT),length(GA),length(lam));
    for j=1:1:length(TT)
        for a=1:1:length(GA)
            for k=1:1:length(lam)
                NblpRC(j,a,k) = BLP(squeeze(log(abs(DensityS(:,1,2,j,a,k)./DensityS(1,1,2,j,a,k)))));
            end
        end 
    end

    %Calculate the RHP non-markovianity measure
    NrhpRC = zeros(length(TT),length(GA),length(lam));
    for j=1:1:length(TT)
        for a=1:1:length(GA)
            for k=1:1:length(lam)
                NrhpRC(j,a,k) = RHP(squeeze(log(abs(DensityS(:,1,2,j,a,k)./DensityS(1,1,2,j,a,k)))));
            end
        end 
    end
end


%Calculate the BLP non-markovianity measure
if m == 1
    NblpWC = zeros(length(TT),length(GA),length(lam));
    for j=1:1:length(TT)
        for a=1:1:length(GA)
            for k=1:1:length(lam)
                NblpWC(j,a,k) = BLP(squeeze(log(abs(PS(:,1,2,j,a,k)./PS(1,1,2,j,a,k)))));
            end
        end 
    end

    %Calculate the RHP non-markovianity measure
    NrhpWC = zeros(length(TT),length(GA),length(lam));
    for j=1:1:length(TT)
        for a=1:1:length(GA)
            for k=1:1:length(lam)
                NrhpWC(j,a,k) = RHP(squeeze(log(abs(PS(:,1,2,j,a,k)./PS(1,1,2,j,a,k)))));
            end
        end 
    end
end