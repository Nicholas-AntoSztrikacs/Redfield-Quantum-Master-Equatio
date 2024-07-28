%compute the steady state density, then current
rhoE0=zeros((nn)^2*m^(2*nbath),1);
rhoE0((nn)^2*m^(2*nbath))=1;
MM=MAT;
MM((nn)^2*m^(2*nbath),1:(nn)^2*m^(2*nbath))=0;

for ic1=1:1:nn*m^nbath
    for ic2=1:1:nn*m^nbath
        if ic1==ic2
            MM((nn)^2*m^(2*nbath),(ic1-1)*(nn*m^nbath)+ic2)=1;
        end
    end
end
rhoEss=reshape(pinv(MM)*rhoE0,nn*m^nbath,nn*m^nbath); %steady-state density in energy basis
rhoEss2(:,:,iP,iK,iQ)=permute(PE(end,:,:,iP,iK,iQ),[2,3,1,4]);
rhoEss1(:,:,iP,iK,iQ)=rhoEss;
rhoSss=Uh*rhoEss*Uh'; %steady-state density in system basis
rhoSss1(:,:,iP,iK,iQ)=rhoSss;

%current from each bath
for iv1=1:1:nbath
    jq(iv1,iP,iK,iQ)=trace(reshape(mat(:,:,iv1)*reshape(rhoEss,(nn)^2*m^(2*nbath),1),(nn)*m^nbath,(nn)*m^nbath)*Dh); %steady-state heat current per bath
end
