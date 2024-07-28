%=====================================
 nB=@(adell,iba) 1./(exp(bb(iba)*adell)-1*(adell~=0)); %neq condition prevents infinities / NaN in evaluation below
 mat=zeros((nn)^2*m^(2*nbath),(nn)^2*m^(2*nbath),nbath);
 coh=zeros((nn)^2*m^(2*nbath),(nn)^2*m^(2*nbath));

 if Spectral_Density==1 %Ohmic sprectral density
     J=@(adell,iba) ga(iba)*adell*exp(-abs(adell)/wc(iba));
     G=@(adell,dell,iba) J(adell,iba)*(nB(adell,iba)+1)*(dell>0) + J(adell,iba)*nB(adell,iba)*(dell<0) + ga(iba)/bb(iba)*(dell==0);
     S=@(adell,dell,iba) 0;
 end
 
 if Spectral_Density==2 %Brownian spectral density
     J=@(adell,iba) (4*(ga(iba)/pi)*Om(iba)^2*la(iba)^2*adell)/((Om(iba)^2 - adell^2)^2 + (2*pi*(ga(iba)/pi)*Om(iba)*adell)^2);
     G=@(adell,dell,iba) pi*J(adell,iba)*(nB(adell,iba)+1)*(dell>0) + pi*J(adell,iba)*nB(adell,iba)*(dell<0) + pi*(4*la(iba)^2*ga(iba)/(pi*bb(iba)*Om(iba)^2))*(dell==0);
     S=@(adell,dell,iba) 0;
 end 
 
 if Spectral_Density==3 %Lorentzian spectral density
     J=@(adell,iba) alpha(ibs)*WC(iba)*adell/(adell^2 + WC(iba)^2);
     G=@(adell,dell,iba) pi*J(adell,iba)*(nB(adell,iba)+1)*(dell>0) + pi*J(adell,iba)*nB(adell,iba)*(dell<0) + pi*(4*la(iba)^2*ga(iba)/(bb(iba)*Om(iba)^2))*(dell==0);
     S=@(adell,dell,iba) 0;
 end 

 
for ii=1:1:nn*m^nbath
    for jj=1:1:nn*m^nbath
        Enn=En(ii)-En(jj);
        coh((ii-1)*(nn*m^nbath)+jj,(ii-1)*(nn*m^nbath)+jj)=-1i*Enn;
        
        %==================================
        % number #1 in Nitzan 10.155
        
        for kk=1:1:nn*m^nbath;
            for ll=1:1:nn*m^nbath;
                dell=En(ll)-En(kk);
                adell=abs(dell);
                for ivs=1:1:nbath
                    R1=VSd(ii,kk,ivs)*VSd(kk,ll,ivs)*(G(adell,dell,ivs)+1i*S(adell,dell,ivs));
                    if Secular==1
                        if ii==jj && ll==jj
                            mat((ii-1)*(nn*m^nbath)+jj,(ll-1)*(nn*m^nbath)+jj,ivs)=-R1 + mat((ii-1)*(nn*m^nbath)+jj,(ll-1)*(nn*m^nbath)+jj,ivs);
                        elseif ii==jj && ll~=jj
                        elseif ii~=jj && ll==jj
                        elseif ii~=jj && ll~=jj
                            mat((ii-1)*(nn*m^nbath)+jj,(ll-1)*(nn*m^nbath)+jj,ivs)=-R1 + mat((ii-1)*(nn*m^nbath)+jj,(ll-1)*(nn*m^nbath)+jj,ivs);
                        end
                    elseif Secular==0
                        mat((ii-1)*(nn*m^nbath)+jj,(ll-1)*(nn*m^nbath)+jj,ivs)=-R1 + mat((ii-1)*(nn*m^nbath)+jj,(ll-1)*(nn*m^nbath)+jj,ivs);
                    end
                end
            end
        end
        
        %===============
        % number #3 in Nitzan 10.155
        for kk=1:1:nn*m^nbath;
            for ll=1:1:nn*m^nbath;
                dell=En(kk)-En(ii);
                adell=abs(dell);
                for ivs=1:1:nbath
                    R2=VSd(ll,jj,ivs)*VSd(ii,kk,ivs)*(G(adell,dell,ivs)+1i*S(adell,dell,ivs));
                    if Secular==1
                        if ii==jj && kk==ll
                            mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs)=R2 + mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs);
                        elseif ii==jj && kk~=ll
                        elseif ii~=jj && kk==ll
                        elseif ii~=jj && kk~=ll
                            mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs)=R2 + mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs);
                        end
                    elseif Secular==0
                        mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs)=R2 + mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs);
                    end
                end
            end
        end
        %===================
        % number #4 in Nitzan
        for kk=1:1:nn*m^nbath;
            for ll=1:1:nn*m^nbath;
                dell=En(ll)-En(jj);
                adell=abs(dell);
                for ivs=1:1:nbath
                    R3=VSd(kk,ii,ivs)*VSd(jj,ll,ivs)*(G(adell,dell,ivs)+1i*S(adell,dell,ivs));
                    if Secular==1
                        if ii==jj && kk==ll
                            mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs)=conj(R3) + mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs);
                        elseif ii==jj && kk~=ll
                        elseif ii~=jj && kk==ll
                        elseif ii~=jj && kk~=ll
                            mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs)=conj(R3) + mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs);
                        end
                    elseif Secular==0
                        mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs)=conj(R3) + mat((ii-1)*(nn*m^nbath)+jj,(kk-1)*(nn*m^nbath)+ll,ivs);
                    end
                end
            end
        end
        %===================
        % number #2 in Nitzan 10.155
        for kk=1:1:nn*m^nbath;
            for ll=1:1:nn*m^nbath;
                dell=En(kk)-En(ll);
                adell=abs(dell);
                for ivs=1:1:nbath
                    R4=VSd(jj,ll,ivs)*VSd(ll,kk,ivs)*(G(adell,dell,ivs)+1i*S(adell,dell,ivs));
                    if Secular==1
                        if ii==jj && ii==kk
                            mat((ii-1)*(nn*m^nbath)+jj,(ii-1)*(nn*m^nbath)+kk,ivs)=-conj(R4) + mat((ii-1)*(nn*m^nbath)+jj,(ii-1)*(nn*m^nbath)+kk,ivs);
                        elseif ii==jj && ii~=kk
                        elseif ii~=jj && ii==kk
                        elseif ii~=jj && ii~=kk
                            mat((ii-1)*(nn*m^nbath)+jj,(ii-1)*(nn*m^nbath)+kk,ivs)=-conj(R4) + mat((ii-1)*(nn*m^nbath)+jj,(ii-1)*(nn*m^nbath)+kk,ivs);
                        end
                    elseif Secular==0
                        mat((ii-1)*(nn*m^nbath)+jj,(ii-1)*(nn*m^nbath)+kk,ivs)=-conj(R4) + mat((ii-1)*(nn*m^nbath)+jj,(ii-1)*(nn*m^nbath)+kk,ivs);
                    end
                end
            end
        end
        %=======================
    end
end
