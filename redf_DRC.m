%Master Equation for Reduced Density Matrix "p" using Redfield Dynamics, as well as energy currents between various baths
tic
clear
%====================================
%Simulation Parameters
tend=100; %number of periods for plots
Spectral_Density=1; %1=ohmic, %2=Brownian, %3=Lorentzian
Model=2; %1=Spin Boson, 2=Pure Dephasing
Secular=0; %1= do secular approx., 0=do full dynamics
%Physical Parameters NOTE LENGTH OF VECTORS TT, wc, ga MUST EQUAL NUMBER OF BATHS, AND LENGTH OF ens MUST EQUAL NUMBER OF SITES

nn=2; %number of sites (2, qubit)
m=5; %number of oscillator levels after truncation (m=1 is no RC)
%=====================================

%=====================================
%Model Parameters
del=1; %Delta term in Hamiltonian (Tunneling rate) 
TT = [0.95]; %Bath temperature
bb = 1./TT; %Inverse temperature
Om = 10*[1]; %RC frequency
wc = pi*[1000]; %Cutoff frequency
ga = 0.0071*[1]; %RC system-bath coupling
la = del*[1]; %RC energy
%======================================

%======================================
%loop parameters
GA = [0.01, 0.1, 1];
lam = [0.1, 1, 5, 10];
Temps = [0.5];
%======================================
for iP=1:1:length(Temps) %do the code over this loop
    for iK=1:1:length(GA)
        for iQ=1:1:length(lam) 
                la = lam(iQ)*[1];
                TT = Temps(iP);
                bb = 1./TT;
                ga = GA(iK);

                Hamiltonian_DRC; %build the Hamiltonian/interaction matrices and express them in the energy basis
                tensor2_DRC; %Calculates Dissipators and coherent evolution
            
                %Dynamics and initial conditions
                p0=zeros(nn*m^nbath,nn*m^nbath);

                %Pure states
                if m == 1                              
                    p0(1,1) = 0.5;
                    p0(2,2) = 0.5;
                    p0(1,2) = 0.5;
                    p0(2,1) = 0.5;                    
                end

                if m > 1
                    %State of the qubit (I can choose this)
                    p = [0.5 0.5;0.5 0.5];
                    %Initialize the RC
                    pRC = zeros(m);
                    %Thermal state of H_RC (canonical)
                    for vv=1:1:m
                       pRC(vv,vv) = exp(-vv*0.5*Om*bb); %Thermal state of H_RC 
                    end
                    pRC = 1./trace(pRC) * pRC;     

                    p0 = kron(p,pRC); %Combine the system and RC
                end
                
                p0d=Uh'*p0*Uh; %define system initial conditions in the eigenbasis %MK double-check
                pv0=reshape(p0d,[nn^2*m^(2*nbath),1]);
                
                MAT=coh;
                for ivs=1:1:nbath
                    MAT=MAT+mat(:,:,ivs);
                end
                [UM,DM]=eig(MAT);
                
                per=tend/1000; % time step (Dynamical feature)
                dt=per;
                tt=0:dt:tend; %length of time for graphs (small effect on time to evaluate)
                      
                pv=@(It) UM*expm(DM*It)/UM*pv0; %vector of density matrix in energy basis
                rhoE=@(It) reshape(pv(It),[nn*m^nbath,nn*m^nbath]); %back-transform to site basis, square matrix
                rhoS=@(It) Uh*reshape(pv(It),[nn*m^nbath,nn*m^nbath])*Uh'; %back-transform to site basis, square matrix
                
                %time evolve
                for it=1:1:length(tt)
                    PE(it,:,:,iP,iK,iQ)=rhoE(tt(it));
                    PS(it,:,:,iP,iK,iQ)=rhoS(tt(it));
                end
                
                %Only useful to get steady state coherences
                SS_current_DRC %script to compute steady-state outputs 
        end
    end 
end 

%Compute final system state and measures of non-Markovianity
measures;

toc