nbath=min(length(bb),length(lam));
%Building the Hamiltonian: 

%Define useful matrices for bulding the Hamiltonian
sigz = [1 0;0 -1]; %Z-pauli matrix
sigx = [0 1;1 0]; %X-pauli matrix
sigy = [0 -1i;1i 0]; %Y-pauli matrix
ata = 0.5*diag([1:2:2*m]); %at*a (number operator in matrix form)
atpa = diag(sqrt([1:1:m-1]),1) + diag(sqrt([1:1:m-1]),-1); %at + a (creation + annihilation) in matrix form

%Spin-Boson 
if Model == 1
    if m == 1 
        HH = 0.5*del*sigz; %Weak-Limit Rotated model;
    end
    
    if m > 1
        HH = 0.5*del*kron(sigz,eye(m)) + la(1)*kron(sigx,atpa) + Om(1)*kron(eye(nn),ata); %'Rotated' Hamiltonian
    end
end

%Dephasing model
if Model == 2
    if m == 1
       HH = 0.5*del*sigz; %Weak-Limit;
    end
    %Dephasing model with RC
    if m > 1
        HH = 0.5*del*kron(sigz,eye(m)) + la(1)*kron(sigz,atpa) + Om(1)*kron(eye(nn),ata); %'Rotated' Hamiltonian

    end
end

%Diagonalize the system hamiltonian
[Uh,Dh]=eig(HH);
    En(1:nn*m^nbath)=diag(Dh);
%=====================================================
%System matrices:
VS=zeros(nn*m^nbath,nn*m^nbath,nbath);
VSd=zeros(nn*m^nbath,nn*m^nbath,nbath);

if Model == 1
    %Spin-Boson
    if m == 1
        VS(:,:,1) = sigx; %system operator 
    end
    if m > 1
        VS(:,:,1) = kron(eye(nn),atpa); %RC System operator
    end
end

if Model == 2
    %Dephasing
    if m == 1
        VS(:,:,1) = sigz; %system operator
    end 
    if m > 1
        %Dephasing with RC
        VS(:,:,1) = kron(eye(nn),atpa); %RC System operator
    end     
end

   for ivs=1:1:nbath
        VSd(:,:,ivs)=Uh'*VS(:,:,ivs)*Uh;
   end