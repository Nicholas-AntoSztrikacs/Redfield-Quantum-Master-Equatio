tic
%Decoherence function exact expression
%Parameters
del = 1; %TLS splitting
la = 1*del*[1]; %RC coupling
Om = 3*del*[1]; %RC frequency
ga = 0.071*[1]; %Residual coupling (assumed small)
T =  0.5*del*[1]; %Temperature of the bath

%Pure state
rho01 = 0.5;
rho00 = 0.5;

bound = 10; %upper limit of frequency integral (lowered it to save time)

TT = [0.5,3];
GA = [0.01,0.1,1,5,10,20];
lam = [0.1,1];
BB = 1./TT;
tt = [0:0.025:25];

%Compute the decoherence function for the pure decoherence model
Deco = zeros(length(tt),length(TT),length(GA),length(lam));
Times = zeros(1001,length(GA));
for pp=1:1:length(TT)
    for aa=1:length(GA)
        for rr=1:1:length(lam)
            
            if lam(rr) == 1
                tt = [0:0.025:25];
            end 
            if lam(rr) == 0.1
                tt = [0:1:1000];
            end 
            
                for qq=1:1:length(tt)
                    JJ=@(w) -4.*(4.*GA(aa)/pi.*Om.^2.*lam(rr).^2.*w).*(coth(BB(pp).*w./2)).*(1 - cos(w*tt(qq)))./(w.^2 .* ((Om.^2 - w.^2)^2 + (2*pi.*GA(aa)/pi.*Om.*w).^2));                     
                    Deco(qq,pp,aa,rr) = integral(JJ,0,bound,'ArrayValued',true);
                end
        end
    end 
end

%Calculate coherences
coherence = abs(rho01.*exp(Deco));

%Calculate the Purity of this state (Analytic expression) 
Purity = 2*abs(rho01.*exp(Deco)).^2 + 2*rho00^2 - 2*rho00 + 1; 

%Calculate the trace distance (Analytic expression)
Dist = rho01*abs(exp(Deco) - exp(Deco(length(tt)))); 

%Calculate the BLP non-markovianity measure
Nblp = zeros(length(TT),length(GA),length(lam));
for j=1:1:length(TT)
    for a=1:1:length(GA)
        for k=1:1:length(lam)
            Nblp(j,a,k) = BLP(squeeze((Deco(:,j,a,k))));
        end
    end 
end

%Calculate the RHP non-markovianity measure
Nrhp = zeros(length(TT),length(GA),length(lam));
for j=1:1:length(TT)
    for a=1:1:length(GA)
        for k=1:1:length(lam)
            Nrhp(j,a,k) = RHP(squeeze((Deco(:,j,a,k))));
        end
    end 
end
toc

%Plot figures
figure(1)
plot(tt,coherence(:,1))
xlabel('$\Delta t$','interpreter','latex')
ylabel('|\rho_{01}|')

figure(2)
plot(tt,Deco(:,1))
xlabel('$\Delta t$','interpreter','latex')
ylabel('$\Gamma(t)$','interpreter','latex')
