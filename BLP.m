function N = BLP(G) %Input a decoherence function
%First, differentiate it
dG = diff((G));

%Get the indices for the increasing intervals
debut = zeros(length(dG));
fin = zeros(length(dG));
for l=1:1:length(dG)-1
    if dG(l) < 0 && dG(l+1) > 0
        debut(l) = l+1;
    end 

    if dG(l) > 0 && dG(l+1) < 0
        fin(l) = l;
    end
end

%Get rid of all the zeros
debut = nonzeros(debut);
fin = nonzeros(fin);

N = 0; %Initialize the measure
for i=1:1:min(length(debut),length(fin))
    N = N + ((exp(G(fin(i))) - (exp(G(debut(i))))));
end
end