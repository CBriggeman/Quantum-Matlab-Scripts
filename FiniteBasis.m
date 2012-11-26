%{
Uses finite basis method to find first four quantum states for a
potential well 6 Angstroms wide with a field of 1 Volt/Angstrom. Uses
symbolic math toolbox. This is based on a problem from "Quantum Mechanics
for Scientists and Engineers" by David A.B. Miller.
%}
clear all
hbar = 1.054571596e-34; %js
me = 9.109381883e-31; %kg
Lx = 6e-10; %meters 
e = 1.602176462e-19; %C
joulestoeV = 6.24150974e18;

%Inifinate well
syms Psi0(n,x);

Psi0(n,x) = sqrt(2/Lx)*sin((n*pi*x)/Lx);

%Potential
syms Hp(x);

Hp(x) = e*1e10*x;

%Finite representation of Hamiltonian
syms Hij(i,j);

Hij(i,j) = (-hbar^2/(2*me))*int(conj(Psi0(i,x))*diff(Psi0(j,x),x,2),x,0,Lx) + ...
    int(conj(Psi0(i,x))*Hp(x)*Psi0(j,x),x,0,Lx);

digits(5);
Hhat = vpa([Hij(1,1) Hij(1,2) Hij(1,3) Hij(1,4); Hij(2,1) Hij(2,2) Hij(2,3) Hij(2,4);...
    Hij(3,1) Hij(3,2) Hij(3,3) Hij(3,4); Hij(4,1) Hij(4,2) Hij(3,4) Hij(4,4)]);
%Sort and process results
[V,D] = eig(Hhat);
D = diag(D);
[C Ind] = sort(D);
Copsi1 = V(:,Ind(1));

%Function space of infinate well states
PsiInf(x) = [Psi0(1,x) Psi0(2,x) Psi0(3,x) Psi0(4,x)]';

%Calculate state and probability functions
Psi1 = sum(Copsi1.*PsiInf(x));
Prob1 = (Psi1)^2;

Copsi2 = V(:,Ind(2));
Psi2 = sum(Copsi2.*PsiInf(x));
Prob2 = (Psi2)^2;

Copsi3 = V(:,Ind(3));
Psi3 = sum(Copsi3.*PsiInf(x));
Prob3 = (Psi3)^2;

Copsi4 = V(:,Ind(4));
Psi4 = sum(Copsi4.*PsiInf(x));
Prob4 = (Psi4)^2;

%Plots probability function of the ground state
h = matlabFunction(Prob1);
fplot(h, [0 Lx])
hold on
xlabel('Position (Meters)')
ylabel('Probability')
title('Probability Wave of An Electron in A Well With Field')
hold off
%h = matlabFunction(Prob2);
%fplot(h,[0 Lx], '--g')

%in eV

E1 = double(D(Ind(1))*joulestoeV);
E2 = double(D(Ind(2))*joulestoeV);
E3 = double(D(Ind(3))*joulestoeV);
E4 = double(D(Ind(4))*joulestoeV);

sprintf('First four eigenenergies are %5f %5f %5f, and %5f eV\n',...
        E1, E2, E3, E4)
