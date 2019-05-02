function dertax = Dirac(phi, efso);

dertax = (1/3.1415)*efso*(1./(efso.*efso + phi.*phi));

end