function  [IncoGreenExp,IncoGreenVar]=IncoGreenFunc(par,tra,rec,users)
% This function generates the incoherent part, i.e., stochastic Green's
% Function
 
%% Generate the \overline{\overline{\mathbf{D}}}
term_kR=par.wavenumber.*par.UETX.distance; 
term_2kR=2*term_kR; 
% expectation of D_{xx}
Exp_Dxx=1/(2*par.cavity) .* (sin(term_kR)./term_kR+ ...
    cos(term_kR)./(term_kR.^2)- sin(term_kR)./(term_kR.^3));
% variance of D_{xx}
Var_Dxx= 9/(128*par.cavity^2) .* (4/3 + 2* sin(term_2kR)./term_2kR...
    + 2*cos(term_2kR)./(term_2kR.^2)-2* sin(term_2kR)./(term_2kR.^3))...
    -Exp_Dxx.^2;

Exp_Dyy=Exp_Dxx;
Var_Dyy=Var_Dxx;

Exp_Dzz=1/par.cavity .* (sin(term_kR)./term_kR ...
    - term_kR.* cos(term_kR))./(term_kR.^3); 
Var_Dzz= 3/(16*par.cavity^2) .* (8/15  - 8* sin(term_2kR)./(term_2kR.^3)...
    -24* cos(term_2kR)./(term_2kR.^4)+ 24* sin(term_2kR)./(term_2kR.^5))...
    -Exp_Dzz.^2;

Var_Dxy= 1/(128*par.cavity^2) * (64/15 + 8 * sin(term_2kR)./term_2kR ...
    + 16* cos(term_2kR)./(term_2kR.^2) - 40 * sin(term_2kR) ./ (term_2kR.^3)...
    - 72* cos(term_2kR)./(term_2kR.^4) + 72 *sin(term_2kR)./(term_2kR.^5));


Var_Dxz=1/(8* par.cavity^2) * (4/15 - 2 * cos(term_2kR)./(term_2kR.^2) ...
    + 8 * sin(term_2kR)./(term_2kR.^3) + 18* cos(term_2kR) ./ (term_2kR.^4)...
    -18 * sin(term_2kR) ./(term_2kR.^5) );

IncoGreenExp.XX=Exp_Dxx;
IncoGreenExp.YY=Exp_Dyy;
IncoGreenExp.ZZ=Exp_Dzz;
IncoGreenExp.XY=zeros(users.num*rec.totalNum,tra.totalNum);
IncoGreenExp.XZ=zeros(users.num*rec.totalNum,tra.totalNum);
IncoGreenExp.YZ=zeros(users.num*rec.totalNum,tra.totalNum);


IncoGreenVar.XX=Var_Dxx;
IncoGreenVar.YY=Var_Dyy;
IncoGreenVar.ZZ=Var_Dzz;
IncoGreenVar.XY=Var_Dxy;
IncoGreenVar.XZ=Var_Dxz;
IncoGreenVar.YZ=Var_Dxz;
  

end