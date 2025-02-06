function [cavity_ki,CavityEigTerm]=GOE_Eigen(num_M,alpha,par)
% Gausian orthogonal ensemble of M random matrices

diagEle_var=1;
offDiagEle_var=0.5;
term_offMatrix=sqrt(offDiagEle_var)*randn(num_M,num_M);
offMatrix=tril(term_offMatrix,-1)+triu(term_offMatrix',0);
GOEMatrix=offMatrix-diag(diag(offMatrix))+diag(randn(1,num_M)); 
CavityEig=eig(GOEMatrix );
CavityEig=sort(CavityEig,'descend');
NorCavityEig=num_M/pi*asin(CavityEig/sqrt(2*num_M))...
    +CavityEig.*sqrt(2*num_M-CavityEig.^2)/(2*pi);
cavity_ki=par.wavenumber^2-NorCavityEig*2*pi^2/(par.wavenumber*par.cavity);
cavity_ki=sqrt(cavity_ki);

lamdaRow=1./(NorCavityEig-1i*alpha);
sumNorLambda=sum(lamdaRow);
CavityEigTerm=par.wavenumber*par.cavity/(2*pi^2)*sumNorLambda;


end