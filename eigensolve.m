function sweep=eigensolve(klis,creal,cimag,model,modelparams,tol,itmax,VERBOSE,sweep)

% sweep=eigensolve(klis,creal,cimag,[model,modelparams,tol,itmax,VERBOSE,sweep])
%
% Calculates the solutions to the eigenvalue problem formulated in Zurita-Gotor and Held (2025)
%
% Input arguments:
% ----------------
%  klis:          list of wavenumbers for which the calculation is performed (a vector)
%  creal:         range of real phase speeds to explore (a vector)
%  cimag:         range of imaginary phase speeds to explore (a vector)
%
% Ideally, the spacing of the vectors creal and cimag should be as fine as the spacing between the real eigenvalues, or
% otherwise some solutions may be lost. The required spacing can be determined by trial and error.
%
% The following input arguments are optional:
%  model:         function calculating the cost for the specified model (a function handle, default: @eadymodel)
%  modelparams:   model parameters (a list of cells, default: {}). A model may not have parameters, or allow defaults
%  tol:           cost function tolerance required for convergence (a real number, default: 1E-5)
%  itmax:         maximum number of refining iterations to achieve convergence (an integer, default: 50)
%  VERBOSE:       level of verbosity in the output (0:no output, 1:progress info only (default), 2:full info during iterations)
%  sweep:         Allows to augment the results of a previous calculation when passed as input (a list of datapoints)
%
% Output argument:
% ----------------
% The function returns as only argument a list of datapoints (sweep), where each of these datapoints correspond to
% one of the wavenumbers provided in klis. The datapoints sweep(:) are structures consisting of the following fields: 
%  sweep(:).k:      wavenumber (a real number)
%  sweep(:).nmodes: number of modes found (an integer)
%  sweep(:).eigen:  eigenvalues of the modes (a list of complex numbers, with length nmodes)
%  sweep(:).cost:   converged cost function for the mode calculation (a list of real numbers, with length nmodes)
%  sweep(:).sol:    eigenmodes (a list of modes, one for each eigenvalue). These are structures with the fields:
%                         sweep(:).sol(:).y: latitude (a list of real values)
%                         sweep(:).sol(:).u: zonal velocity (a list of complex values)
%                         sweep(:).sol(:).v: meridional velocity (a list of complex values)
%                         sweep(:).sol(:).z: geopotential (a list of complex values)
%
% Example of use:
% ---------------
% The following call produces the dispersion relation (at reduced resolution) for the control simulation in the paper
%
% sweep=eigensolve([0.01:0.01:1],[-4:0.1:8],[-2:0.05:2],@eadymodel,{3.5,-0.5});
%
% The results from this calculation are also provided in the file control-lr.mat
%
% After running this command or loading the data, we can plot the phase speed as follows:
%
%    for i=1:length(sweep); 
%       for k=1:sweep(i).nmodes; 
%          plot(sweep(i).k,real(sweep(i).eigen(k)),'k.','MarkerSize',3); hold on; 
%    end; end
%
% and we can contour the geopotential for the most unstable mode at k=0.2 as follows:
%
%    k=20; [~,i]=max(imag(sweep(k).eigen));  %indices of most unstable mode with k=0.2  
%    x=linspace(0,2*pi);  
%    y=sweep(k).sol(i).y;    
%    z=sweep(k).sol(i).z;
%    contourf(x,y,real(z*exp(sqrt(-1)*x)),15);

% Read parameters and set defaults
if nargin<3 error('Not enough input arguments in call to function eigensolve'); end
if ~exist('VERBOSE') VERBOSE=1; disp('Using basic verbosity by default'); end
if ~exist('model')  model=@eadymodel; disp('Using Eady model by default'); end
if ~exist('modelparams') modelparams={}; disp('Using the default model parameters'); end
if ~exist('tol') tol=1e-5; 
  if(VERBOSE>0) disp('Default tolerance (1E-5) used'); end; end
if ~exist('itmax') itmax=50; 
  if(VERBOSE>0) disp('Default maximum number of iterations allowed (50) used'); end; end  
if(exist('sweep')) ii=length(sweep)+1; else ii=1; end

for i=1:length(klis); K=klis(i); 
    sweep(ii).k=K;
    nmodes=0;
    cost=costfunction(creal,cimag,K,model,modelparams);
    M=findlocalmin(abs(cost));   
    for k=1:size(M,1)
       crealfine=linspace(creal(M(k,1)-1),creal(M(k,1)+1),5);
       cimagfine=linspace(cimag(M(k,2)-1),cimag(M(k,2)+1),5);
       [cr ,ci]=refinesol(crealfine,cimagfine,K,model,modelparams,tol,itmax,VERBOSE);       
       for j=1:length(cr)
           if(isnan(cr(j))) continue; end   %Not converged
           nmodes=nmodes+1;
           c=cr(j)+sqrt(-1)*ci(j);
           sweep(ii).eigen(nmodes)=c;
           [sweep(ii).cost(nmodes),sweep(ii).sol(nmodes)]=model(K,c,modelparams);
       end
    end
    sweep(ii).nmodes=nmodes;
    if(VERBOSE>0) disp(i+"/"+length(klis)+"   "+nmodes+" modes found"); end
    ii=ii+1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cr,ci]=refinesol(creal,cimag,K,model,modelparams,tol,itmax,VERBOSE)
% Recursive function to succesively refine the solution into an identified minima

if (itmax==0) cr=NaN; 
  if(VERBOSE>1) warning('No convergence refining. Try increasing itmax'); end; return; end
C=costfunction(creal,cimag,K,model,modelparams);
M=findlocalmin(abs(C));
if isempty(M) cr=NaN; return; end
if(size(M,1)>1)     %More than one minimum after refining
    m=inf; i=1; for k=1:size(M,1); if abs(C(M(k,1),M(k,2)))<m i=k; m=abs(C(M(k,1),M(k,2))); end; end
    M=M(i,:);  if(VERBOSE>1) warning('Possibly skiping some minima. Try a finer phase speed range'); end
end

if(abs(C(M(1),M(2)))<tol)    %A single minimum, already converged
    cr=creal(M(1)); ci=cimag(M(2));
    if abs(ci)< cimag(2)-cimag(1); ci=0; end
    if(VERBOSE>1) disp("Done refining: C="+abs(C(M(1),M(2)))); end
else     %A single minimum, needs refining
    if(VERBOSE>1) disp("Pending "+itmax+" refining steps"); end
    creal=linspace(creal(M(1)-1),creal(M(1)+1),5);
    cimag=linspace(cimag(M(2)-1),cimag(M(2)+1),5);
    [cr,ci]=refinesol(creal,cimag,K,model,modelparams,tol,itmax-1,VERBOSE);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,sol]=costfunction(creal,cimag,K,model,modelparams)
%Evaluate cost function for the given model over the full range of phase speeds

NR=length(creal); NI=length(cimag); 
for i=1:NR; 
    for j=1:NI;     
        c=creal(i)+sqrt(-1)*cimag(j);
        C(i,j)=model(K,c,modelparams);
    end
end        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind=findlocalmin(A)
%Find indices of all local minima (with values smaller than at all 8 adjacent points) in matrix A

[NX, NY]=size(A);
ind=[];
for i=2:NX-1; for j=2:NY-1;
  if A(i,j)>A(i,j-1) continue; end
  if A(i,j)>A(i,j+1) continue; end
  if A(i,j)>A(i-1,j) continue; end
  if A(i,j)>A(i+1,j) continue; end
  if A(i,j)>A(i-1,j-1) continue; end
  if A(i,j)>A(i+1,j-1) continue; end
  if A(i,j)>A(i-1,j+1) continue; end
  if A(i,j)>A(i+1,j+1) continue; end
  ind=[ind; i j]; 
end; end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


