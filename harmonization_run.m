function normdata=harmonization_run(data,batch,covariates,outpath,method)
if nargin<2
fprintf('--------harmonize multi-batches(centers) data based on Combat---------\n');
fprintf('Format: normdata=harmonization_run(data,batch,covariates,outpath,method)\n');
fprintf('Input:\n--data:\tdata need harmonization,numeric matrix [Nvariable,Nsample]\n');
fprintf('--batch:\tbatch (center) effects that need be removed,vector[1,Nsample]\n');
fprintf('--covariates:\tcovariates of Biological information that need preserved,[Ncovariates,Nsample]\n');
fprintf('--outpath:\toutput file path,can be .mat or csv format\n');
fprintf('--method:\t ''parametric''(default) or ''non-parametric''\n');
fprintf('-------Written by QinWen 20200707--------\n');
fprintf('This script depend on ''ComBatHarmonization''(by Jfortin1) at github:\nhttps://github.com/Jfortin1/ComBatHarmonization/tree/master/Matlab\n');
return
end
if nargin<3;covariates=[];end
if nargin<4;outpath='';end
if nargin<5;method='parametric';end
if strcmp(method,'parametric');method=1;end
[nvar,nsamp]=size(data);
[m,n]=size(batch);
[x,y]=size(covariates);

if n~=nsamp&&m==nsamp
    batch=batch';
    [m,n]=size(batch);
end
if n~=nsamp
    error('The dimension of batch [%d %d] is not equal to that of the data[%d %d]\n',m,n,nvar,nsamp);    
end

if x~=nsamp&&y==nsamp
    covariates=covariates';
    [x,y]=size(covariates);
end
if x~=nsamp
    if isempty(covariates)
        covariates=[];
    else
    error('The dimension of covariates [%d %d] is not equal to that of the data[%d %d]\n',x,y,nvar,nsamp);   
    end
end
    
normdata = combat(data,batch,covariates,method);

if ~isempty(outpath)
[~,~,ext]=fileparts(outpath);
switch ext
    case '.csv'
        csvwrite(outpath,normdata);
    case '.mat'
        save(outpath,'normdata');
    otherwise
        return
end
end
