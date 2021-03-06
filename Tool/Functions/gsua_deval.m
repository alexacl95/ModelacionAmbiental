function y = gsua_deval(par,Table2,xdata)
try
    Table=Table2.Properties.CustomProperties;
catch
    Table=load('ATable.mat');
    Table=Table.Table2;
end
kind=Table.Kind;   
  
if kind==1
    [y,~]=sens_montecarlo(Table2,par,xdata,false,[0,size(par,1)]);
    if size(y,1)==1
        y=squeeze(y)';
    end
    return
end
ni=Table.NumVars;
fun=Table.Solver; 
out=Table.output;
try
switch kind
    case {2,6}
        y=fun(par);
    case 3
        y=deval(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);
    case 4
        y=gsua_intrp(fun(par(:,1:ni),par(:,ni+1:end)),xdata,out);
    case 5
        domain=Table.Domain;
        opt=Table.copt;
        y=gsua_intrp(fun(par,domain,opt),xdata,out,2);
end
catch ME
    warning(ME.message)
    y=inf(length(out),max(size(xdata))); 
end
end