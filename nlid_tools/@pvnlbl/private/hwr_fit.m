function [alpha,yhat,u,Gu_n,y] = hwr_fit(uMin,uMax,slope,threshold,order,polynomType,plotFlag)

u = (uMin:0.0001:uMax)';
u_n = u/uMax;

y = slope*(max(u,threshold) - threshold);

switch polynomType
    case 'power'
        aPolynom = polyfit(u,y,order);
        yhat = polyval(aPolynom,u);
        alpha = aPolynom; %.coeff;
        Gu_n = '';
    case 'chebychev'
        Gu_n = multi_tcheb(u_n,order);
        alpha = lscov(Gu_n,y);
        yhat = Gu_n*alpha;
end

if plotFlag
    figure;
    plot(u,y,u,yhat,'r'); xlabel('Input'); ylabel('NL Output')
    legend('HWR','Fitted Polynom')
end


