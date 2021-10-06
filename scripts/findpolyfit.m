function [pred,bestOrder,stats] = findpolyfit(x,y,maxPoly)
% Fit up to 'maxPoly' order polynomials, minimising BIC (using mean-squared error term)

% fit all polynomials
aic = nan(maxPoly,1);
bic = nan(maxPoly,1);
params = cell(1,maxPoly);
for i = 1:maxPoly
    fit = polyfit(x,y,i);
    vals = polyval(fit,x);
    mse = sqrt(mean((y-vals).^2));
    LL = log(mse);
    aic(i,1) = length(y)*LL + 2*i;
    bic(i,1) = length(y)*LL + log(length(y))*i;
    params{i} = fit;
end

stats = [];
stats.aic = aic;
stats.bic = bic;
stats.params = params;

% use BIC to determine best fit
[bic,bestOrder] = sort(bic);
diffbic = [diff(bic); NaN] > 10;
bestOrder = bestOrder(find(diffbic,1,'first'));
if isempty(bestOrder)
    bestOrder=1;
end

% generate predicted values
pred = polyval(polyfit(x,y,bestOrder),x);

end