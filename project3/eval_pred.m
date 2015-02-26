function [confidence, rank] = eval_pred(likelihoods)
%output a metrics for confidence
[vals, rank] = sort(likelihoods, 'descend');
confidence = 1/(vals(1)/vals(2));

end

