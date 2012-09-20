function yhat = constraintCB(param, model, x, y)
% slack resaling: argmax_y delta(yi, y) (1 + <psi(x,y), w> - <psi(x,yi), w>)
% margin rescaling: argmax_y delta(yi, y) + <psi(x,y), w>

  % the kernel is linear, get a weight vector back
  if size(model.svPatterns, 2) == 0
    w = zeros(size(x)) ;
  else
    w = [model.svPatterns{:}] * (model.alpha .* [model.svLabels{:}]') / 2 ;
  end
  if dot(y*x, w) > 1, yhat = y ; else yhat = -y ; end
  if param.verbose
    fprintf('yhat = violslack([%8.3f,%8.3f], [%8.3f,%8.3f], %3d) = %3d\n', ...
            w, x, y, yhat) ;
  end
end
