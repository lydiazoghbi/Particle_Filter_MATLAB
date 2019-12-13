function data = solverFunSingle(point,thres)
data = log(thres./point(1,:))./point(2,:);
end