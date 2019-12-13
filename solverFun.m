function data = solverFun(point,thres)
syms t
func = point(1)*exp(point(2)*t)+point(3)*exp(point(4)*t) == thres;
data = vpasolve(func,t);
end