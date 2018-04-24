clear all;
close all;
clc;

%% fitting method 拟合方法
% step 1. gradient descent for all parameters in one iter
% step 2. change one parameter at a time when backtracking cannot find a
%   step. This is because simulator might suffer from accuracy problem
% step 3. random fitting when step 2 is done (try 10000 points around the opt point and find the true opt)

%% old gradient descent
err = 1.0;
lamda = 1e-4;
max_iter = 300;
iter = 1;
all_error = [];
all_y = [];

init_x = 0;
step = 0.1;% when step = 1, init_x = 0, weird thing happens, f(x_new) - f(x) = 0 while none of them are opt
while (err > lamda && iter <= max_iter)
    [y, gradient] = func(init_x);
    all_y = [all_y, y];
    % fixed step
    new_x = init_x - gradient*step;
    [new_y, new_gradient] = func(new_x);
    new_err = abs(new_y - y)^2; % 容忍度是两次y的变化量
    all_error = [all_error, new_err]; % 将容忍度存起来
    
    err = new_err; % 更新容忍度
    init_x = new_x; % 更新x
end

plot(all_y,'bo-');
hold on;


%% backtracing line search 回溯线搜索
% 在搜索方向上，先设置一个初始步长，如果步长太大，则缩减步长，直到合适为止
% still have to change the alpha beta for multiple times to see if x_final
% is the right one!
err = 1.0;
lamda = 1e-4;
max_iter = 300;
iter = 1;
all_error = [];
all_y = [];
% para for backtracing
% when alpha = 0.5, beta = 0.5, init_x = 0, only one iter will take to opt
alpha = 0.2; % 0 - 0.5
beta = 0.7;

init_x = 0;

while (err > lamda && iter <= max_iter)
    [y, gradient] = func(init_x);
    all_y = [all_y, y];
    % step is caled in the inner loop
    step = 1;
    % from convex optimization chap 9
    while (func(init_x - step*gradient) > func(init_x) - alpha*step*(gradient')*gradient)
        f_x = func(init_x - step*gradient)
        f_alpha_x = func(init_x) - alpha*step*(gradient')*gradient
        
        step = step * beta;
    end
    iter
    step
    new_x = init_x - gradient*step;
    [new_y, new_gradient] = func(new_x);
    new_err = abs(new_y - y)^2;
    all_error = [all_error, new_err];
  
    err = new_err;
    init_x = new_x;
end

plot(all_y,'ro-');

%% exact line search


%% newton method
% http://www.eecs.berkeley.edu/~wainwrig/ee227a/hw6soln_fa09.pdf
% convex optimization p484
clear f;

MAXITS = 500; % Maximum number of iterations
BETA = 0.5; % Armijo parameter
SIGMA = 0.1; % Armijo parameter

GRADTOL = 1e-7; % Tolerance for gradient

% load xinit.ascii;
% load A.ascii;
% load b.ascii

%x = xinit; % x is a vector
x = [0.1;0.2;0.3];
% m = size(A,1);
% n = size(A,2);

for iter=1:MAXITS,
    val = x'*log(x); % current function value
    f(iter) = val;
    grad = 1 + log(x); % current gradient
    hess = diag(1./x); % hessian of the function
    % !!! in newton method, hess >= 0, otherwise, use the gradient step, other
    % than the newton step
    assert
    msg = sprintf('iter: %d\n func val: %f\t grad: %f %f %f\n',iter, val, grad(1), grad(2), grad(3));
    disp(msg);
    % newt = -hess^(-1)*grad; descmag = grad'*(-newt) Page:487 in convex
%     temp = -[hess A'; A zeros(m,m)] \ [grad; zeros(m,1)];
%     newt = temp(1:n);
%     primal_lambda = temp(n+1:(n+m));
%     descmag = grad'*newt; % Check magnitude of descent
    %newt = -hess \grad;
    newt = - inv(hess)*grad;
    descmag = grad'*newt; % Check magnitude of descent

    %stopping criterion
    if (abs(descmag) < GRADTOL) 
        break; 
    end;
    
    % find step by Armijo
    t = 1;
    %while (min(x + t*newt) <= 0) t = BETA*t; end;
    % f(x+d_x) > f(x) + sigma*t*grad'*newt
    while ( ((x+t*newt)'*log(x+t*newt)) - val >= SIGMA*t*descmag)
        t = BETA*t;
    end;
    
    x = x + t*newt;
end;
% gradviol = norm(A*x -b,2);
% pstar = val;

%% with constraints
% http://www.cs.berkeley.edu/~pabbeel/cs287-fa11/slides/NonlinearOptimizationForOptimalControl.pdf

% constraints with Ax=b
% method 1: f(x) Ax = b => min f(z) = min f(Fz + x^) where F and x^ are the
% allfine feasible solutions, However, we need to find F and x^

% method 2: extended Newton Method
% prereq: 1. init point must be on Ax = b 
%         2. direction newt must be A*newt = 0;

