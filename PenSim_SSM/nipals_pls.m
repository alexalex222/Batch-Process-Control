function [w_star, c, t, u, beta, p, E, F, w] = nipals_pls(x, y, nPC)
rows = size(x,1);
if(rows ~= size(y,1))
    error('x and y must have the same number of observations');
end

col_x = size(x,2);
col_y = size(y,2);

% allocate variables
w = zeros(col_x, nPC);
t = zeros(rows, nPC);
u = zeros(rows, nPC);
c = zeros(col_y, nPC);

E = zeros(rows, col_x, nPC);
F = zeros(rows, col_y, nPC);

x_deflated = x;
y_deflated = y;


for a = 1:nPC
    % select a random initial guess for ua
    u(:,a) = rand(rows, 1);
    
    ua_old = zeros(rows, 1);
    
    i = 0;
    while (i<800 && (i ==0 || sqrt((ua_old - u(:,a))'*(ua_old - u(:,a)))>eps))
        % regress x onto u (remember terminoligy regress y onto x)
        w(:,a) = (1/(u(:,a)'*u(:,a))*u(:,a)'*x_deflated)';
        
        %normalize the weights
        w(:,a) = 1/(sqrt(w(:,a)'*w(:,a)))*w(:,a);
        
        % regress x onto wa' to get scores ta
        t(:,a) = 1/(w(:,a)'*w(:,a))*x_deflated*w(:,a);
        
        %regress columns in y onto ta to get loadings ca
        c(:,a) = (1/(t(:,a)'*t(:,a))*t(:,a)'*y_deflated)';
        
        % store old ua for use in convergnece testing
        ua_old = u(:,a);
        
        % regrss rows of Y onto ca to calcualte scores ua
        u(:,a) = 1/(c(:,a)'*c(:,a))*y_deflated*c(:,a);
        
        % iteration counter safty 
        i = i+1;
    end
    
    % calculate loadings matrix from converged scores by regressing x onto
    % ta
    p(:,a) = 1/(t(:,a)'*t(:,a))*x_deflated'*t(:,a);
    
    % deflate x
    E(:,:,a) = x_deflated - t(:,a)*p(:,a)';
    x_deflated = E(:,:,a);
    
    % deflate y
    F(:,:,a) = y_deflated - t(:,a)*c(:,a)';
    y_deflated = F(:,:,a);
end

w_star = w*(p'*w)^-1;

beta = w_star*c';
end