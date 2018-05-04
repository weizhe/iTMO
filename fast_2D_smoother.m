% reference: Fast Global Image Smoothing Based on Weighted Least Squres, Dongbo Min
function [ u ] = fast_2D_smoother( f, g, lambda, T )
% separable global smoother for 2D image smoothing
% linear smoothing with the exponential kernal(x,y directions) multiple times 
% approximates the Gaussian kernal
% T: iteration num
% g: guide image, g = f for image filtering, g ~= f for joint filtering

if nargin ~=4
    error('Wrong parameters');
end

f = double(f);

% lambda = 30, sigma_c = 0.03(for image rage 0~1)

u = f;
[H,W] = size(f);

% if image filtering, padded_f, otherwise padded_g
padded_g = padarray(g, [1 1], 'replicate', 'post');

vec_h = padded_g(1:(end-1), 2:end) - g;
vec_v = padded_g(2:end,1:(end-1)) - g;
alpha = 1.2;
wx_mat = 1./(abs(vec_h).^alpha+0.0001);
wx_mat(:,end) = 0;
wy_mat = 1./(abs(vec_v).^alpha+0.0001);
wy_mat(end,:) = 0;
% wx_mat = exp(-abs(vec_h)/0.01);
% wy_mat = exp(-abs(vec_v)/0.01);

for t=1:1:T
%     lambda_t = 1.5 * power(4,T-t) / (power(4,T)-1) * lambda *2;
    lambda_t = 1.5 * power(4,T-t) / (power(4,T)-1);
%     lambda_t = T / t;
    
    for y=1:1:H
        fh = u(y,:)';
        wx1 = [0, wx_mat(y, 1:(end-1))]';
        wx2 = [wx_mat(y, 1:(end-1)), 0]';
        
        
        a = -wx1*lambda*lambda_t;
        b = 1 + (wx1+wx2)*lambda*lambda_t;
        c = -wx2*lambda*lambda_t;
        
        uh = one_d_smoother([a b c], fh);    
        u(y,:) = uh';
    end
  
    
    for x=1:1:W
        fv = u(:,x);
        wy1 = [0;wy_mat(1:(end-1), x)];
        wy2 = [wy_mat(1:(end-1), x);0];   
        
        a = -wy1*lambda*lambda_t;
        b = 1 + (wy1+wy2)*lambda*lambda_t;
        c = -wy2*lambda*lambda_t;        
        uv = one_d_smoother([a b c], fv);
        
        u(:,x) = uv;
    end    
end


end


function [ u ] = one_d_smoother( A, f )
% 1D Fast Global Smoother
% solve equation A * u = f, A is a tridiagnal matrix
% u and f are vectors of the same size

% A should be square or [a b c]
if size(A,1)~=size(A,2)
    if size(A,2)~=3
        error('A should be square or [a b c]');
    else
        a = A(:,1);
        b = A(:,2);
        c = A(:,3);
    end
else
    a = [0; diag(A, -1)];
    b = diag(A);
    c = [diag(A, 1); 0];
end

if size(f,2)~=1
    error('Input f should be a vector.');
end


n = size(f,1); % f -> n x 1 vector
u = zeros(n,1);


% 
%         % normalize first two rows
%         % first two rows of a are 0
%         ch(1:2) = ch(1:2) ./ bh(1:2);
%         fh(1:2) = fh(1:2) ./ bh(1:2);
%         bh(1:2) = 1;
%         for i=3:1:width
%             fh(i) = (fh(i) - ah(i)*fh(i-2)) / (bh(i) - ah(i)*ch(i-2));
%             ch(i) = ch(i) / (bh(i) - ah(i)*ch(i-2));
%         end
%         
%         uh = zeros(width, 1);
%         uh(width) = fh(width);
%         uh(width-1) = fh(width-1);
%         for i=(width-2):-1:1
%             uh(i) = fh(i) - ch(i)*fh(i+2);
%         end
% %         u(y,:) = (vec_h(y,:)' + uh)';
%         u(y,:) = uh';



% make b1 = 1
c(1) = c(1) / b(1);
f(1) = f(1) / b(1);
b(1) = 1;


% forward 
for i=2:1:(n-1)
    c(i) = c(i) / (b(i)-c(i-1)*a(i));
    f(i) = (f(i)-f(i-1)*a(i)) / (b(i)-c(i-1)*a(i));
end
f(n)=(f(n)-f(n-1)*a(n)) / (b(n)-c(n-1)*a(n));

% backward
u(n) = f(n);
for i=(n-1):-1:1
    u(i) = f(i) - c(i)*u(i+1);
end

end

