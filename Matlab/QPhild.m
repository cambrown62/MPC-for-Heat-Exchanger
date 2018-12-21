function eta=QPhild(H,f,A_cons,b);
clear all; clc
%E=H
%F=f
%M=A_cons
%gamma=b
%eta=x

H = [2 -1;-1 1];
f = [-1;0];
A_cons = [-2 0;0 -1;3 3];
b = [0;0;-4];

% H = [2 2;3 1];
% f = [-2;-2];
% A_cons = [1 0;-2 0;0 1;0 -1];
% b = [1;0;3;0];

[n1,m1] = size(A_cons);
eta = -H\f;
kk=0;

for i=1:n1
    if (A_cons(i,:)*eta>b(i));
        kk=kk+1;
    else
        kk=kk+0;
    end
end

if kk==0
    return
end

P = A_cons*(H\A_cons');
%P
d = (A_cons*(H\f)+b);
%d
[n,m] = size(d);
x_ini = zeros(n,m);
lambda = x_ini;
al = 10;

for km=1:400
%find the elements in the solution vector one by one
%km could be larger if the Langrange multiplier has a slow
%convergence rate
    lambda_p=lambda;
    for i=1:n
        %i
        w = P(i,:)*lambda-P(i,i)*lambda(i,1);
        w = w+d(i,1);
        %w
        la = -w/P(i,i);
        %la
        lambda(i,1) = max(0,la);
    end
    %lambda
    %lambda_p
    al = (lambda-lambda_p)'*(lambda-lambda_p);
    %al
    if (al<10e-8)
        break
    end
end
%km

eta = -H\f - H\A_cons'*lambda;