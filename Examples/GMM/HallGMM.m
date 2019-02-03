hall = importdata('hall.m');
n = size(hall, 1);
options = optimset('TolX',1e-15,'TolFun',1e-12);

% growth rate of consumption, and 2 lags
c = hall(3:n,1);
c1 = hall(2:n-1,1);
c2 = hall(1:n-2,1);

% equally weighted returns, and 2 lags
r = hall(3:n,2);
r1 = hall(2:n-1,2);
r2 = hall(1:n-2,2);

%{
% value weighted returns, and 2 lags
r = hall(3:n,3);
r1 = hall(2:n-1,3);
r2 = hall(1:n-2,3);
%}

% sample size after computing lags
n = size(c,1); 

% start values for parameters: discount and utility
thetastart = [0.5; 0.5];

% Euler equation error
e = @(theta) theta(1,:)*r.*c.^(theta(2,:) - 1) -1;

% instruments
inst = [ones(n,1) c1 r1];
%inst = [ones(n,1) c1 r1 c2 r2];

% initial weight matrix
W = eye(size(inst,2));

% moments
m = @(theta) (1/n)*inst'*e(theta); % average moments for obj. fun.
ms = @(theta) diag(e(theta)) * inst; % moment contribs, to compute weight matrix

% initial consistent
[thetahat, obj] = fminunc(@(theta) n*m(theta)'*W*m(theta), thetastart);  % note the n in there, obj function is scaled to converge to chi^2.
fprintf('initial consistent\n');
thetahat
obj

% compute efficient weight matrix
mm = ms(thetahat);
omega = mm'*mm/n;
W = inv(omega);
% efficient
[thetahat, obj] = fminunc(@(theta) n*m(theta)'*W*m(theta), thetahat);
fprintf('two step\n');
thetahat
obj

% an iteration
mm = ms(thetahat);
omega = mm'*mm/n;
W = inv(omega);
% efficient
[thetahat, obj] = fminunc(@(theta) n*m(theta)'*W*m(theta), thetahat);
fprintf('iterated\n');
thetahat
obj

% what about Newey-West, even though theory suggests it's not needed?
mm = ms(thetahat);
lags = floor(n^0.25); % choose lag length
omega = NeweyWest(mm, lags);
W = inv(omega); 
[thetahat, obj] = fminunc(@(theta) n*m(theta)'*W*m(theta), thetahat);
fprintf('using NW\n');
thetahat
obj

% iterate NW
mm = ms(thetahat);
lags = floor(n^0.25); % choose lag length
omega = NeweyWest(mm, lags);
W = inv(omega); 
[thetahat, obj] = fminunc(@(theta) n*m(theta)'*W*m(theta), thetahat);
fprintf('NW iterated\n');
thetahat
obj
