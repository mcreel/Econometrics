% These are the candidate auxiliary statistics for ABC estimation of
% the simple DSGE model of Creel and Kristensen (2013)

function Z = aux_stat(data)

    function bad = check_bad_data(data)
        bad = sum(any(isnan(data)))||sum(any(isinf(data)))||sum(any(std(data)==0));
    end 

    bad_data = false;
    % check for bad inputs
    if check_bad_data(data) || sum(any(data<=0))
        bad_data = true;
    else
        % variables coming in are y c n r w (output cons hours intrate wages)
        % bound data: final result is between 1e-6 and 1000
        test = data > 100;
        data = test*100+ (1-test).*data;
        logdata = log(data);
        lagdata = lags(data,1); % lag one time
        laglogdata = lags(logdata,1); % lag one time
        n = rows(data);
        % drop missing first obs after adding lags
        data = data(2:n,:);
        logdata = logdata(2:n,:);
        lagdata = lagdata(2:n,:);
        laglogdata = laglogdata(2:n,:);

        % check for bad inputs
        checkme = [data lagdata logdata laglogdata];
        if check_bad_data(checkme)
            bad_data = true;
        else
            output = data(:,1);     % output
            cons = data(:,2);       % consumption
            hours = data(:,3);      % hours
            intrate = data(:,4);    % intrate
            wages = data(:,5);      % wages

            logoutput = logdata(:,1);   % output
            logcons = logdata(:,2);     % consumption
            loghours = logdata(:,3);    % hours
            logintrate = logdata(:,4);  % intrate
            logwages = logdata(:,5);    % wages

            lagoutput = lagdata(:,1);   % output
            lagcons = lagdata(:,2);     % consumption
            laghours = lagdata(:,3);    % hours
            lagintrate = lagdata(:,4);  % intrate
            lagwages = lagdata(:,5);    % wages

            nobs = rows(data);

            % alpha, rho1, sig1 (production function)
            % use a MA in investment as proxy for capital
            investment = output - cons;
            capitalproxy = investment(1:nobs-11,:) +...
                investment(2:nobs-10,:)+...
                investment(3:nobs-9,:)+...
                investment(4:nobs-8,:)+...
                investment(5:nobs-7,:)+...
                investment(6:nobs-6,:)+...
                investment(7:nobs-5,:)+...
                investment(8:nobs-4,:)+...
                investment(9:nobs-3,:)+...
                investment(10:nobs-2,:);
            logcapitalproxy = log(capitalproxy);
            x = [logcapitalproxy loghours(12:nobs,:)];
            x = [ones(rows(x),1) x];
            y = logoutput(12:nobs,:);
            w = laglogdata(12:nobs,:);
            w = [ones(rows(w),1) w];
            xhat = w*ols(x,w); % generalized IV estimator 
            b = ols(y,xhat);
            e = y-x*b;
            junk = [e lag(e,1)];
            junk = junk(2:rows(junk),:);
            y = junk(:,1);
            x = junk(:,2);
            rho1 = corr(x,y);
            e = y-x*rho1;
            sig1 = e'*e/nobs;
            Z = [1-b(3,:); rho1; sig1];

            % gam, rho2, sig2 (1/MRS=wage)
            x = [ones(nobs,1) logcons];
            y = logwages;
            w = laglogdata;
            w = [ones(rows(w),1) w];
            xhat = w*ols(x,w);
            b = ols(y,xhat); % generalized IV estimator
            e = y-x*b;
            junk = [e lag(e,1)];
            junk = junk(2:rows(junk),:);
            y = junk(:,1);
            x = junk(:,2);
            rho2 = corr(y,x);
            e = y-x*rho2;
            sig2 = e'*e/nobs;
            Z = [Z; b; rho2; sig2];

            % standard devs. and correlations
            data = data(:,1:5);
            [data, m, s] = st_norm(data); % keep means and std. devs., the VAR uses standardized and normalized

            % AR(1)s
            maxlag = 1;
            data = [data lags(data, 1)]; % add lags
            data = data(2:end,:); % drop rows with missing
            y = data(:,1:5);
            x = data(:,6:end);
            n = rows(y);
            rhos = zeros(5,1);
            es = zeros(n,5);
            for i = 1:5
                [rho, junk, e] = ols(y(:,i),x(:,i));
                rhos(i,:) = rho;
                es(:,i) = e;
            end        
            varv = vech(cov(es)); % AR(1) error covariance elements 
            % ratios
            s1 = mean(cons./output);
            s2 = mean(intrate./wages);
            s3 = mean(cons./hours);
            Z = [Z; m(:); s(:); rhos(:); varv(:); s1; s2; s3];
            Z = real(Z);
            if check_bad_data(Z)
                bad_data = true;
            end
        end
    end
    if bad_data
        Z = -1000*ones(7 + 5 + 5  + 5 + 15 + 3,1);
    end
end


