function cue()
    # example of CUE GMM: draws from N(0,1)
    y = randn(100,1)
    # 3 moment conditions
    moments = theta -> [y-theta[1] (y.^2.0).-theta[2] (y.-theta[1]).^3.0]
    weight = theta -> inv(NeweyWest(moments(theta)))
    obj = theta -> vec(mean(moments(theta),1))'*weight(theta)*vec(mean(moments(theta),1))
    theta = [0.0, 1.0]
    results = fminunc(obj, theta)
end    


