using Econometrics, Test, Random
function main()
    # Utilities
    # sortbyc
    @testset "sortbyc" begin
        @test sortbyc()
    end    
    # eye
    @testset "eye" begin
        @test eye()
    end    
    # dstats
    @testset "dstats" begin
        @test dstats()
    end    
    # vech
    @testset "vech" begin
        @test vech()
    end
    # LinearRegression
    # OLS
    x = [ones(10,1) (1:10)]
    y = [10.,9.,8.,7.,6.,5.,4.,3.,2.,1.]
    b = ols(y,x,silent=true)[1]
    @testset "ols" begin
        @test b[1,1] / 11.0 ≈ 1.0 atol=1e-5
    end    
    # ML
    @testset "mle" begin
        b = mle()
        @test b[1,1] / 0.49466773000112035 ≈ 1.0 atol = 1e-5
    end    
    # GMM
    @testset "gmm" begin
        b = gmm()
        @test b[1,1] / 0.4946717833359474 ≈ 1.0  atol = 1e-5
    end    
    # Optimization
    # samin
    xopt, fopt = samin(1)
    @testset "samin" begin
        @test xopt[1] ≈ 0.0 atol=1e-5
        @test fopt ≈ 2.0 atol=1e-5
    end
    # fmincon
    @testset "fmincon" begin
        @test fmincon()
    end    
    # fminunc
    @testset "fminunc" begin
        @test fminunc()
    end    
    #=
    # this is ill-advised, need to test with non-random data,
    # so it doesn't break when julia changes
    # npreg
    @testset "npreg" begin
        yhat = npreg(1);
        @test yhat[1,1] ≈ 0.02080984791302587
    end
    =#
end
main()

