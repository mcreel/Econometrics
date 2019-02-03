using Test, Random
function main()
    # OLS
    x = [ones(10,1) (1:10)]
    Random.seed!(1)
    y = rand(10,1)
    b = ols(y,x,silent=true)[1]
    @testset "ols" begin
        @test b[1,1] ≈ 0.07655377367377458
    end    
    # samin
    xopt, fopt = samin()
    @testset "samin" begin
    @test xopt[1] ≈ 0.0 atol=1e-5
    @test fopt ≈ 2.0 atol=1e-5
    end
    # npreg
    @testset "npreg" begin
        yhat = npreg();
        @test yhat[1,1] ≈ -0.8108308248103517
    end    
end
main()

