using .SBI
using Test


@test_throws DomainError SBI.blinding_index(-50, 50, 50, 50)
@test_throws DomainError SBI.blinding_index(0, 0, 50, 50)
@test_throws DomainError SBI.blinding_index(50, 50, 0, 0)

@testset "Perfect blinding." begin
    @test round(SBI.blinding_index(50, 50, 50, 50)[1, "est"]; digits=4) == 0
    @test round(SBI.blinding_index([50 50; 50 50])[1, "est"]; digits=4) == 0
    @test round(SBI.blinding_index(50, 50, 50, 50)[1, "p_value"]; digits=4) == 1
end

@testset "Zero blinding." begin
    @test round(SBI.blinding_index(50, 50, 50, 50)[1, "est"]; digits=4) == 0
    @test round(SBI.blinding_index(50, 50, 50, 50)[1, "p_value"]; digits=4) == 1
end

@testset "Example from Newcombe 1998, PMID: 9595617, Table II, Method 10, Column (a)." begin
    @test round(SBI.blinding_index(56, 14, 48, 32)[1, "est"]; digits=4) == 0.2
    @test round(SBI.blinding_index([56 48; 14 32])[1, "est"]; digits=4) == 0.2
    @test round(SBI.blinding_index(56, 14, 48, 32)[1, "lwr_ci"]; digits=4) == 0.0524
    @test round(SBI.blinding_index(56, 14, 48, 32)[1, "upr_ci"]; digits=4) == 0.3339
end

@testset "Example from Newcombe 1998, PMID: 9595617, Table II, Method 10, Column (b)." begin
    @test round(SBI.blinding_index(9, 1, 3, 7)[1, "est"]; digits=4) == 0.6
    @test round(SBI.blinding_index(9, 1, 3, 7)[1, "lwr_ci"]; digits=4) == 0.1705
    @test round(SBI.blinding_index(9, 1, 3, 7)[1, "upr_ci"]; digits=4) == 0.8090
end

@testset "Example from Newcombe 1998, PMID: 9595617, Table II, Method 10, Column (c)." begin
    @test round(SBI.blinding_index(6, 1, 2, 5)[1, "est"]; digits=4) == 0.5714
    @test round(SBI.blinding_index(6, 1, 2, 5)[1, "lwr_ci"]; digits=4) == 0.0582
    @test round(SBI.blinding_index(6, 1, 2, 5)[1, "upr_ci"]; digits=4) == 0.8062
end

@testset "Example from Newcombe 1998, PMID: 9595617, Table II, Method 10, Column (d)." begin
    @test round(SBI.blinding_index(5, 51, 0, 29)[1, "est"]; digits=4) == 0.0893
    @test round(SBI.blinding_index(5, 51, 0, 29)[1, "lwr_ci"]; digits=4) == -0.0381
    @test round(SBI.blinding_index(5, 51, 0, 29)[1, "upr_ci"]; digits=4) == 0.1926
end

@testset "Example from Newcombe 1998, PMID: 9595617, Table II, Method 10, Column (e)." begin
    @test round(SBI.blinding_index(0, 10, 0, 20)[1, "est"]; digits=4) == 0
    @test round(SBI.blinding_index(0, 10, 0, 20)[1, "lwr_ci"]; digits=4) == -0.1611
    @test round(SBI.blinding_index(0, 10, 0, 20)[1, "upr_ci"]; digits=4) == 0.2775
end

@testset "Example from Newcombe 1998, PMID: 9595617, Table II, Method 10, Column (f)." begin
    @test round(SBI.blinding_index(0, 10, 0, 10)[1, "est"]; digits=4) == 0
    @test round(SBI.blinding_index(0, 10, 0, 10)[1, "lwr_ci"]; digits=4) == -0.2775
    @test round(SBI.blinding_index(0, 10, 0, 10)[1, "upr_ci"]; digits=4) == 0.2775
end

@testset "Example from Newcombe 1998, PMID: 9595617, Table II, Method 10, Column (g)." begin
    @test round(SBI.blinding_index(10, 0, 0, 20)[1, "est"]; digits=4) == 1
    @test round(SBI.blinding_index(10, 0, 0, 20)[1, "lwr_ci"]; digits=4) == 0.6791
    @test round(SBI.blinding_index(10, 0, 0, 20)[1, "upr_ci"]; digits=4) == 1
end

@testset "Example from Newcombe 1998, PMID: 9595617, Table II, Method 10, Column (h)." begin
    @test round(SBI.blinding_index(10, 0, 0, 10)[1, "est"]; digits=4) == 1
    @test round(SBI.blinding_index(10, 0, 0, 10)[1, "lwr_ci"]; digits=4) == 0.6075
    @test round(SBI.blinding_index(10, 0, 0, 10)[1, "upr_ci"]; digits=4) == 1
end