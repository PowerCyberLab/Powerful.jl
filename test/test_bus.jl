using Test
using PowerFlowData
using Powerful.Models

case = PowerFlowData.parse_network(joinpath(@__DIR__, "..", "cases", "ieee14.raw"))
buses = case.buses

@testset "BusInput" begin
    bus_input = BusInput(case.buses)
    @test length(bus_input.i) == 14

    @test get_online_count(bus_input) == 14
end
