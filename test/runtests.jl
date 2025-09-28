using PartonShower
using Test

@testset "Momentum Conservation" begin
   myEvents, sigma, error = generateEvents(1E4, 1E6)
   showeredEvents = []
   for ev in (myEvents)
        newEvent = showerEvent(ev, pTmin, aSover)
        push!(showeredEvents, newEvent)
        totmom = [0, 0, 0, 0]
        for p in newEvent
            if p.status ==1
                totmom += [p.px, p.py, p.pz, p.E]
            end
        end
        @test totmom ≈ [0, 0, 0, 0] atol=0.0001
   end

end
