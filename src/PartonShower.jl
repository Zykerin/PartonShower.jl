module PartonShower


export generateEvents
export WriteToLHE
export showerEvent
export pTmin
export aSover
export ShowerLHE

include("Constants.jl")
include("SplittingFunctions.jl")
include("LHEWriter.jl")
include("Structures.jl")
include("Shower.jl")
include("HardProccessEventGenerator.jl")




using LHEF
using ProgressBars


# Function to read a lhe file using LHEF and then shower then events
function ShowerLHE(inputFile::String)
    events::Vector{Event} = []


    # Read the lhe file and turn the LHEF events into events that the parton shower can work with
    lheevents = parse_lhe(inputFile)

    #Get the center of mass energy 
    ECM = lheevents[1].header.scale

    for ev in lheevents
        newEvent = Event([], [])
        for p in ev.particles
            newP = Particle(p.id, p.status, 0, 1, (p.m), 0, p.px, p.py, p.pz, p.e, 0, [0, 0, 0, 0], p.color1, p.color2, 1, 0, 0, true, "", [])
            push!(newEvent.Jets, newP)
        end
        push!(events, newEvent)
    end

    # Shower the events
    showeredEvents::Vector{Event} = []

    for (i, ev) in tqdm(enumerate(events))
        print(ev)
        newEvent = showerEvent(ev, pTmin, aSover)
        push!(showeredEvents, newEvent)
    end

    return showeredEvents, ECM
end



# Function to take in a list of showered events and then write them to a lhe file
function WriteToLHE(showeredEvents::Vector{Event}, outputFile::String, ECM::Float64, sigma::Float64, error::Float64)
    showeredEV = []
    # Turn the shower format into one that is readable by the LHEWriter
    for ev in showeredEvents
        sParts = []
        for p in ev.AllParticles
            push!(sParts, [p.id, p.status, p.px, p.py, p.pz, p.E, p.m, p.color, p.antiColor])
        end
        push!(showeredEV, sParts)
    end


    print("Writing to " * outputFile * " \n")
    writeLHE(outputFile, showeredEV, ECM^2, ECM, sigma, error)

end

eventhz, energyhz = ShowerLHE("ee_hz_ggmumu.lhe.gz")



end
