module PartonShower

export generateEvents
export writeLHE
export showerEvent
export pTmin
export aSover


include("Constants.jl")
include("SplittingFunctions.jl")
include("LHEWriter.jl")
include("Structures.jl")
include("Shower.jl")
include("HardProccessEventGenerator.jl")
using LHEF
using ProgressBars



# Function to shower each event in a given list
function ShowerEvents(events::Vector{Event})
    showeredEvents = []
    print("Showering " * string(length(events))* " events \n")
    for (i, ev) in tqdm(enumerate(events))
        newEvent = showerEvent(ev, pTmin, aSover)
        push!(showeredEvents, newEvent)
    end
    return showeredEvents
end

# Function to take in a list of showered events and then write them to a lhe file
function WriteToLHE(showeredEvents::Vector{Event}, outputFile::String)
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


end