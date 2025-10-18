module PartonShower

export generateEvents
export showerEvent
export pTmin
export aSover

export WriteToLHE
export ShowerLHE

include("Constants.jl")
include("SplittingFunctions.jl")
include("LHEWriter.jl")
include("Shower.jl")

using LHEF
using EzXML
using ProgressBars


# Function to read a lhe file using LHEF and then shower then events
function ShowerLHE(inputFile::String)
    #Vector to hold the events in
    events::Vector{Event} = []

    #String to hold the event information
    eventinfo::String = ""

    # Get the information about the type of events this file contains
    open(EzXML.StreamReader, inputFile) do reader 

        #Loop to iterate the file until the event info is found
        while true
        #Check if the reader has the events info and then break
        if reader.name == "init"
            eventinfo = reader.content
            break
        end
            iterate(reader)

        end
    end

    
    #Remove the 'please cite' part of the init text
    eventinfo = split(eventinfo, "please cite")[1]
    
    # Read the lhe file and turn the LHEF events into events that the parton shower can work with
    lheevents = parse_lhe(inputFile)

    # List to hold each event's info
    eventsHeader = []

    for ev in lheevents
        newEvent = Event([], [])
        for p in ev.particles
            newP = Particle(p.id, p.status, 0, 1, (p.m), 0, p.px, p.py, p.pz, p.e, 0, [0, 0, 0, 0], p.color1, p.color2, 1, 0, 0, true, "", [])
            push!(eventsHeader, [ev.header.scale, ev.header.pid, ev.header.weight, ev.header.aqed, ev.header.aqcd])
            push!(newEvent.Jets, newP)
        end
        push!(events, newEvent)
    end

    # Shower the events
    showeredEvents::Vector{Event} = []
    for (i, ev) in tqdm(enumerate(events))
        newEvent = showerEvent(ev, pTmin, aSover)
        push!(showeredEvents, newEvent)
    end

    return showeredEvents, eventinfo, eventsHeader
end


# Function to take in a list of showered events and then write them to a lhe file
function WriteToLHE(showeredEvents::Vector{Event}, outputFile::String, eventInfo::String, eventHeaders::Vector{})
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
    writeLHE(outputFile, showeredEV, eventInfo, eventHeaders)

end

end

