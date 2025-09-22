include("Constants.jl")
include("SplittingFunctions.jl")
include("LHEWriter.jl")
include("Structures.jl")
include("Shower.jl")
include("HardProccessEventGenerator.jl")
using LHEF
using ProgressBars
using ArgParse

# Variable to decide whether or not to generate the events from the hard proccess as well or read from an already generated lhe file 
# Options are: "generate" or "lhe"
whichEvents::String = "lhe"
global numGen::Int64 = 1E6

arguments = ArgParseSettings()

# The table to implement command line arguments
@add_arg_table arguments begin
    "--generate"
        help = "Generate the hard proccess events then shower them"
        action = :store_true
    "--lhefile"
        help = "Option to shower events from a lhe file"
        action = :store_true
    "--infile"
        help = "File to get hard proccess events from"
        arg_type = String
        default = ""
        required = "--lhefile" in ARGS
    "-N"
        help = "Number of events to generate"
        arg_type = Int64
        required = "--Generate" in ARGS
    "outfile"
        help = "The file to output the showered events"
        arg_type = String
        required = true
      
end


args = parse_args(arguments)

outputFile::String = args["outfile"]

# Parsing the inputted arguments
if args["generate"]
    global whichEvents = "generate"
    global numGen = args["N"]
elseif args["lhefile"]
    global whichEvents = "lhe"
    global outputFile = args["infile"]
end

outputFile::String = args["outfile"]
if last(outputFile, 4) != ".lhe"
      outputFile *= ".lhe"
end

error::Float64 = 0.1
sigma::Float64 = 1.2
myEvents = []


# Case to read an lhe file
if whichEvents == "lhe"
    print("LHE file selected. \n")
    error = 0.1
    sigma = 1.2
    inputfile::String = args["infile"]
    
    print("Reading Events from " * inputFile * "\n")
    events = parse_lhe(inputFile)

    for ev in events
        newEvent = Event([], [])
        for p in ev.particles
            newP = Particle(p.id, p.status, 0, 1, (p.m), 0, p.px, p.py, p.pz, p.e, 0, [0, 0, 0, 0], p.color1, p.color2, 1, 0, 0, true, "", [])
            push!(newEvent.Jets, newP)
        end
        push!(myEvents, newEvent)

    end
# Case to generate the hard proccess events (e+e- -> qqbar)
elseif whichEvents == "generate"
    print("Generate events selected\n")
    print("Now generating " * string(numGen)* " events\n")
    myEvents, sigma, error = generateEvents(numGen, 1E6)
end
showeredEvents = []


print("Showering " * string(length(myEvents))* " events \n")
for (i, ev) in tqdm(enumerate(myEvents))
    newEvent = showerEvent(ev, pTmin, aSover)
    push!(showeredEvents, newEvent)
end


showeredEV = []
# Turn my format into one that is readable by the LHEWriter
for ev in showeredEvents
    sParts = []
    for p in ev.AllParticles
        push!(sParts, [p.id, p.status, p.px, p.py, p.pz, p.E, p.m, p.color, p.antiColor])
    end
    push!(showeredEV, sParts)
end


print("Writing to " * outputFile * " \n")
writeLHE(outputFile, showeredEV, ECM^2, ECM, sigma, error)
