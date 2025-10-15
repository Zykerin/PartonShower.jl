
function writeLHE(infile::String, events::Vector{}, eventinfo::String, eventHeaders::Vector{})
    

    open(infile, "w") do file
    
    # Write the header for the lhe file
    write(file, "<LesHouchesEvents version =\"1.0\">\n")
    write(file, "<!--\n")
    write(file, "File generated with lhe Julia writer\n")
    write(file, "-->\n")
    write(file, "<init>")
    write(file, "\t" * eventinfo)
    write(file, "</init>\n")

    # Iterate through the events
    for (eventno, event) in enumerate(events)
        status = Int64[]
        momenta = []
        flavours = Int64[]
        colours = []
        anticolours = []
        helicities = []
        relations = []
        masses = []
        # Iterate through the particles in each event
        for p in event
            push!(momenta, [p[3], p[4], p[5], p[6]])
            append!(status, Int(p[2]))
            push!(masses, p[7])
            # Check if the particle is a final state particle
            if p[2] != 1
                push!(relations, [0,0])
            elseif p[2] == 1
                push!(relations, [1, 2])
            end
            append!(flavours, Int(p[1]))
            append!(helicities, 1)
            # Check if the particle is neither a quark nor a gluon and add no color to it
            if (abs(p[1]) < 0 || abs(p[1]) > 6) && abs(p[1]) !=21
                append!(colours, 0)
                append!(anticolours, 0)
            else
                append!(colours, Int(p[8]))
                append!(anticolours, Int(p[9]))
            end
            
        end
            # Print this event
            # Write this event's header info
            write(file, "<event>\n")
            write(file, string(length(momenta)) * "\t "* string(eventHeaders[eventno][2])* "\t " * string(eventHeaders[eventno][3]) * "\t" * string(eventHeaders[eventno][1])* "\t " * string(eventHeaders[eventno][4]) * " \t " * string(eventHeaders[eventno][5]) *"\n")
            # Write the actual particles info
            for i in range(1, length(momenta))
                p = momenta[i]
                particlestring = string(flavours[i]) * "\t" * string(status[i]) * "\t" * string(relations[i][1]) * "\t" * string(relations[i][2]) * "\t" * string(colours[i]) * "\t" * string(anticolours[i]) * "\t" * string(p[1]) * "\t" * string(p[2]) * "\t" * string(p[3]) *"\t" * string(p[4])* "\t" * string(masses[i]) * "\t0\t"* string(helicities[i]) * "\n" 
                write(file, particlestring)
                
            end
            write(file, "</event> \n")
        



    end
    end


end
