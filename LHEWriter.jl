


function writeLHE(infile, events, shat, ECM, sigma, stddev)

    open(infile, "w") do file
    write(file, "<LesHouchesEvents version =\"1.0\">\n")
    write(file, "<!--\n")
    write(file, "File generated with lhe Julia writer\n")
    write(file, "-->\n")
    write(file, "<init>\n")
    write(file, "\t-11\t 11\t" * string(ECM/2) * "\t" * string(ECM/2) * "\t 0 \t 0 \t 7\t 7 \t 1 \t 1\n")
    write(file, "\t" * string(sigma) * "\t" * string(stddev) * "\t1.00000 \t9999\n")
    write(file, "</init>\n")

    for (eventno, event) in enumerate(events)
        ng = 0
        status = Int64[]
        momenta = []
        flavours = Int64[]
        colours = []
        anticolours = []
        helicities = []
        relations = []
        for p in event
            push!(momenta, [p[3], p[4], p[5], p[6]])
            append!(status, Int(p[2]))
            if p[2] == -1
                push!(relations, [0,0])
            elseif p[2] == 1
                push!(relations, [1, 2])
            end
            append!(flavours, Int(p[1]))
            append!(helicities, 1)
            if abs(p[1]) == 11
                append!(colours, 0)
                append!(anticolours, 0)
            else
                append!(colours, Int(p[8]))
                append!(anticolours, Int(p[9]))
            end
            
            
            #=
            if abs(p[1]) > 0 && abs(p[1]) < 6 # q or qbar -> this only works for the specific process e+e- -> qqbar
                
                if p[1] < 0
                    append!(colours, 0)
                    append!(anticolours, 501)
                elseif p[1] > 0
                    append!(colours, 501)
                    append!(anticolours, 0)
                end

            end
            if p[1] == 21 # Gluons are singlets for now
                ng = ng + 1
                append!(colours, 500 + 2 * ng)
                append!(anticolours, 500+ 2 * ng)
            end
            =#
        end
        
        
       

            write(file, "<event>\n")
            write(file, string(length(momenta)) * "\t 9999\t 1.000000\t" * string(sqrt(shat))* "\t 0.0078125 \t 0.1187\n")
            for i in range(1, length(momenta))
                p = momenta[i]
                mass = 0
                particlestring = string(flavours[i]) * "\t" * string(status[i]) * "\t" * string(relations[i][1]) * "\t" * string(relations[i][2]) * "\t" * string(colours[i]) * "\t" * string(anticolours[i]) * "\t" * string(p[1]) * "\t" * string(p[2]) * "\t" * string(p[3]) *"\t" * string(p[4])* "\t" * string(mass) * "\t0\t"* string(helicities[i]) * "\n" 
                write(file, particlestring)
                
            end
            write(file, "</event> \n")
        



    end
    end


end
