immutable PetriNet
    g::GraphViz.Graph
end

function petrinet(model::Network; rankdir=:LR)
    reactions = model.reactions
    b = IOBuffer()
    modelid = string(model.id)
    modelid = replace(modelid, r"(\W\s?|_)", "")
    modelid = replace(modelid, " ", "_")

    # graph attributes
    write(b,
    """
    digraph $(modelid) {
    // Graph Attributes
    layout=dot;
    rankdir=$(rankdir)
    """)
    # add reactions
    for r in values(reactions)
        add_reaction!(b, r)
    end
    write(b, "}")
    g = GraphViz.Graph(takebuf_string(b))
    return PetriNet(g)
end

function add_reaction!(b::IOBuffer, r)
    rname     = r.id
    reactants = r.reactants
    products  = r.products

    # add reaction node style
    write(b, "\"", r.id, "\"", " [shape=box, style=filled];\n")

    # add each reactant
    for x in keys(reactants)
        write(b, "\"", x, "\"", " -> ", "\"", r.id, "\"", ";\n")
    end

    # add each product
    for x in keys(products)
        write(b, "\"", r.id, "\"", " -> ", "\"", x, "\"", ";\n")
    end
    write(b, "\n")
    return b
end

# From Gadfly

function default_mime()
    "image/svg+xml"
end

import Base: display
import Base.REPL: REPLDisplay
import Base.Multimedia: @try_display, xdisplayable

function display(p::PetriNet)
    displays = Base.Multimedia.displays
    for i = length(displays):-1:1
        m = default_mime()
        if xdisplayable(displays[i], m, p)
             @try_display return display(displays[i], m, p)
        end

        if xdisplayable(displays[i], p)
            @try_display return display(displays[i], p)
        end
    end
    invoke(display,(Any,),p)
end

function display(d::REPLDisplay, ::MIME"image/svg+xml", p::PetriNet)
    filename = string(tempname(), ".svg")
    output = open(filename, "w")
    GraphViz.writemime(output, "image/svg+xml", p.g)
    close(output)
    open_file(filename)
end

function open_file(filename)
    if OS_NAME == :Darwin
        run(`open $(filename)`)
    elseif OS_NAME == :Linux || OS_NAME == :FreeBSD
        run(`xdg-open $(filename)`)
    elseif OS_NAME == :Windows
        run(`$(ENV["COMSPEC"]) /c start $(filename)`)
    else
        warn("Showing petri nets is not supported on OS $(string(OS_NAME))")
    end
end
