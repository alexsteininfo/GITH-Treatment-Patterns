#####################################################
### Functions to analyse simulations of tree data ###
#####################################################

# author: alexander stein
# comment: functions to evaluate the tree data as given by
#           Jessie Renton's SomaticEvolution package
#           using the data format SimpleTreeCell


# Copmute the branch widths from the tree
# Input:
# Output: 
function get_branchwidths(root)
    # Initilization
    branchwidths = Int32[]
    # Iterate over all nodes in the tree
    for cellnode in AbstractTrees.PreOrderDFS(root)
        width = treebreadth(cellnode)
        push!(branchwidths, width)
    end
    return branchwidths
end

# Compute the site frequencies from the tree
# Input:
# Output: 
function get_sitefreqs(root)
    # Initilization
    sitefreqs = Int32[]
    # Iterate over all nodes in the tree
    for cellnode in AbstractTrees.PreOrderDFS(root)
        width = treebreadth(cellnode)
        # add all mutations specific for this node
        for i in 1:cellnode.data.mutations
            push!(sitefreqs, width)
        end
    end
    return sitefreqs
end

# Compute the site frequencies from the tree
# Input:
# Output: 
function get_sitefreqs_selection(root)
    # Initilization
    sitefreqs_clone1 = Int32[]
    sitefreqs_clone2 = Int32[]
    numnodes = 0
    # Iterate over all nodes in the tree
    for cellnode in AbstractTrees.PreOrderDFS(root)
        width = treebreadth(cellnode)
        numnodes += 1
        # add all mutations specific for this node
        for i in 1:cellnode.data.mutations
            if(cellnode.data.clonetype == 1)
                push!(sitefreqs_clone1, width)
            else
                push!(sitefreqs_clone2, width)
            end
        end
    end
    #println("Number of nodes = ", numnodes)
    return sitefreqs_clone1, sitefreqs_clone2
end

# Compute the site frequencies from the tree
# Input:
# Output: 
function get_sitefreqs(root, t_detec)
    # Initilization
    sitefreqs_pre = Int32[]
    sitefreqs_new = Int32[]
    # Iterate over all nodes in the tree
    for cellnode in AbstractTrees.PreOrderDFS(root)
        width = treebreadth(cellnode)
        if(cellnode.data.birthtime <= t_detec)
            # add all mutations specific for this node
            for i in 1:cellnode.data.mutations
                push!(sitefreqs_pre, width)
            end
        else
            # add all mutations specific for this node
            for i in 1:cellnode.data.mutations
                push!(sitefreqs_new, width)
            end
        end
    end
    return sitefreqs_pre, sitefreqs_new
end

# Compute the division distribution from the tree
# Input:
# Output: 
function get_division_dist(root)
    # Initilization
    division_dist = Int32[]
    # Iterate over all leaves (aka living cells)
    for cellnode in Leaves(root)
        # Each cell has at least undergone one division (unless the tree is trivial)
        divisions = 1
        # This leave has a parent
        parent = AbstractTrees.parent(cellnode)
        # If the parent is not root, we add a division and go to the next parent
        while parent != nothing
            divisions += 1
            parent = AbstractTrees.parent(parent)
        end
        push!(division_dist, divisions)
    end
    return division_dist
end

# Compute the scMB from the tree
# Input:
# Output: 
function get_scMBs(root)
    # Initilization
    scMBs = Int32[]
    # Iterate over all leaves (aka living cells)
    for cellnode in Leaves(root)
        # Each cell has at least undergone one division (unless the tree is trivial)
        mutations = cellnode.data.mutations
        # This leave has a parent
        parent = AbstractTrees.parent(cellnode)
        # If the parent is not root, we add a division and go to the next parent
        while parent != nothing
            mutations += parent.data.mutations
            parent = AbstractTrees.parent(parent)
        end
        push!(scMBs, mutations)
    end
    return scMBs
end

# Create a histogram of integer-valued measurements
# Input: a list of integer-valued measurements (e.g. [1,2,1,2,1,0,0,5,2,4,13,15,...])
# Output: xaxis and yaxis of the histogram going from mininum to maximum
function create_histo(measurements)
    # Define the xaxis
    min = minimum(measurements)
    max = maximum(measurements)
    xaxis = Int32[i for i in min:max]
    # Initilization
    yaxis = Int32[0 for i in min:max]
    # Iterate over measurements and count occurances
    for i in measurements
        yaxis[i] += 1
    end
    return xaxis, yaxis
end

# Create histogram of measurements with predefined xaxis
# Input: a list of integer-valued measurements (e.g. [1,2,1,2,1,0,0,5,2,4,13,15,...])
#           and the relevant xaxis
# Output: yaxis of the histogram going from mininum to maximum
function create_histo(measurements, xaxis)
    # Initilization
    yaxis = Int32[0 for i in xaxis]
    # Iterate over measurements and count occurances
    min = xaxis[1]
    for i in measurements
        yaxis[i-min+1] += 1     # minimum does not need to start at 1
    end
    return yaxis
end

# Coarse grain scatter plot into Histogram
# Input: xaxis and yaxis of a scatter plot
# Output: xaxis and yaxis of a histogram with defined binwidth
function ScatterToHisto(xaxis, yaxis, binwidth)
    # Initilization
    xaxis_histo = Float64[]
    yaxis_histo = Float64[]

    new_x = 0
    new_y = 0

    # Compute new values and add to the list
    for (i,(x,y)) in enumerate(zip(xaxis,yaxis))
        if(mod(i,binwidth)!=0)
            new_x += x
            new_y += y
        else
            new_x += x
            new_y += y
            # Add to the new lists
            push!(xaxis_histo, new_x/binwidth)
            push!(yaxis_histo, new_y/binwidth)
            # Reset
            new_x = 0
            new_y = 0
        end
    end
    return xaxis_histo, yaxis_histo
end


################################################
### Analysing the exponent of simulated data ###
################################################

# Estimate alpha from raw data
function alpha_estimated(xvector, xmin)
    xvector_truncated = deleteat!(xvector, findall(x->x<xmin,xvector))
    n = length(xvector_truncated)
    sum = 0
    for xi in xvector_truncated
        sum += log(xi/(xmin-0.5))
    end
    return 1+n/sum
end

# Estimate alpha from sorted SFS
function get_exponent(SFS, kmin)
    n = 0
    sum = 0
    for (k,Sk) in enumerate(SFS)
        if(k>=kmin && Sk > 0)
            sum += Sk*log(k/(kmin-0.5))
            n += Sk
        end
    end
    return 1+n/sum
end

