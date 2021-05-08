# Sensitivity Study

This is a follow-along tutorial for a basic sensitivity study.  You can find the template for this tutorial in the templates directory at the root level of this repository ([here's a link](https://github.com/BYU-Aeronautics-Club/RCDesignSuite.jl/blob/master/templates/sensitivity.jl)).

!!! warning "Be Advised"
    This is **not** the only way to do a sensitivity study; it is simple one of many possible options.


## What is a sensitivity study?
Put simply, a sensitivity study is a study where we explore the sensitivities of the outputs to the inputs of a system.


## How do I start a sensitivity study?
A sensitivity study can start at either the input or output end.  We'll start at both and work our way toward the middle.
At the output side of things, we'll pick some performance criteria.  In this case, let's go with the Mission Scoring Criteria from the 2020-2021 AIAA DBF Competition Rules:

**Flight Mission 1:**
Failure = 0, Success = 1

**Flight Mission 2:**
$M2 = 1 + [N_{(\#containers/time)}/Max_{(\#containers/time)}]$
where $Max_{(\#containers/time)}$ is the highest $\#containers/time$ of all competing teams, and time is the time to complete 3 laps.

**Flight Mission 3:**
$M3 = 2 + [N_{(\#laps \times sensor length \times sensor weight)} / Max_{(\#laps \times sensor length \times sensor weight)}]$
where $Max_{(\#laps \times sensor length \times sensor weight)}$ is the highest $\#laps \times sensor length \times sensor weight$ score of all competing teams.

**Ground Mission**
$GM = [Min_{time} / N_{time}]$
where $Min_{time}$ is the fastest time of all competing teams.

**Design Report**
Range from 0 to 100 based on rubric.


For the sake of our example here, we'll chose a few design variables (inputs) to which the design will  likely be sensitive.  These are:

- Weight, $W$,
- Wing Area, $S$,
- Zero Lift Drag Coefficient, $C_{D_0}$,
- Available Battery Energy, $eb$,
- Total Propulsive Efficiency, $\eta$, and
- Available Thrust at Propeller, $T$.

Noting some of the specific payload items in the Flight Mission scores, we'll also probably need to include
- the number of containers, $ncon$, and
- sensor length times sensor weight, $slxsw$.

Take note that we are *not* simply choosing the direct inputs to the scoring functions.  We want to calculate things like time and number of laps based on our design to see how the scores are sensitive to the *design* rather than the scoring function inputs.


## Do I have to include all the inputs and outputs in my study?
For the purposes of this example sensitivity study, we'll assume success on the first flight mission, and we might as well assume a perfect score on the design report, but either way, these metrics are not sensitive to our design, so we exclude them.  We will also exclude the Ground Mission in this example since it is difficult to quantify *a priori* what the score might be, as physical testing is most likely needed.

That leaves us with Flight Missions 2 and 3 as our sensitivity metrics. That is to say, we are going to look at how changing our basic design variables affects the values for M2 and M3.

As for inputs, our choices are already very basic.  Each of these is comprised of more detailed variables. For example, the weight is a sum of all the component weights, the wing area is a function of span and chord, and the propulsive efficiency is a combination of battery, wire, ESC, motor, and propeller efficiencies.  For the sake of the sensitivity study, if a change in one of the detailed components yields the same change in output as another, then you should probably choose a more general input.  This is precisely why we chose the quantity sensor length times sensor weight rather than each of those separately.

!!! tip "Tip: How to know what to include"
    Another way to think about which inputs to include is to look at how they affect the derivative of the function. If they have the same effect, then include them together.

## Connecting the inputs and the outputs

This is where we arrive at the hard part of the sensitivity study.  In order to connect the inputs and the outputs, we need to create some sort of model.  To do so, we need to decide how simple or complex we want to be and what assumptions we are willing to make.  For the sake of this example, we are going to make some big assumptions and keep things quite simple.

### The Objective Function

We're going to start with a function we'll call ```objective```. Why? Mostly because how we have chosen to set things up in this example resembles an optimization problem, and that's a common name for this kind of function.  This is the function that takes in the basic inputs and outputs the overall "objective," which is the sum of the mission scores in our case.

```@
function objective(design_variables)

    # Do some stuff

    return M2 + M3

end
```

But what do we put in the middle?  How do we get from inputs to outputs?  There are a lot of options, but let's try the following:

```@example
import RCDesignSuite; rcds = RCDesignSuite;

function objective(design_variables; return_all=false)

    ### STEP 1: UNPACK DESIGN VARIABLES ###
    # Here is where you unpack the relevant design variables.
    W               = design_variables[1] # Aircraft Weight
    S               = design_variables[2] # Wing Area
    CD0             = design_variables[3] # Zero Lift Drag Coeff of Aircraft
    eta             = design_variables[4] # Total Propulsive Efficiency
    eb              = design_variables[5] # Available Battery Energy
    T               = design_variables[6] # Available Thrust at Prop.
    ncon            = design_variables[7] # Number of Containers
    slxsw           = design_variables[8] # Length x weight of sensor
    # Normalization Factors are the last 2 inputs in this example
    M2_norm_factor  = design_variables[end-1] # Flight Mission 2 Normalization Factor
    M3_norm_factor  = design_variables[end] # Flight Mission 3 Normalization Factor



    ### STEP 2a: CALCULATE FLIGHT MISSION 2 SCORE ###
    # Here is where you set up the inputs for flight mission 2
    M2_inputs = design_variables[[1:3;6:7]]

    # Call the mission2() function here (to be defined below)
    M2 = mission2(M2_inputs, M2_norm_factor)


    ### STEP 2b: CALCUALATE FLIGHT MISSION 3 SCORE ###
    # Here is where you set up the inputs for flight mission 3
    M3_inputs = design_variables[[4;5;8]]

    # Call the mission3() function here (to be defined below)
    M3 = mission3(M3_inputs, M2_norm_factor)



    ### STEP 3: SUM UP OBJECTIVES ###

    # Sum them up
    totalpoints = M2 + M3

    # Return the outputs (sometimes we might want to see them individually)
    if return_all
        return totalpoints, M2, M3
    else
        return totalpoints
    end

end
```

You'll notice that we still need to define some auxiliary functions, ```mission2()``` and ```mission3()```.  We could have just put everything in one function, but it's already a little messy, so we'll break it up into bite-sized chunks.

### The Mission Functions

For Flight Mission 2, remember our objective is $1 + [N_{(\#containers/time)}/Max_{(\#containers/time)}]$.  Number of containers is a direct input, and for the sake of this example, we're going to select a lap distance that seems reasonable (with the right models, you could calculate the lap distance).  We're also going to use one of the functions found in RCDesignSuite.jl, namely, ```maxvelocity()``` to find the velocity and then calculate the time to complete the mission.

```@example
function mission2(M2_inputs, M2_norm_factor, lap_distance=3000)

    W       = M2_inputs[1] # Aircraft Weight
    S       = M2_inputs[2] # Wing Area
    CD0     = M2_inputs[3] # Zero Lift Drag Coeff of Aircraft
    T       = M2_inputs[4] # Available Thrust
    ncon    = M2_inputs[5] # Number of Containers

    # Get Velocity
    V = rcds.maxvelocity(T, W, S, CD0)

    # Calculate Time
    time = 3*lap_distance/V

    # Return M2 Score
    return 1 + (ncon/time) / M2_norm_factor

end
```

For Flight Mission 3, the metric is $2 + [N_{(\#laps \times sensor length \times sensor weight)} / Max_{(\#laps \times sensor length \times sensor weight)}]$.  So we'll calculate the total number of laps using the ```endurance_time()``` function. Note that we are just making up some of the inputs for the endurance function for simplicity of this example.

!!! note "Note"
    Several of the inputs to the ```endurance_time()``` function could have been included as full design variables, but we're takeing some liberties here to keep things simple.


```@example
function mission3(M3_inputs, M2_norm_factor, lap_distance=3000)

    eta     = M3_inputs[1] # Total Propulsive Efficiency
    eb      = M3_inputs[2] # Available Battery Energy
    slxsw   = M3_inputs[3] # Length x weight of sensor

    # make up L and D to get an L/D of 10, which is reasonable.
    L = 10
    D = 1

    g = 9.81 #fact

    # make up a reasonable velocity
    V = 30 # fast enough to fly, but not win.

    # make up a take off and battery mass that makes sense.
    mb = 1
    mto = 4

    # find endurance time
    time = rcds.endurance_time(eb, eta, L, mb, g, V, D, mto)

    # calculate number of laps
    distance = V*time
    nlaps = distance/lap_distance

    # Return M3 score
    return 2 + (nlaps*slxsw) / M3_norm_factor

end
```

Now we should have all the functions we need to have a working objective function (remembering that this is NOT a real sensitivity study, and should NOT be used for any sort of competition or design, but rather as a template for how you might put together a real one).

Now we can move on to the actual sensitivity study part of things.



## Setting up the Sensitivity Study

This is the easy, but tedious, part.  All we need to do is to loop through each of the design variables and sweep them through a range of values, keeping all the rest constant, and then plot how the objective changes with the changes in the design variables.

But you're probably wondering where the norm_factor variables come into play, right? (Maybe not, it depends on how much you've been paying attention.)  For this example, we're actually going to run through the whole loop twice.  The first time, we're going to save the maximum values of the individual scores and use those as the normalization factors.  This isn't a great option, but we had to pick something.

Before we get to that, however, let's set up the framework for our study:

```@example
function run_sensitivity()

    # Make up some numbers, but in real  life, actually find some reasonable ones.
    W               = 5.0 # Aircraft Weight
    S               = 0.5 # Wing Area
    CD0             = 0.05 # Zero Lift Drag Coeff of Aircraft
    eta             = 0.65 # Total Propulsive Efficiency
    eb              = 10000.0 # Available Battery Power
    T               = 5.0 # Available Thrust at Prop.
    ncon            = 5.0 # Number of Containers
    slxsw           = 5.0 # Length x weight of sensor

    # Define input array for objective function.
    design_variables        = zeros(10)
    design_variables[1]     = W     # Aircraft Weight
    design_variables[2]     = S     # Wing Area
    design_variables[3]     = CD0   # Zero Lift Drag Coeff of Aircraft
    design_variables[4]     = eta   # Total Propulsive Efficiency
    design_variables[5]     = eb    # Available Battery Energy
    design_variables[6]     = T     # Available Thrust at Prop.
    design_variables[7]     = ncon  # Number of Containers
    design_variables[8]     = slxsw # Length x weight of sensor
    # Normalization Factors are the last inputs in this example
    design_variables[end-1] = 1.0 # Initialize M2_norm_factor to 1.0
    design_variables[end]   = 1.0 # Initialize M3_norm_factor to 1.0


    # Number of steps in the range
    N = 50 #number of points in ranges

    # Define range of percentages to do objective study over.
    r = range(-50,stop=50,length=N)*1e-2


    # Get objective Values
    obj, obj0, dobj = sensitivity(design_variables,r,N)


    ## Plot Sensitivities Separately
    intermediate_plots(design_varibles, obj, obj0, dobj;
    labels = [# TODO: Create a vector for the labels that match your design variables, so you know what the following plots are for.
            "W";
            "S";
            "CD0";
            "eta";
            "P";])


    ## Create Final Plots for Reports
    final_plots(obj, obj0, dobj, r, labels,;
            # TODO: identify the index of the variables you want
            toplot = [1;2;3;4;5],
            scalefactor = 1e3, # If needed, add a scale factor that helps the plot be more clear
            colors = ["C0";"C1";"C2";"C0";"C1";"C2"], # Pick your colors
            styles = ["-";"-";"-";"--";"--";"--"],  # Pick your styles.
            fs = (4.5,4.5*3/4), #figure size
            save_path = "../figs/", #save path
            )
end
```

Now for the details of running the sensitivity function.

```@example
function sensitivity(design_variables,r,N)


    # Initialize the objective function outputs array
    obj = zeros(length(design_variables),N)
    # Initialize the objective function derivatives array
    dobj = zeros(length(design_variables),N)


    ### GET NORMALIZATION FACTORS ###
    # For each input (other than normalization factors), loop through defined range (r) from the nominal

    # Initialize normalization factors.
    GM_norm_factor = 0.0
    M2_norm_factor = 0.0
    M3_norm_factor = 0.0

    # Run once to get maxima (or mean, or whatever you want to normalize by) for mission normalization factors.
    for i = 1:length(design_variables)-3 #(Don't include norm factors)

        # Create a copy of the original design variables so you can isolate each one in turn.
        design_variablescopy = copy(design_variables)

        for j = 1:N

            # Loop through each of the design variables across the defined relative range.
            design_variablescopy[i] = (1.0 + r[j])*design_variables[i]

            #try catch in case you get a combination that doesn't make sense or has a negative square root or something.
            try

                #throw away the total, we don't need it, we're just getting the maximum components for normalization now.
                _, GM_temp, M2_temp, M3_temp = objective(design_variablescopy)

                # Find Maximum normalization values.
                if GM_temp > GM_norm_factor
                    GM_norm_factor = GM_temp
                end
                if M2_temp > M2_norm_factor
                    M2_norm_factor = M2_temp
                end
                if M3_temp > M3_norm_factor
                    M3_norm_factor = M3_temp
                end

            catch

                #printl the cases that don't work so you can debug and make adjustments as needed.
                println("i: $i, j: $j")

            end #try/catch
        end #for range
    end #for number of design variables


    ## Put your Normalization Factors in your Design Variables Array:
    # Normalization Factors are the last 3 inputs in this example
    design_variables[end-2] = GM_norm_factor
    design_variables[end-1] = M2_norm_factor
    design_variables[end] = M3_norm_factor



    ### GET SENSITIVITIES ###

    # Run a second time to get actual values for the objectives, normalized by the max GM, M2, and M3 values just found.
    for i = 1:length(design_variables)-3 #(Don't include norm factors)

        # Same copying as before
        design_variablescopy = copy(design_variables)

        for j = 1:N

            # Same loop over range as before.
            design_variablescopy[i] = (1.0 + r[j])*design_variables[i]

            try

                # Don't return everything this time.
                obj[i,j] = objective(design_variablescopy)

            catch

                # Set failures to NaN so they don't show up in the plot.
                obj[i,j] = NaN

            end #try/catch
        end #for range

        # Spline the outputs of each design variable range.
        so = FLOWMath.Akima(r, obj[i,:])

        # Get the derivatives relative to the percent of input.
        for j = 1:N
            dobj[i,j] = FLOWMath.derivative(so,r[j])
        end

    end # for number of variables


    ## Calculate the objective value for the nominal case for getting the percent differences later.
    obj0 = objective(design_variables)


    ## Return the
    return obj, obj0, dobj

end
```

And finally, we probably want to first plot everything separately, thus the intermediate plotting function.  This will tell us if we selected the variables wisely and if any end up being identically sensitive.

```@example
function intermediate_plots(nvar, obj, obj0, dobj;
                        labels = [
                                "W";
                                "S";
                                "CD0";
                                "eta";
                                "P";])

    # Plot each variable on it's own plot so you can clearly see its sensitivity.
    for i=1:nvar-3
        figure()
        clf()
        plot(r,(obj[i,:].-obj0)./obj0,label=labels[i])
        legend()
    end

end
```

And for the final plot, we'll want everything on the same plot so we can fit it into a design proposal or report.

```@example
function final_plots(obj, obj0, dobj, r, labels,;
            # TODO: identify the index of the variables you want
            toplot = [1;2;3;4;5],
            scalefactor = 1e3, # If needed, add a scale factor that helps the plot be more clear
            colors = ["C0";"C1";"C2";"C0";"C1";"C2"], # Pick your colors
            styles = ["-";"-";"-";"--";"--";"--"],  # Pick your styles.
            fs = (4.5,4.5*3/4), #figure size
            save_path = "../figs/", #save path
            )


    # Get the plot labels you want to show for those variables
    toplotlabels = labels[toplot]


    #### PLOT FINAL FIGURE FOR REPORTS ####
    # Initialize the figure so that it will be about the right size for the paper (4x3 ratio), but also such that the legend doesn't overlap things weird.
    figure(figsize=fs)

    xlabel("Percent Change in Design Variable (from nominal)")
    ylabel("Normalized Change in Score (x $scalefactor)")

    # Plot the interesting things.
    for i=1:length(toplot)
        plot(r*1e2,scalefactor*(obj[toplot[i],:].-obj0)./obj0,label=toplotlabels[i],linestyle=styles[i],color=colors[i])
    end

    legend()

    #save the figure.
    savefig(save_path*"sensitivityobj.png",bbox_inches="tight")



    #### MAYBE PLOT THE DERIVATIVES INSTEAD, IF THAT MAKES MORE SENSE ####
    #initialize the figure so that it will be about the right size for the paper (4x3 ratio), but also such that the legend doesn't overlap things weird.
    figure(figsize=(4.5,4.5*3/4))

    xlabel("Percent Change in Design Variable (from nominal)")
    ylabel("Differential Change in Score (dScore/dx)")

    #plot the interesting things.
    for i=1:length(toplot)
        plot(r*1e2,dobj[toplot[i],:],label=toplotlabels[i],linestyle=styles[i],color=colors[i])
    end

    legend()

    #save the figure.
    savefig(save_path*"sensitivitydobj.png",bbox_inches="tight")

end
```

And if all goes well, you should have a nice plot of sensitivities after running the run_sensitivity function.
