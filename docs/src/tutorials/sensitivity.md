# Sensitivity Study

This is a follow-along tutorial for a basic sensitivity study.  You can find the template for this tutorial in the templates directory at the root level of this repository ([here's a link](../../../templates/sensitivity.jl)).


## What is a sensitivity study?
Put simply, a sensitivity study is a study where we explore the sensitivities of the outputs to the inputs of a system.


## How do I start a sensitivity study?
A sensitivity study can start at either the input or output end.  We'll start at both and work our way toward the middle.
At the output side of things, we'll pick some performance criteria.  In this case, let's go with the Mission Scoring Criteria from the 2020-2021 AIAA DBF Competition Rules:

**Flight Mission 1:**
Failure = 0, Success = 1

**Flight Mission 2:**
$M2 = 1 + [N_{(\#containers/time)}/Max_{(\#containers/time)}]$
where $Max_{(\#containers/time)}$ is the highest $\#containers/time$ of all competing teams.

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
- Available Battery Power, $P$, and
- Total Propulsive Efficiency, $\eta$.

Noting some of the specific payload items in the Flight Mission scores, we'll also probably need to include
- the number of containers, $ncon$, and
- sensor length times sensor weight, $slxsw$.

Take note that we are *not* simply choosing the direct inputs to the scoring functions.  We want to calculate things like time and number of laps based on our design to see how the scores are sensitive to the *design* rather than the scoring function inputs.


## Do I have to include all the inputs and outputs in my study?
For the purposes of this example sensitivity study, we'll assume success on the first flight mission, and we might as well assume a perfect score on the design report, but either way, these metrics are not sensitive to our design, so we exclude them.  We will also exclude the Ground Mission in this example since it is difficult to quantify *a priori* what the score might be, as physical testing is most likely needed.

That leaves us with Flight Missions 2 and 3 as our sensitivity metrics. That is to say, we are going to look at how changing our basic design variables affects the values for M2 and M3.

As for inputs, our choices are already very basic.  Each of these is comprised of more detailed variables. For example, the weight is a sum of all the component weights, the wing area is a function of span and chord, and the propulsive efficiency is a combination of battery, wire, ESC, motor, and propeller efficiencies.  For the sake of the sensitivity study, if a change in one of the detailed components yields the same change in output as another, then you should probably choose a more general input.  This is precisely why we chose the quantity sensor length times sensor weight rather than each of those separately.

!!! tip

    Another way to think about which inputs to include is to look at how they affect the derivative of the function. If they have the same effect, then include them together.

## Connecting the inputs and the outputs

This is where we arrive at the hard part of the sensitivity study.  In order to connect the inputs and the outputs, we need to create some sort of model.  To do so, we need to decide how simple or complex we want to be and what assumptions we are willing to make.  For the sake of this example, we are going to make some big assumptions and keep things quite simple.

