# Big Wave simulation

Example of Big Wave movement simulation in Cartesian coordinates by Gerstner waves distribution.

**Algorithm:**

In this model, the following distribution in the time step was used (no symmetrization of the solution was carried out, we are only interested in the front of Big Wave movement along OZ):

<pre>
<a href="https://github.com/NuclearRazor/gerstner_waves_distribution/blob/master/common/img/distribution.png"><img src="https://github.com/NuclearRazor/gerstner_waves_distribution/blob/master/common/img/distribution.png" align="middle">
</a>
</pre>

𝜔1 − 𝑤𝑎𝑣𝑒_𝑠𝑖𝑧𝑒
𝜔2 − 𝑤𝑖𝑛𝑑_𝑎𝑙𝑖𝑔𝑚𝑒𝑛𝑡
𝜔3 − 𝑤𝑎𝑣𝑒_𝑠𝑝𝑒𝑒𝑑
𝐴 ∈ [𝑚𝑖𝑛𝐴, 𝑚𝑎𝑥𝐴] − wave amplitude (gaussian magnitude)

The last series component specifies the distribution of small waves, subject to the boundary conditions for the Herterstern wave, this function is periodic.

The minimum algorithm steps is:

1. Set initial parameters (wind speed, grid size, initial wave amplitude, minimum of wave size, shift wave value).
2. Scale coordinates by grid step.
3. Phases projections calculation of the coordinates of the Gerstner wave by the Fourier method.
4. Gradient projections processing of the coordinates of the Gerstner wave.
5. Calculate fz at the current step (this solution infinitely oscillates).
6. Update current step, increase (or decrease) shift, repeat from 3-nd state.

**Input**:
- wind speed (-w)
- surface dim (-d)
- initial amplitude (-a)
- wave size (-l)

**Run**:
Run main script with next parameters for example:

main.py -w 1.5 -d 100 -a 0.60 -l 0.06

**Examples**:

<pre>
<a href="https://github.com/NuclearRazor/gerstner_waves_distribution/blob/master/common/img/view_1.png"><img src="https://github.com/NuclearRazor/gerstner_waves_distribution/blob/master/common/img/view_1.png" align="middle">
</a>
</pre>

wind speed = 1.5, surface dim = 100, amplitude = 0.5, wave size = 0.006

<pre>
<a href="https://github.com/NuclearRazor/gerstner_waves_distribution/blob/master/common/img/view_2.png"><img src="https://github.com/NuclearRazor/gerstner_waves_distribution/blob/master/common/img/view_2.png" align="middle">
</a>
</pre>

wind speed = 1.5, surface dim = 100, amplitude = 0.5, wave size = 0.06

<pre>
<a href="https://github.com/NuclearRazor/gerstner_waves_distribution/blob/master/common/img/view_3.png"><img src="https://github.com/NuclearRazor/gerstner_waves_distribution/blob/master/common/img/view_3.png" align="middle">
</a>
</pre>

wind speed = 1.5, surface dim = 100, amplitude = 0.5, wave size = 0.003
