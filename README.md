# 2D-Rocket-Trajectory

The RK4 numerical method of integration is used to time-step through rocket trajectories around single-body and stationary two-body (Earth-Moon) gravity wells.
On the former, the user has the choice of launching from Earth, Mars or Jupiter.
Note that, on the latter, the introduction of another body (moon) leads the simulation to be particularly sensitive to the initial velocity parameter. 
The input prompt statement is designed to aid in choosing parameters that result in a useful simulation.

Using the principle of energy conservation, the rocket’s initial total energy should be equal to the rocket’s final total energy.
After completing the simulation, a percentage uncertainty is printed to the screen (calculated by considering the difference between initial and final total energies)

Tested on Python 3.9.12
