# Notes on optimizing the diffusivity integrator

### June 26, 2020

I'm using ShowParticlesInSEAK.m to do a profile.

Apparently the accessor methods with polshape objects are really slow. So I need to save this down differently on initialization.
