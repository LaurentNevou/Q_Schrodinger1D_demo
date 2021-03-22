# Schroedinger1D_demo
This code solves the time independent Schroedinger equation in 1D with a constant mass.

![image](https://user-images.githubusercontent.com/35040499/111983702-b48ceb00-8b0a-11eb-8b66-944631713179.png)


It uses 4 different algorithms that can be switched ON/OFF:

-> The FDM: Finite Difference Method

-> The Scanning/Shooting method using the Euler approach

-> The PWE: Plane Wave Expansion method that solves the equation in the Fourier space

-> The TMM: Tranfer Matrix Method using the Scanning/Shooting method

Three different potentials are proposed:

-> The multilayer where the height and thickness of each layer can be adjusted

-> The parabolic

-> The sinus

Any potential can be loaded with homogeneous grid z [m] and V [eV] and be solved with the 3 first methods. For the TM method, a sequence of layer must be provided.

=> If you like it, don't forget the star!
