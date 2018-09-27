# Schroedinger1D_demo
This code solves the time independent Schroedinger equation in 1D with a constant mass.

It uses 4 different algorithms that can be switched ON/OFF:

-> The FEM: Finite Elements Method

-> The Scaning method using the Euler approch

-> The PWE: Plane Wave Expansion method that sloves the equation in the Fourier space

-> The TMM: Tranfer Matrix Method

Three different potentials are proposed:

-> The multilayer where the height and thickness of each layer can be adjusted

-> The parabolic

-> The sinus

Any potential can be loaded with homogeneous grid z [m] and V [eV] and be solved with the 3 first methods. For the TM method, a sequence of layer must be provided.
