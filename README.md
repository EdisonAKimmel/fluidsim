# Fluid Simulation
Credit to Jos Stam for the techniques used.
This project implements a stable fluid dynamics solver using a semi-Lagrangian advection scheme and implicit diffusion. 
The simulation runs in real time, making it suitable for educational use or integration into game engines.
Left-click to inject fluid into the sim,
Right-click + drag to move around fluid

Getting started
-----------------------
1. Install Rust and Cargo from [rustup.rs](https://rustup.rs/)
2. Run this code:
```bash
git clone https://github.com/EdisonAKimmel/fluidsim.git
cd fluidsim
cargo run --release
