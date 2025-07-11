#Binary Distillation Simulation: Ethanol–Water System with Reflux Optimization
**Overview**
This project presents a dynamic simulation of binary distillation for an ethanol–water mixture using MATLAB. The system is modeled to explore transient effects due to time-varying feed composition and optimize energy consumption through reflux ratio control. The McCabe–Thiele method is used to track theoretical stage requirements, and energy duties for the condenser and reboiler are computed over a 10-hour run.

The simulation captures realistic process behavior by incorporating:

Fitted vapor–liquid equilibrium (VLE) data

Analytical and numerical reflux ratio estimation

Real-time stage tracking using staircase construction

Material and energy balance calculations

#Key Features
Feed Conditions:

Feed rate: 100 mol/h (normalized to 1 kg/s in energy calculations)

Feed composition x_f(t): sinusoidally varied in the range 0.28–0.42 mole fraction

Distillate composition target: ≥ 95 mol% ethanol

Bottom composition: fixed at 5 mol% ethanol

Equilibrium Modeling:

Vapor–liquid data fitted using a rational curve model:
y = (a * x) / (1 + b * x + c * x^2)

Parameters fitted to known ethanol–water equilibrium data

Reflux Ratio Control:

Minimum reflux ratio (R_min) estimated analytically using the operating line slope at feed composition

Actual reflux ratio enforced: R = max(R_min, 1.20), with adaptive increase up to 2.5 if needed to meet stage constraints

Simulation kept R = 1.20 throughout, ensuring minimal energy use without violating purity constraints

Stage Tracking:

Theoretical stage count computed dynamically using a custom staircase algorithm based on the McCabe–Thiele graphical method

Minimum of 10 theoretical stages enforced at all time points

#Energy Calculations:
<img width="400" height="344" alt="image" src="https://github.com/user-attachments/assets/3057f491-2020-4456-b74d-4db7f5f34d5c" />

Reboiler and condenser duties computed from material balance:
Q_R = B * L, Q_C = D * L
where L = 8.55e5 J/kg (latent heat), and B, D are bottom and distillate flow rates respectively

Total energy duties over the 10-hour run:
Reboiler duty: 2.85 MJ
Condenser duty: 5.70 MJ

#Energy Optimization:
Reflux optimization led to an estimated 18 percent reduction in reboiler duty compared to a fixed reflux baseline
Energy savings achieved without increasing stage count or sacrificing product purity
