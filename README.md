# Compositional Model

This project implements a compositional reservoir model in Python to simulate the phase behavior of multicomponent hydrocarbon systems. The model explicitly tracks fluid composition and phase equilibrium, providing a more detailed and rigorous description of reservoir fluids than black-oil formulations. It is intended for educational purposes and as a foundation for future compositional reservoir simulation development.

## Requirements
- Python 3.10+ (It is not recommended to use the latest version)
- Dependencies listed in `requirements.txt` (numpy, math, barril)

## Installation and Usage
Clone the repository and create a virtual environment:

```bash
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```
Then run the main script:

```bash
python main.py
```
## Project Structure
- `main.py` – main execution script
- `requirements.txt` – project dependencies
- `venv/` – virtual environment

## Thermodynamic Model
Based on classical thermodynamic principles, the model solves the phase-equilibrium problem through a compositional flash calculation. The flash formulation determines phase compositions and phase fractions by enforcing equality of component fugacities between phases, using a cubic equation of state with the Soave–Redlich–Kwong correlation. Phase stability and equilibrium are assessed through Gibbs free energy criteria inherent to the fugacity-based formulation. For the implicits nonlinear equations we solved using the Newton–Raphson method. Once equilibrium compositions and phase fractions are obtained, volumetric and thermophysical properties are calculated from the equation-of-state results.


