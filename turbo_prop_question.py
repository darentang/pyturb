from design_point import Engine, FlightCon, EngineInputs, TurboJet
import numpy as np
from pprint import pprint

n_shaft = 33000
gearing_ratio = 10
eta_gear = 0.95
D_prop = 1.6

# Inputs to the problem
inputs = EngineInputs()
inputs.pi_c = 7
inputs.tit = 1200
inputs.tot = 800
inputs.eta_c = 0.85
inputs.eta_t = 0.85
inputs.eta_m = 0.9
inputs.eta_cc = 0.99
inputs.lambda_cc = 0.05
inputs.m_dot = 2.75
inputs.eta_type = 'polytropic'

# Unknowns
inputs.symbolise('fuel_air_ratio')
inputs.symbolise('power_offtake')

# Flight conditions
flightcon = FlightCon(1828.8, 0.3)

# Set up cycle calcs and solve
engine = TurboJet(inputs, flightcon)
engine.solve()
engine.thrust()

pprint(engine.inputs.__dict__)
pprint(flightcon.__dict__)
pprint(engine.stations.T)
pprint(engine.stations.P)
pprint(engine.__dict__)

n_prop = n_shaft/gearing_ratio/60
omega_prop = n_prop*2*np.pi
Q = eta_gear*engine.inputs.power_offtake/omega_prop
J = flightcon.V/n_prop/D_prop
print(Q)
Cq = Q/(flightcon.rho * n_prop**2 * D_prop**5)
Cp = 2*np.pi*Cq
Cs = (flightcon.rho * flightcon.V**5/flightcon.P/n_prop**2)**(1/5)
print('Cq:', Cq)
print('Cp:', Cp)
print('Cs:', Cs)
print('J:', J)