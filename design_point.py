import numpy as np

class FlightCon:
     alt = None
     P = None
     T = None
     M = None

class EngineStations:
    T = dict()
    P = dict()

class EngineInputs:
    n_spools = 1
    pi_f = 0
    pi_c = 0
    eta_f = 1
    eta_c = 1
    lambda_cc = 1
    eta_cc = 0.992
    ram_recovery = 1
    bypass_ratio = 1

    gamma_c = 0
    gamma_h = 0

    eta_type = 'isentropic'

class TurboJet:
    def __init__(self, inputs, flightcon):
        self.inputs = inputs
        self.flightcon = flightcon
        self.stations = EngineStations()

    def solve(self):
        exp = (self.inputs.gamma / (self.inputs.gamma - 1))

        # Station 0
        self.stations.T['0'] = self.flightcon.T*(1+0.5*(self.inputs.gamma_c -1)*self.flightcon.M**2)
        self.stations.P['0'] = self.flightcon.P*(1+0.5*(self.inputs.gamma_c -1)*self.flightcon.M**2)**exp

        # Station 2
        self.stations.T['2'] = self.stations.T['0']
        self.stations.P['2'] = self.stations.P['0']*self.inputs.ram_recovery

        # Station 3
        self.stations.P['3'] = self.stations.P['2']*self.inputs.pi_c
        if self.inputs.eta_type == 'polytropic':
            self.stations.T['3'] = self.stations.T['2']*(self.stations.P['3']/self.stations.P['2'])**(1/exp/self.inputs.eta_c)
        else:
            self.stations.T['3'] = self.stations.T['2']*(1 + (self.inputs.pi_c**(1/exp/self.inputs.eta_c)))

        # Station 4
        self.stations.P['4'] = self.stations.P['3']*(1-self.inputs.lambda_cc)
        