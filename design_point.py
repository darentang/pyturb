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
    lhv = 0

    fuel_air_ratio = 0

    gamma_c = 0
    gamma_h = 0

    cph = 0
    cpc = 0

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

        # Combustion chamber
        cc_exit_enthalpy = self.stations.T['4']*(1+self.inputs.fuel_air_ratio)*self.inputs.cph
        cc_entry_enthalpy = self.statins.T['3']*self.inputs.cpc
        cc_excess_enthalpy = self.inputs.eta_cc*self.inputs.lhv


        # Compressor
        compressor_power = self.inputs.cpc*(self.stations.T['3'] - self.stations.T['2'])

        # Turbine
        turbine_power = (1 + self.inputs.fuel_air_ratio)*self.inputs.cph*(self.stations.T['4']-self.stations.T['5'])

        if self.inputs.eta_type == 'polytropic':
            self.stations.P['5'] = 