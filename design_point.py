import numpy as np
from sympy.solvers.solveset import nonlinsolve, linsolve
from sympy import Symbol
from collections import defaultdict
from functools import partial

def named_symbol(name, key):
    return Symbol(name + key)

class KeyDict(defaultdict):
    def __init__(self, f):
        super().__init__(None)
        self.f = f

    def __missing__(self, key):
        self[key] = self.f(key)
        return self[key]

class FlightCon:
    alt = 1000
    P = 10
    T = -40+273.15
    M = 0.88

class EngineStations:
    T = KeyDict(partial(named_symbol, 'T'))
    P = KeyDict(partial(named_symbol, 'P'))

class EngineInputs:
    """
    Default inputs
    """
    n_spools = 1
    pi_f = 0
    pi_c = 0
    eta_f = 0.9
    eta_c = 0.9
    eta_t = 0.9
    eta_m = 0.99
    lambda_cc = 1
    eta_cc = 0.992
    ram_recovery = 0.99
    bypass_ratio = 1
    lhv = 43.1e6
    tit = 1625
    tot = None
    fuel_air_ratio = 0.02
    gamma_c = 1.4
    gamma_h = 1.33
    cph = 1156.7
    cpc = 1004.5
    power_offtake = 0
    m_dot = 1
    eta_type = 'isentropic'

class TurboJet:
    def __init__(self, inputs, flightcon):
        self.inputs = inputs
        self.flightcon = flightcon
        self.stations = EngineStations()

    def compressor(self, station_in, station_out, pi):
        """
        Calculate compressor outlet temp and pressure based on inlet conditions and pressure ratio
        """
        self.stations.P[station_out] = self.stations.P[station_in]*self.inputs.pi_c
        if self.inputs.eta_type == 'polytropic':
            self.stations.T[station_out] = self.stations.T[station_in]*(self.stations.P[station_out]/self.stations.P[station_in])**(1/self.exp_c/self.inputs.eta_c)
        else:
            self.stations.T[station_out] = self.stations.T[station_in]*(1 + (self.inputs.pi_c**(1/self.exp_c)-1)/self.inputs.eta_c)

    def turbine(self, station_in, station_out):
        """
        Calculate turbine outlet pressure based on temp difference
        """
        
        if self.inputs.eta_type == 'polytropic':
            self.stations.P[station_in] = self.stations.P[station_in]*(self.stations.T[station_out]/self.stations.T[station_out])**(self.exp_h/self.inputs.eta_t)
        else:
            self.stations.T[station_out+'s'] = self.stations.T[station_in] - (self.stations.T[station_in] - self.stations.T[station_out])/self.inputs.eta_t
            self.stations.P[station_out] = self.stations.P[station_in]*(self.stations.T[station_out+'s']/self.stations.T[station_out])**self.exp_h


    def solve(self):
        """
        Solve design point calcs
        """
        self.exp_c = (self.inputs.gamma_c / (self.inputs.gamma_c - 1))
        self.exp_h = (self.inputs.gamma_h / (self.inputs.gamma_h - 1))

        # Station 0
        self.stations.T['0'] = self.flightcon.T*(1+0.5*(self.inputs.gamma_c -1)*self.flightcon.M**2)
        self.stations.P['0'] = self.flightcon.P*(1+0.5*(self.inputs.gamma_c -1)*self.flightcon.M**2)**self.exp_c

        # Station 2
        self.stations.T['2'] = self.stations.T['0']
        self.stations.P['2'] = self.stations.P['0']*self.inputs.ram_recovery

        # Station 3
        self.compressor('2', '3', self.inputs.pi_c)

        # enthalpy balance
        cc_exit_enthalpy = (self.inputs.tit-288)*(1+self.inputs.fuel_air_ratio)*self.inputs.cph
        cc_entry_enthalpy = (self.stations.T['3']-288)*self.inputs.cpc
        cc_excess_enthalpy = self.inputs.eta_cc*self.inputs.lhv*self.inputs.fuel_air_ratio
        cc_enthalpy_balance = cc_entry_enthalpy + cc_excess_enthalpy - cc_exit_enthalpy

        # Compressor
        compressor_power = self.inputs.m_dot*self.inputs.cpc*(self.stations.T['3'] - self.stations.T['2'])
        # Turbine
        turbine_power = self.inputs.m_dot*self.inputs.eta_m*(1 + self.inputs.fuel_air_ratio)*self.inputs.cph*(self.inputs.tit-self.inputs.tot)
        compressor_turbine_balance = self.inputs.power_offtake + compressor_power - turbine_power

        # Solve the equations
        symbols = cc_enthalpy_balance.free_symbols.union(compressor_turbine_balance.free_symbols)
        sol = linsolve([cc_enthalpy_balance, compressor_turbine_balance], list(symbols))
        for sym, val in zip(symbols, list(sol)[0]):
            exec(f'self.{str(sym)}={val}')

        # Station 4
        self.stations.T['4'] = self.inputs.tit
        self.stations.P['4'] = self.stations.P['3']*(1-self.inputs.lambda_cc)
        
        # Station 5
        self.stations.T['5'] = self.inputs.tot


        self.turbine('4', '5')




if __name__ == '__main__':
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
    inputs.fuel_air_ratio = Symbol('inputs.fuel_air_ratio')
    inputs.power_offtake = Symbol('inputs.power_offtake')

    # Flight conditions
    flightcon = FlightCon()
    flightcon.T = 276.26
    flightcon.P = 81.2
    flightcon.M = 0.3
    
    # Set up cycle calcs and solve
    engine = TurboJet(inputs, flightcon)
    engine.solve()

    print(engine.stations.T)
    