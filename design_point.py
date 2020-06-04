import numpy as np
from sympy.solvers.solveset import nonlinsolve, linsolve
from sympy.solvers import solve
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
    """
    Takes in altitude and speed, turn it into static temp, pressure and density
    """
    def __init__(self, alt, M):
        self.alt, self.M = alt, M
        mat = np.loadtxt('isa.txt', delimiter='\t')
        elevation = mat[:, 0]
        mat = mat[:, 1:]
        self.T = np.interp(alt, elevation, mat[:, 0])
        self.P = np.interp(alt, elevation, mat[:, 1]) * 100
        a = np.interp(alt, elevation, mat[:, -1])
        self.V = M*a

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

    def symbolise(self, var_name):
        """
        Make an input a sympy symbol
        """
        exec(f'self.{var_name} = Symbol("inputs.{var_name}")')

class Engine:
    def __init__(self, inputs, flightcon):
        self.inputs = inputs
        self.flightcon = flightcon
        self.stations = EngineStations()

        # Specific heat exponents for hot and cold
        self.exp_c = (self.inputs.gamma_c / (self.inputs.gamma_c - 1))
        self.exp_h = (self.inputs.gamma_h / (self.inputs.gamma_h - 1))

    def compressor(self, station_in, station_out, pi):
        """
        Calculate compressor outlet temp and pressure based on inlet conditions and pressure ratio
        """
        self.stations.P[station_out] = self.stations.P[station_in]*self.inputs.pi_c
        if self.inputs.eta_type == 'polytropic':
            self.stations.T[station_out] = self.stations.T[station_in]*(self.stations.P[station_out]/self.stations.P[station_in])**(1/self.exp_c/self.inputs.eta_c)
        else:
            self.stations.T[station_out] = self.stations.T[station_in]*(1 + (self.inputs.pi_c**(1/self.exp_c)-1)/self.inputs.eta_c)

    def turbine(self, station_in, station_out, out_temp):
        """
        Calculate turbine outlet pressure based on temp difference
        """
        
        self.stations.T[station_out] = out_temp

        if self.inputs.eta_type == 'polytropic':
            self.stations.P[station_out] = self.stations.P[station_in]*(self.stations.T[station_out]/self.stations.T[station_out])**(self.exp_h/self.inputs.eta_t)
        else:
            self.stations.T[station_out+'s'] = self.stations.T[station_in] - (self.stations.T[station_in] - self.stations.T[station_out])/self.inputs.eta_t
            self.stations.P[station_out] = self.stations.P[station_in]*(self.stations.T[station_out+'s']/self.stations.T[station_out])**self.exp_h

    def inlet(self, station_in, station_out, ram_recovery):
        """
        Calculate the total pressure and temp at inlet, account for ram recovery
        """
        
        self.stations.T[station_in] = self.flightcon.T*(1+0.5*(self.inputs.gamma_c -1)*self.flightcon.M**2)
        self.stations.P[station_in] = self.flightcon.P*(1+0.5*(self.inputs.gamma_c -1)*self.flightcon.M**2)**self.exp_c
        
        self.stations.T[station_out] = self.stations.T[station_in]
        self.stations.P[station_out] = self.stations.P[station_in]*ram_recovery

    def nozzle(self, station_in, station_out, gamma=None, m_dot=None):
        """
        Adiabatic nozzle
        """
        if not gamma:
            gamma = self.inputs.gamma_h

        if not m_dot:
            m_dot = self.inputs.m_dot
        
        self.stations.T[station_out] = self.stations.T[station_in]
        npr = self.stations.P[station_in]/self.flightcon.P
        pi_crit = (1 + (gamma-1)/2)**(gamma/(gamma-1))

        if npr > pi_crit:
            # Choked
            self.stations.P[station_out+'_static'] = self.stations.P[station_in]/pi_crit
            self.stations.T[station_out+'_static'] = self.stations.T[station_out]/(1 + (gamma-1)/2)
            M = 1
        else:
            # Adapted
            self.stations.P[station_out+'_static'] = self.flightcon.P
            self.stations.T[station_out+'_static'] = self.stations.T[station_in]*(1/npr)**((gamma-1)/gamma)
            M = np.sqrt(2/(gamma-1)*(npr)**((gamma-1)/gamma)-1)
        
        rho = self.stations.P[station_out+'_static']/(287*self.stations.T[station_out+'_static'])
        V = M*np.sqrt(gamma*287*self.stations.T[station_out+'_static'])
        A = m_dot/(rho*V)

        return A, V

    def afterburner(self):
        raise NotImplementedError

    def combustion_chamber(self, station_in, station_out, tit, lambda_cc, far,
                           gamma_in=None, gamma_out=None):
        """
        Enforce combustion chamber enthalpy balance
        """
        if not gamma_in:
            gamma_in = self.inputs.gamma_c
            cp_in = gamma_in/(gamma_in-1)*287
        
        if not gamma_out:
            gamma_out = self.inputs.gamma_h
            cp_out = gamma_out/(gamma_out-1)*287

        self.stations.T[station_out] = tit
        self.stations.P[station_out] = self.stations.P[station_in]*(1-lambda_cc)

        cc_exit_enthalpy = (tit-288)*(1+far)*cp_out
        cc_entry_enthalpy = (self.stations.T[station_in]-288)*cp_in
        cc_excess_enthalpy = self.inputs.eta_cc*self.inputs.lhv*far
        cc_enthalpy_balance = cc_entry_enthalpy + cc_excess_enthalpy - cc_exit_enthalpy

        self.enforce(cc_enthalpy_balance)

    def compressor_turbine_balance(self, comp_station_in, comp_station_out, turb_temp_in, turb_temp_out, comp_m_dot=None, turb_m_dot=None):
        """
        Calculate power difference between turbine and compressor
        """

        if not comp_m_dot:
            comp_m_dot = self.inputs.m_dot

        if not turb_m_dot:
            turb_m_dot = self.inputs.m_dot

        # Compressor
        compressor_power = comp_m_dot*self.inputs.cpc*(self.stations.T[comp_station_out] - self.stations.T[comp_station_in])
        # Turbine
        turbine_power = turb_m_dot*self.inputs.eta_m*(1 + self.inputs.fuel_air_ratio)*self.inputs.cph*(turb_temp_in - turb_temp_out)
        compressor_turbine_balance = self.inputs.power_offtake + compressor_power - turbine_power
        
        self.enforce(compressor_turbine_balance)

    def enforce(self, equation):
        """
        Enforce equation to be 0, save the numerical solution
        """
        sols = solve(equation, dict=True)[0]
        for k, v in sols.items():
            exec(f'self.{k}={v}')

class TurboJet(Engine):
    def solve(self):
        # Station 0 and 2
        self.inlet('0', '2', self.inputs.ram_recovery)

        # Station 3
        self.compressor('2', '3', self.inputs.pi_c)
        
        # Enforce combustion chamber balance
        self.combustion_chamber('3', '4', self.inputs.tit, self.inputs.lambda_cc, self.inputs.fuel_air_ratio)
        
        # Enforce power balance
        self.compressor_turbine_balance('2', '3', self.inputs.tit, self.inputs.tot)

        # Station 4 and 5        
        self.turbine('4', '5', self.inputs.tot)

        self.A9, self.V9 = self.nozzle('5', '9')
class TurboFanSepExhaust(Engine):
    def solve(self):
        bypass_mdot = self.inputs.m_dot*self.inputs.bypass_ratio

        # Station 0 and 2
        self.inlet('0', '2', self.inputs.ram_recovery)

        # Fan
        self.inlet('2', '25', self.inputs.pi_f)

        # HPC
        self.compressor('25', '3', self.inputs.pi_c)
        
        # Enforce combustion chamber balance
        self.combustion_chamber('3', '4', self.inputs.tit, self.inputs.lambda_cc, self.inputs.fuel_air_ratio)
        
        # Enforce HP Power Balance
        self.compressor_turbine_balance('25', '3', self.inputs.tit, self.inputs.tot)

        # Enforce LP Power Balance
        self.compressor_turbine_balance('2', '25', self.inputs.tot, self.inputs.tot2, comp_m_dot=self.inputs.m_dot + bypass_mdot)

        # LPT      
        self.turbine('4', '45', self.inputs.tot)

        # HPT
        self.turbine('45', '5', self.inputs.tot2)

        # Core nozzle
        self.A9, self.V9 = self.nozzle('5', '9')

        # Bypass Nozzle
        self.A19, self.V19 = self.nozzle('5', '9', gamma=self.inputs.gamma_c, m_dot=bypass_mdot)

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
    inputs.symbolise('fuel_air_ratio')
    inputs.symbolise('power_offtake')

    # Flight conditions
    flightcon = FlightCon(1828.8, 0.3)
    
    # Set up cycle calcs and solve
    engine = TurboJet(inputs, flightcon)
    engine.solve()

    print(engine.stations.T)
    print(engine.stations.P)
    print(engine.inputs.power_offtake)