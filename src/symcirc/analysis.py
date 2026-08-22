import os
import copy, time
from typing import List, Union

from symcirc import laplace
from symcirc.circuit import Circuit
from symcirc.component import Component, CurrentSource, VoltageSource, ANALYSIS_DEPENDENT
from symcirc.utils import *
from symcirc.matrix_constructors import TableauBuilder, TGNBuilder


class Analysis:
    """
    Main SymCirc class.
    When initialized it parses input netlist, and conducts the desired analysis which is then stored as a set
    of equations in the equation_matrix variable.

    :param Circuit circuit: Loaded netlist file which contains the circuit description.
        If you intend to load from a file, use the utils.load_file() function.
    :param str analysis_type: Analysis type identifier: "dc", "ac", "tf", "tran".
    :param bool symbolic: False if you want your results evaluated with numerical values from the netlist.

    :raise ValueError: If the analysis_type argument is invalid.
    """
    SOLVERS = ['gauss', 'lu', 'ddd']
    METHODS = ["tableau", "two_graph_node"]
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):

        self.analysis_type: str
        if use_symengine:
            os.environ["USE_SYMENGINE"] = "1"

        method = method.lower()
        if method not in self.METHODS:
            raise ValueError(f"Nonexistent analysis method: {method}\n Pick from: {self.METHODS}")

        if solver not in self.SOLVERS:
            raise ValueError(f"Nonexistent solver: {solver}\nPick from: {self.SOLVERS}")

        self.is_symbolic: bool = symbolic
        self.auto_eval = auto_eval
        self.precision: int = precision
        self.method: str = method
        self.sympy_ilt: bool = sympy_ilt
        self.solver = solver

        self.sbg = True

        self.circuit: Circuit = circuit
        self.node_voltage_identities: list = []

        self.symbol_dict = {}

        self.node_voltage_symbols: List[sympy.Symbol] = self.circuit.get_node_symbols()
        self.node_dict = self.circuit.get_node_dict()
        self.c_count = self.circuit.count_ports()
        self.node_count = self.circuit.count_nodes()

        self.eqn_matrix, self.solved_dict, self.symbols = self._analyse()  # solved_dict: {sympy.symbols(<vaviable_name>): <value>}
        self.symbol_dict = self.get_symbols()  # format: {<symbol_name> : <Symbol>}

    def v(self, name: str) -> sympy.Expr:
        """
        Returns the specified voltage
        """
        return self.component_voltage(name)

    def i(self, name: str) -> sympy.Expr:
        """
        Returns the specified current
        """
        return self.component_current(name)

    def get_all_results(self):
        results = self.node_voltages()
        results.update(self.component_values("all"))
        return results

    def get_symbols(self) -> Dict[str, sympy.Symbol]:
        #t0 = time.time()
        symbol_dict = {"s": s, "t": t, "f": f}
        for comp in self.circuit.components.values():
            if comp.type in ANALYSIS_DEPENDENT:
                sym_val = self._choose_source_val(comp, form="sym")
            else:
                sym_val = comp.sym_value

            if sym_val is None:
                continue
            elif hasattr(sym_val, "name"):
                symbol_dict[sym_val.name] = sym_val
            else:
                try:
                    for sym in sym_val.free_symbols:
                        symbol_dict[sym.name] = sym
                except AttributeError:
                    pass
        return symbol_dict

    def get_node_voltage_symbol(self, node: str) -> sympy.Symbol:
        return self.node_voltage_symbols[self.node_dict[node]]

    def component_voltage(self, name: str) -> sympy.Symbol:
        """Has to be implemented in child class"""
        ret = None
        v = f"v({name})"

        if self.method == "tableau":
            ret = self.solved_dict[sympy.symbols(v)]

        elif self.method == "two_graph_node":
            c = self.circuit.components[name]
            n1 = c.node1
            n2 = c.node2
            if n1 == "0":
                vn1 = 0
            else:
                vn1 = self.get_node_voltage(n1, force_s_domain=True)
            if n2 == "0":
                vn2 = 0
            else:
                vn2 = self.get_node_voltage(n2, force_s_domain=True)

            ret = sympy.cancel((vn1 - vn2))
        if ret is None:
            ret = sympy.Symbol(v, )
        return ret

    def component_current(self, name: str) -> sympy.Symbol:
        """Has to be implemented in child class"""
        ret = None
        i = f"i({name})"
        i_symbol = sympy.symbols(i)

        if self.method == "tableau":
            ret = self.solved_dict[i_symbol]

        if self.method == "two_graph_node":
            if i_symbol in self.solved_dict:
                return self.solved_dict[i_symbol]

            c = self.circuit.components[name]
            if c.type in ["r", "l", "c"]:
                n1 = c.node1
                n2 = c.node2
                if n1 == "0":
                    vn1 = 0
                else:
                    vn1 = self.get_node_voltage(n1, force_s_domain=True)
                if n2 == "0":
                    vn2 = 0
                else:
                    vn2 = self.get_node_voltage(n2, force_s_domain=True)

                if self.is_symbolic:
                    val = c.sym_value
                else:
                    val = c.value

                if c.type == "r":
                    ret = sympy.cancel((vn1 - vn2) / val)
                if c.type == "l":
                    ret = sympy.cancel((vn1 - vn2) / (val * s))
                if c.type == "c":
                    ret = sympy.cancel((vn1 - vn2) * val * s)
        return ret

    def get_node_voltage(self, node: str, force_s_domain=False) -> Union[sympy.Expr, None]:
        func = None
        try:
            func = self.solved_dict[self.node_voltage_symbols[self.node_dict[node]]]
        except KeyError:
            for identity in self.node_voltage_identities:
                if node in identity:
                    if "0" in identity:
                        return 0
                    else:
                        for n in identity:
                            try:
                                func = self.solved_dict[self.node_voltage_symbols[self.node_dict[n]]]
                            except:
                                pass
        return func

    def value_map(self, is_float=False):
        value_map = {}
        for comp in self.circuit.components.values():
            if is_float:
                if comp.type in ANALYSIS_DEPENDENT:
                    symbol = self._choose_source_val(comp, form="sym")
                    value_map[symbol] = self._choose_source_val(comp, form="float")
                elif comp.value_float is not None:
                    value_map[comp.sym_value] = comp.value_float
            else:
                if comp.type in ANALYSIS_DEPENDENT:
                    symbol = self._choose_source_val(comp, form="sym")
                    value_map[symbol] = self._choose_source_val(comp, form="rat")
                elif comp.value is not None:
                    value_map[comp.sym_value] = comp.value
        return value_map

    def _choose_source_val(self, c: Component, form="rat") -> sympy.Expr:
        """
        Has to be implemented in child class
        :param form: can be 'rat', 'float', 'sym'
        """
        pass

    def component_values(self, name: str = "all", default_python_datatypes: bool=False) -> Dict[str, sympy.Expr]:
        """
          Takes a string containing a single component name and returns a dictionary containing the voltage and current
            of the input component.

          :param str name: component id
          :param bool default_python_datatypes:
          :return dict ret: in format {"v(name)" : value, "i(name)" : value}
        """
        ret = {}
        if name == "all":
            ret = self.all_component_values(default_python_datatypes)
        else:
            v = f"v({name})"
            i = f"i({name})"
            ret[v] = self.component_voltage(name)
            ret[i] = self.component_current(name)
        return ret

    def all_component_values(self, default_python_datatypes=False) -> Dict[str, sympy.Expr]:
        """
          Returns a dictionary of all relevant voltages and currents in the circuit.

          :return dict ret: in format {"v(id1)" : value, "i(id2)" : value, ...}
        """
        ret = {}
        for key in self.circuit.components:
            if self.circuit.components[key].type == "k":
                pass
            elif self.circuit.components[key].name[-3:] == "_IC":
                pass
            else:
                if default_python_datatypes:
                    name = self.circuit.components[key].name
                    elem_dict = self.component_values(name)
                    try:
                        for elem_key in elem_dict:
                            ret[elem_key] = float(elem_dict[elem_key])
                    except TypeError:
                        ret.update(elem_dict)
                else:
                    name = self.circuit.components[key].name
                    elem_dict = self.component_values(name)
                    ret.update(elem_dict)
        return ret

    def node_voltages(self):
        """
          Returns a dictionary of all node voltages in the circuit.

          :return dict ret: in format {"v(node1)" : value, ...}
        """
        ret = {}
        for node in self.node_dict:
            if node.find("*ctrl") == -1:
                ret[f"v({node})"] = self.get_node_voltage(node)
            else: # Filter nodes added by virtual current sensors for CCC(V)S
                pass
        return ret

    def transfer_function(self, node1, node2):
        """
          Takes names of two nodes and returns their transfer function

          :param str node1: node id
          :param str node2: node id
          :return sympy_object ret: resulting transfer function
        """
        try:
            if node1 in [0, "0"]:
                voltage1 = 0
            else:
                v1 = sympy.Symbol(f"v({node1})")
                voltage1 = self.solved_dict[v1].simplify()
        except KeyError:
            raise ValueError(f"Node {node1} doesn't exist.")

        try:
            if node2 in [0, "0"]:
                voltage2 = 0
            else:
                v2 = sympy.Symbol(f"v({node2})")
                voltage2 = self.solved_dict[v2].simplify()
        except KeyError:
            raise ValueError(f"Node {node2} doesn't exist.")
        tf = (voltage2/voltage1)
        return tf


    def abstract_matrix(self, eqn_matrix, memo):
        i = 0
        for elem in eqn_matrix:
            new_elem = simplify_coeffs(elem, s, f"M{i}", memo)
            eqn_matrix[i] = new_elem
            i += 1

    def unabstract_results(self, solved_dict, memo):
        # Resubstitute
        inverted_memo = {v: k for k, v in memo.items()}
        for key in solved_dict:
            solved_dict[key] = solved_dict[key].subs(inverted_memo)

    def _analyse(self):
        """
          Implementation of all types of supported analysis.

            eqn_matrix, solved_dict, symbols
          :return sympy.Matrix eqn_matrix: matrix of the system equations
          :return dict solved_dict: dictionary of eqn_matrix solve results
          :return list symbols: list of all used sympy.symbol objects
        """

        # TODO: state machine to decide whether to use gauss, LU or DDD

        solved_dict = {}
        if self.solver == "gauss":
            memo = {}
            eqn_matrix, symbols = self._build_system_eqn()
            #self.abstract_matrix(eqn_matrix, memo)
            solved_dict = self._gauss_solve(eqn_matrix, symbols)
            #self.unabstract_results(solved_dict, memo)
        elif self.solver == "lu":
            eqn_matrix, symbols = self._build_system_eqn()
            solved_dict = self._lu_solve(eqn_matrix, symbols)
            for key in solved_dict:
                solved_dict[key] = general_simplify(solved_dict[key])
        elif self.solver == "ddd":
            # TODO: rework this hack into a more professional solution
            memo = {}
            tmp_symbolic = self.is_symbolic
            self.is_symbolic = True
            eqn_matrix, symbols = self._build_system_eqn()
            #self.abstract_matrix(eqn_matrix, memo)
            self.is_symbolic = tmp_symbolic
            solved_dict = self._ddd_solve(eqn_matrix, symbols)
            for key in solved_dict:
                solved_dict[key] = general_simplify(solved_dict[key])
            #self.unabstract_results(solved_dict, memo)

        t1 = time.time()
        #print(f"Matrix solve time by {self.solver}: {t1-t0}")

        # TODO: test whether this is a speedup or not and neccessity for readability
        #for key in solved_dict:
        #    solved_dict[key] = solved_dict[key].factor()
        #print(f"Gauss solve time: {time.time()-t0}")
        #print(solved_dict)
        # TODO: auto_eval implementation in expression retrieval; not directly here, causes issues with laplace
        #if self.auto_eval:
        #    solved_dict = {key: evalf(func, precision=self.precision) for key, func in solved_dict.items()}
        return eqn_matrix, solved_dict, symbols

    def _gauss_solve(self, eqn_matrix, symbols):
        solution_set = sympy.linsolve(eqn_matrix, *symbols)
        # EmptySet means there is no solution
        if solution_set is sympy.EmptySet:
            return {}
        # Extract the tuple from the FiniteSet
        solution = next(iter(solution_set))
        return dict(zip(symbols, solution))

    def _lu_solve(self, eqn_matrix, symbols):
        return sympy.solve_linear_system_LU(eqn_matrix, symbols)

    def _ddd_solve(self, eqn_matrix, symbols):
        from symcirc.ddd_solve import solve_cramer_ddd, solve_cramer_ddd_semi
        if self.is_symbolic:
            return solve_cramer_ddd(eqn_matrix, symbols, nested=True)
        else:
            return solve_cramer_ddd_semi(eqn_matrix, symbols, self.value_map(is_float=True))

    def _build_system_eqn(self):
        if self.method == "tableau":
            matrix_builder = TableauBuilder(self.circuit, self.analysis_type, self.is_symbolic)
            equation_matrix, symbols = matrix_builder.build()
        elif self.method == "two_graph_node":
            matrix_builder = TGNBuilder(self.circuit, self.analysis_type, self.is_symbolic)
            equation_matrix, symbols, node_voltage_identities = matrix_builder.build()
            self.node_voltage_identities = node_voltage_identities
        return equation_matrix, symbols


class DC(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        self.analysis_type = "dc"
        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)

    def _analyse(self):
        eqn_matrix, solved_dict, symbols = super()._analyse()
        for sym in symbols:
            try:
                solved_dict[sym] = sympy.limit(solved_dict[sym], s, 0)
            except KeyError:
                pass
            except TypeError:
                pass
        return eqn_matrix, solved_dict, symbols

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource], form="rat") -> sympy.Expr:
        return c.dc_val(form)

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = super().component_current(name)
        if ret is None:
            c = self.circuit.components[name]
            if c.type == "i":
                if self.is_symbolic:
                    ret = c.dc_val(form="sym")
                else:
                    ret = c.dc_val(form="rat")
        else:
            ret = sympy.limit(ret, s, 0)
        return ret


class TF(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        self.analysis_type = "tf"
        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)

    def _analyse(self):
        eqn_matrix, solved_dict, symbols = super()._analyse()
        return eqn_matrix, solved_dict, symbols

    def component_current(self, name: str) -> sympy.Symbol:
        ret = super().component_current(name)
        if ret is None:
            c = self.circuit.components[name]
            if c.type == "i":
                if self.is_symbolic:
                    ret = c.ac_val(form="sym")
                else:
                    ret = c.ac_val(form="rat")
        return ret

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource], form="rat") -> sympy.Expr:
        return c.tf_val(form)


class AC(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        self.analysis_type = "ac"
        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)

    def _analyse(self):
        eqn_matrix, solved_dict, symbols = super()._analyse()
        for sym in symbols:
            try:
                solved_dict[sym] = solved_dict[sym].subs(s, 2 * pi * f * j)
            except:
                pass
        return eqn_matrix, solved_dict, symbols

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource], form="rat") -> sympy.Expr:
        return c.ac_val(form)

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = super().component_current(name)
        if ret is None:
            c = self.circuit.components[name]
            if c.type == "i":
                if self.is_symbolic:
                    ret = c.ac_val(form="sym")
                else:
                    ret = c.ac_val(form="rat")
        else:
            ret = ret.subs(s, 2*pi*f*j)
        return ret


class ACNumeric(AC):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        self.analysis_type = "ac_num"
        if use_symengine:
            os.environ["USE_SYMENGINE"] = "1"

        method = method.lower()
        if method not in ["tableau", "two_graph_node"]:
            raise ValueError(f"Nonexistent analysis method: {method}")

        self.is_symbolic: bool = symbolic
        self.auto_eval = auto_eval
        self.precision: int = precision
        self.method: str = method
        self.sympy_ilt: bool = sympy_ilt

        self.sbg = True

        self.circuit: Circuit = circuit
        self.node_voltage_identities: list = []

        self.symbol_dict = {}

        self.node_voltage_symbols: List[sympy.Symbol] = self.circuit.get_node_symbols()
        self.node_dict = self.circuit.get_node_dict()
        self.c_count = self.circuit.count_ports()
        self.node_count = self.circuit.count_nodes()

        self.eqn_matrix, self.symbols = self._build_system_eqn()


    def run(self, freqs):
        import numpy as np

        ac_eqn_matrix = self.eqn_matrix.subs(s, 2*j*pi*f)
        A_sym = ac_eqn_matrix[:, :-1]  # all columns except last
        b_sym = ac_eqn_matrix[:, -1]  # last column

        fA = sympy.lambdify(f, A_sym, "numpy")
        fb = sympy.lambdify(f, b_sym, "numpy")

        symbols_str = [str(sym) for sym in self.symbols]
        results = {sym: [] for sym in symbols_str}

        for freq in freqs:
            x = np.linalg.solve(fA(freq), fb(freq)).flatten()

            for sym, val in zip(symbols_str, x):
                results[sym].append(val)

        return results


class TRAN(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        self.analysis_type = "tran"
        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)

    def get_node_voltage(self, node: str, force_s_domain=False) -> Union[sympy.Expr, None]:
        """
        Old way to return a node voltage, will be deprecated soon
        """
        func = super().get_node_voltage(node)
        if func is None or force_s_domain:
            return func
        else:
            res = laplace.iLT(func, self.sympy_ilt)
            return res

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource], form="rat") -> sympy.Expr:
        return c.tran_val(form)

    def component_voltage(self, name: str) -> sympy.Expr:
        """
        Old way to return a component voltage, will be deprecated soon
        """
        ret = super().component_voltage(name)
        if ret is not None:
            ret = laplace.iLT(ret, self.sympy_ilt)
            ret = sympy.factor_terms(ret)
        return ret

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = super().component_current(name)
        if ret is None:
            c = self.circuit.components[name]
            if c.type == "i":
                if self.is_symbolic:
                    ret = c.tran_val(form="sym")
                else:
                    ret = c.tran_val(form="rat")
        if ret is not None:
            ret = laplace.iLT(ret, self.sympy_ilt)
            ret = sympy.factor_terms(ret)
        return ret


def AnalyseCircuit(netlist: str, analysis_type: str = "DC", method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, operating_point: Union[Dict[str, float], None] = None,
                   solver: str = "gauss") -> Analysis:

    circuit = Circuit(netlist, operating_point)
    analysis_type = analysis_type.lower()
    if analysis_type == "dc":
        analysis = DC(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)
    elif analysis_type == "ac":
        analysis = AC(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)
    elif analysis_type == "tf":
        analysis = TF(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)
    elif analysis_type == "tran":
        analysis = TRAN(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)
    else:
        raise ValueError(f"Nonexistent analysis type: {analysis_type}")
    return analysis