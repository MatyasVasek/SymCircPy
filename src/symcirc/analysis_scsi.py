from symcirc.analysis import Circuit, Analysis
import os
import operator
import copy, time
from typing import List, Set, Dict, Union, Any
from collections import defaultdict

from sympy import Symbol
from symcirc.utils import *

class SwitchedAnalysis(Analysis):
    def __init__(self, circuit: Circuit, phases, mode: str = "sc", method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool =False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):

        # TODO: clean up
        self.analysis_type = "TF"
        self.phases, self.frequency = self.parse_phases(phases)
        self.mode = mode
        self.method = method
        self.initial_values()

        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine)

    def initial_values(self):
        """
        Checks if __init__ parameters are valid for SC/SI.
        """
        if self.phases != "undefined" and self.method != "two_graph_node":
            raise ValueError("Valid method for SCSI circuits is: 'two_graph_node'")

    def parse_phases(self, phases): # used by SCSI analysis
        """
        Parameter phases should be in a correct syntax:
        'P=integer', P=(integer,frequency),'P=[...]' or P=([...],frequency)
        Returns phases list with number of phases (phases[0]) and
        phase fractions which signify how long does a phase last (phases[1] to phases[n]).
        Also returns frequency.
        Implemented by Fishecko
        """
        frequency = 1
        if phases != "undefined":
            phase_definition = []
            phase_def_syntax_error = \
                (
                    "Invalid phase definition syntax, use 'P=integer', P=(integer,frequency),'P=[...]' or P=([...],frequency), \n"
                    "where the integer is the number of phases "
                    "and the list contains lengths of phases written as a fraction (fractions must add up to 1)")
            if phases.startswith("P="):
                phases = phases.replace("P=", "")
                if phases.startswith("(") and phases.endswith(")"):
                    phases = phases.replace("(", "")
                    phases = phases.replace(")", "")
                    if phases.startswith("[") and "]" in phases:
                        phases = phases.replace("[", "")
                        phases = phases.replace("]", "")
                        phase_definition = sympy.sympify(phases.split(','))
                        frequency = phase_definition[-1]
                        del phase_definition[-1]
                        phase_sum = sum(phase_definition)
                        if phase_sum != 1:
                            raise ValueError("The sum of phase lengths must be 1")
                        else:
                            phase_definition.insert(0, len(phase_definition))
                    elif len(phases.split(',')) == 2:
                        phases = sympy.sympify(phases.split(','))
                        if int(phases[0]) < 2:
                            raise ValueError("The number of phases can't be less than 2")
                        else:
                            phase_definition.append(phases[0])
                            frequency = phases[-1]
                            for i in range(phase_definition[0]):
                                phase_definition.append(sympy.sympify("1/" + str(phase_definition[0])))
                    else:
                        raise SyntaxError(phase_def_syntax_error)
                elif phases.startswith("[") and phases.endswith("]"):
                    phases = phases.replace("[", "")
                    phases = phases.replace("]", "")
                    phase_definition = sympy.sympify(phases.split(','))
                    phase_sum = sum(phase_definition)
                    if phase_sum != 1:
                        raise ValueError("The sum of phase lengths must be 1")
                    else:
                        phase_definition.insert(0, len(phase_definition))
                elif "," not in phases:
                    if int(phases) < 2:
                        raise ValueError("The number of phases can't be less than 2")
                    else:
                        phase_definition.append(int(phases))
                        for i in range(phase_definition[0]):
                            phase_definition.append(sympy.sympify("1/" + str(phase_definition[0])))
                else:
                    raise SyntaxError(phase_def_syntax_error)

            else:
                raise SyntaxError(phase_def_syntax_error)
            phases = phase_definition
        return phases, frequency

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
            temp = {}
            temp.update(self.SCSI_component_values(name))
            ret = self.SCSI_result_formatter(temp)
            # TODO result_formatter isn't working well, type error
        return ret

    def all_component_values(self, default_python_datatypes=False):
        """
          Returns a dictionary of all relevant voltages and currents in the circuit.

          :return dict ret: in format {"v(id1)" : value, "i(id2)" : value, ...}
        """
        ret = {}
        temp = {}
        for key in self.circuit.components:
            temp.update(self.SCSI_component_values(self.circuit.components[key].name))
        ret = self.SCSI_result_formatter(temp)
        return ret

    def node_voltages(self):
        """
          Returns a dictionary of all node voltages in the circuit.

          :return dict ret: in format {"v(node1)" : value, ...}
        """
        ret = {}
        num_of_phases = self.phases[0]
        for phase in range(1, num_of_phases + 1):
            for node in self.node_dict:
                ret[sympy.symbols(f"v({node})_{phase}")] = self.SCSI_get_node_voltage(node, phase)
        ret = self.SCSI_result_formatter(ret)
        return ret

    def _analyse(self):
        # TODO numerical eval for SCSI - isn't too good, also leave exponents in rational form
        eqn_matrix, symbols = self._build_system_eqn()
        solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
        self.SCSI_symbol_z_factor(solved_dict)
        self.SCSI_z_pow_inv_sub(solved_dict)
        # if self.analysis_type == "TF":
        if not self.is_symbolic:
            for sym in symbols:
                try:
                    for name in self.circuit.components:
                        c = self.circuit.components[name]
                        if c.type in ["v", "i"]:
                            ac_val = c.ac_val("rat")
                            if ac_val:
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, ac_val)
                        else:
                            if c.value:
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                    solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                except KeyError:
                    pass
        '''elif self.analysis_type == "AC":
            f = self.frequency
            if self.is_symbolic:
                for sym in symbols:
                    try:
                        solved_dict[sym] = solved_dict[sym].subs(z, sympy.exp(1) ** (2 * sympy.pi * f * j))
                    except KeyError:
                        pass
            if not self.is_symbolic:
                for sym in symbols:
                    try:
                        solved_dict[sym] = solved_dict[sym].subs(z, sympy.exp(1) ** (2 * sympy.pi * f * j))
                        for name in self.components:
                            c = self.components[name]
                            if c.type in ["v", "i"]:
                                if c.ac_value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                            else:
                                if c.value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                        solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                    except KeyError:
                        pass
        elif self.analysis_type == "tran":
            raise ValueError("Transient analysis not yet implemented.")
            # TODO transient analysis - z-transform and quantization (floor function)'''

        return eqn_matrix, solved_dict, symbols

    def _build_system_eqn(self):

        # TODO switched currents, it's not working now
        sc_value_error = "Resistors, inductors, current sources, CCVSs and VCCSs aren't permitted in switched capacitor mode"
        num_of_phases = self.phases[0]
        symbols = []
        v_graph_collapses = []
        i_graph_collapses = []
        v_graph_nodes = []
        i_graph_nodes = []
        matrix_col_expand = 0

        for key in self.circuit.components:
            c = self.circuit.components[key]
            if c.type == "h":
                raise ValueError("Current controlled voltage sources aren't permitted in SC/SI mode")
            if c.type == "l":
                raise ValueError("Inductors are not permitted in SC/SI mode")
            if self.mode == "scideal":
                if c.type in ["r", "i", "g"]:
                    raise ValueError(sc_value_error)
            # else:
            #     if c.type == "v":
            #         raise ValueError("Inductors, voltage sources and CCVSs aren't permitted in switched current mode")


        for key in self.circuit.components:
            c = self.circuit.components[key]
            if c.type == "c":
                for phase in range(1, num_of_phases + 1):
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
            if c.type == "v":
                for phase in range(1, num_of_phases + 1):
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_collapse_tgn(i_graph_collapses, c.node1 + '_' + str(phase), c.node2 + '_' + str(phase))
            if c.type == "s":
                self.SCSI_graph_append_tgn(c.node1 + '_' + str(c.phase), v_graph_nodes)
                self.SCSI_graph_append_tgn(c.node2 + '_' + str(c.phase), v_graph_nodes)
                self.SCSI_graph_append_tgn(c.node1 + '_' + str(c.phase), i_graph_nodes)
                self.SCSI_graph_append_tgn(c.node2 + '_' + str(c.phase), i_graph_nodes)
                self.SCSI_collapse_tgn(i_graph_collapses, c.node1 + '_' + str(c.phase), c.node2 + '_' + str(c.phase))
            if c.type == "a":
                for phase in range(1, num_of_phases + 1):
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_collapse_tgn(i_graph_collapses, c.node1 + '_' + str(phase), c.node2 + '_' + str(phase))
                    self.SCSI_collapse_tgn(v_graph_collapses, c.node3 + '_' + str(phase), c.node4 + '_' + str(phase))
            if c.type == "e":
                for phase in range(1, num_of_phases + 1):
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_collapse_tgn(i_graph_collapses, c.node1 + '_' + str(phase), c.node2 + '_' + str(phase))
            if c.type == "f": # TODO: fix after bug #8 fix
                c_v = self.circuit.components[c.current_sensor]
                for phase in range(1, num_of_phases + 1):
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c_v.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c_v.shorted_node + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c_v.node2 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c_v.shorted_node + '_' + str(phase), i_graph_nodes)
                    self.SCSI_collapse_tgn(v_graph_collapses, c_v.node2 + '_' + str(phase),
                                           c_v.shorted_node + '_' + str(phase))
                matrix_col_expand += num_of_phases
            if c.type == "r":
                for phase in range(1, num_of_phases + 1):
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
            if c.type == "i":
                for phase in range(1, num_of_phases + 1):
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
            if c.type == "g":
                for phase in range(1, num_of_phases + 1):
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), i_graph_nodes)

        i_graph_collapses = self.SCSI_collapse_list_collapser(i_graph_collapses)
        v_graph_collapses = self.SCSI_collapse_list_collapser(v_graph_collapses)

        for collapse_list in i_graph_collapses:
            i = 0
            tmp_i_graph_nodes = copy.copy(i_graph_nodes)
            for n in tmp_i_graph_nodes:
                if n in collapse_list:
                    if any(node.startswith('0') for node in collapse_list):
                        i_graph_nodes.remove(n)
                    else:
                        i_graph_nodes[i] = min(collapse_list)
                i += 1
        for collapse_list in v_graph_collapses:
            i = 0
            tmp_v_graph_nodes = copy.copy(v_graph_nodes)
            for n in tmp_v_graph_nodes:
                if n in collapse_list:
                    if any(node.startswith('0') for node in collapse_list):
                        v_graph_nodes.remove(n)
                    else:
                        v_graph_nodes[i] = min(collapse_list)
                i+=1
        self.node_voltage_identities = v_graph_collapses

        v_graph_nodes = list(set(v_graph_nodes))
        i_graph_nodes = list(set(i_graph_nodes))

        v_graph_nodes = self.SCSI_node_list_sort(v_graph_nodes)
        i_graph_nodes = self.SCSI_node_list_sort(i_graph_nodes)

        m_size = len(v_graph_nodes) + matrix_col_expand
        M = sympy.Matrix(sympy.zeros(m_size))
        S = sympy.Matrix(sympy.zeros(m_size, 1))
        index_row = 0
        index_col = 0
        symbols_to_append = []

        for key in self.circuit.components:
            c = self.circuit.components[key]
            if c.type == "c":
                self.SCSI_add_capacitor_tgn(M, v_graph_nodes, i_graph_nodes,
                                            c, i_graph_collapses, v_graph_collapses,
                                            num_of_phases, matrix_col_expand)
            if c.type == "v":
                index_row = self.SCSI_add_voltage_source_tgn(M, S, v_graph_nodes, i_graph_nodes,
                                                             c, index_row, v_graph_collapses, num_of_phases)
            if c.type == "s":
                index_row = self.SCSI_add_switch_tgn(M, v_graph_nodes, i_graph_nodes,
                                                     c, index_row, v_graph_collapses, num_of_phases)
            if c.type == "a":
                pass
            if c.type == "e":
                index_row = self.SCSI_add_VVT_tgn(M, v_graph_nodes, i_graph_nodes,
                                                  c, index_row, v_graph_collapses, num_of_phases)
            if c.type == "f":
                pass
            if c.type == "r":
                self.SI_add_resistor_tgn(M, v_graph_nodes, i_graph_nodes,
                                         c, i_graph_collapses, v_graph_collapses, num_of_phases)
            if c.type == "i":
                self.SI_add_current_source_tgn(S, i_graph_nodes, c, i_graph_collapses, num_of_phases)
            if c.type == "g":
                self.SI_add_VCT_tgn(M, v_graph_nodes, i_graph_nodes,
                                    c, i_graph_collapses, v_graph_collapses, num_of_phases)

        for key in self.circuit.components:
            c = self.circuit.components[key]
            if c.type == "f":
                index_col = self.SCSI_add_QQT_tgn(M, v_graph_nodes, i_graph_nodes,
                                                  c, index_col, i_graph_collapses, num_of_phases)
                for phase in range(1, num_of_phases + 1):
                    if self.mode == "scideal":
                        symbols_to_append.append(sympy.Symbol(f"q({c.current_sensor})_{phase}"))
                    else:
                        symbols_to_append.append(sympy.Symbol(f"i({c.current_sensor})_{phase}"))

        equation_matrix = M.col_insert(m_size, S)

        for node in v_graph_nodes:
            symbols.append(sympy.Symbol(f"v({node[0:-2]}){node[-2:]}"))
        for symb in symbols_to_append:
            symbols.append(symb)

        return equation_matrix, symbols

    def SCSI_node_voltage_symbols(self):
        num_of_phases = self.phases[0]
        voltage_symbol_list = []
        for phases in range(num_of_phases):
            for node in self.node_dict:
                if node != "0":
                    voltage_symbol_list.append(sympy.Symbol("v({node})_{phases}".format(node=node, phases=phases + 1)))
        return voltage_symbol_list

    def SCSI_symbol_list_order(self, node_symbols, charge_symbols):
        charge_length = len(charge_symbols)
        node_length = len(node_symbols)
        num_of_phases = self.phases[0]
        symbols = []
        for phase in range(1, num_of_phases + 1):
            for symbol in range(node_length):
                if int(str(node_symbols[symbol])[-1]) == phase:
                    symbols.append(node_symbols[symbol])
            for symbol in range(charge_length):
                if int(str(charge_symbols[symbol])[-1]) == phase:
                    symbols.append(charge_symbols[symbol])
        return symbols

    def SCSI_symbol_z_factor(self, symbols):
        num_of_phases = self.phases[0]
        phase_helper_array = []
        for phase in range(num_of_phases):
            phase_helper_array.append(self.phases[phase + 1])
        for symbol in symbols:
            symbol_phase = int(str(symbol)[-1])
            temp = 0
            for phase in range(symbol_phase - 1):
                temp += phase_helper_array[phase]
            if symbol_phase != num_of_phases:
                symbols[symbol] *= z ** (- temp)
            else:
                symbols[symbol] *= z ** phase_helper_array[-1]

    def SCSI_z_pow_inv_sub(self, solved_dict):
        for expression in solved_dict:
            solved_dict[expression] = solved_dict[expression].subs(self.symbol_dict)

    def SCSI_z_power_substitution(self, symbol):
        dict_length = len(self.symbol_dict)
        x = sympy.Symbol("x_{index}".format(index=dict_length))
        self.symbol_dict[x] = symbol
        return x

    def SCSI_graph_append_tgn(self, node, graph):
        if node.startswith('0'):
            pass
        elif node not in graph:
            graph.append(node)

    def SCSI_collapse_tgn(self, graph_collapses, node1, node2):
        collapsed = False
        for collapse_list in graph_collapses:
            if node1 in collapse_list:
                collapsed = True
                collapse_list.append(node2)  # set node2 to be collapsed into node1 on the current graph
            elif node2 in collapse_list:
                collapsed = True
                collapse_list.append(node1)  # set node1 to be collapsed into node2 on the current graph
        if not collapsed:
            graph_collapses.append([node1, node2])  # set node2 to be collapsed into node1 on the current graph

    def SCSI_collapse_list_collapser(self, graph_collapses):
        for collapse_list in graph_collapses:
            for node in collapse_list:
                for cl in graph_collapses:
                    if node in cl and cl!=collapse_list:
                        for n in cl:
                            if n not in collapse_list:
                                collapse_list.append(n)
        for collapse_list in graph_collapses:
            collapse_list.sort()

        temp = []
        for collapse_list in graph_collapses:
            if collapse_list not in temp:
                temp.append(collapse_list)

        return temp

    def SCSI_node_list_sort(self, nodes):
        helper = []
        i = 0
        for element in nodes:
            helper.append([])
            helper[i].append(element[:element.index('_')])
            helper[i].append(element[element.index('_') + 1:])
            i += 1
        helper.sort(key=operator.itemgetter(1, 0))
        nodes_new = []
        for element in helper:
            nodes_new.append(element[0] + '_' + element[1])
        return nodes_new

    def SCSI_index_tgn(self, nodes, node, collapses):
        try:
            return nodes.index(node)
        except ValueError:
            for collapse_list in collapses:
                if node in collapse_list:
                    if any(n.startswith('0') for n in collapse_list):
                        return None
                    else:
                        return nodes.index(min(collapse_list))
                elif node.startswith("0"):
                    return None

    def SCSI_add_capacitor_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses, num_of_phases, matrix_col_expand):
        phase_offset = int((M.shape[0] - matrix_col_expand) / num_of_phases)
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(phase), v_graph_collapses)
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
            if n1v is not None:
                if n1i is not None:
                    if phase == num_of_phases:
                        M[n1i, n1v] += +val*z
                    else:
                        M[n1i, n1v] += +val
                    if phase == 1:
                        M[n1i, n1v - phase_offset - matrix_col_expand] += -val
                    else:
                        M[n1i, n1v - phase_offset] += -val
                if n2i is not None:
                    if phase == num_of_phases:
                        M[n2i, n1v] += -val*z
                    else:
                        M[n2i, n1v] += -val
                    if phase == 1:
                        M[n2i, n1v - phase_offset - matrix_col_expand] += +val
                    else:
                        M[n2i, n1v - phase_offset] += +val
            if n2v is not None:
                if n1i is not None:
                    if phase == num_of_phases:
                        M[n1i, n2v] += -val*z
                    else:
                        M[n1i, n2v] += -val
                    if phase == 1:
                        M[n1i, n2v - phase_offset - matrix_col_expand] += +val
                    else:
                        M[n1i, n2v - phase_offset] += +val
                if n2i is not None:
                    if phase == num_of_phases:
                        M[n2i, n2v] += +val*z
                    else:
                        M[n2i, n2v] += +val
                    if phase == 1:
                        M[n2i, n2v - phase_offset - matrix_col_expand] += -val
                    else:
                        M[n2i, n2v - phase_offset] += -val

    def SCSI_add_voltage_source_tgn(self, M, S, v_nodes, i_nodes, c, index_row, v_graph_collapses, num_of_phases):
        if self.is_symbolic:
            val = c.sym_value
        else:
            # self.analysis_type is always "TF" for SwitchedAnalysis
            val = c.tf_val("rat")
        phase_helper_array = []
        z_symbol_array = []
        temp = 0
        for f in range(num_of_phases):
            phase_helper_array.append(temp)
            temp += self.phases[f + 1]
            z_symbol_array.append(z ** sympy.sympify(phase_helper_array[f]))
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(phase), v_graph_collapses)
            row = len(i_nodes) + index_row
            if n1v is not None:
                if phase == num_of_phases:
                    M[row, n1v] += z
                else:
                    M[row, n1v] += 1
            if n2v is not None:
                if phase == num_of_phases:
                    M[row, n2v] += -z
                else:
                    M[row, n2v] += -1
            S[row, 0] += val * self.SCSI_z_power_substitution(z_symbol_array[phase - 1])
            index_row += 1
        return index_row

    def SCSI_add_switch_tgn(self, M, v_nodes, i_nodes, c, index_row, v_graph_collapses, num_of_phases):
        n1v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(c.phase), v_graph_collapses)
        n2v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(c.phase), v_graph_collapses)
        row = len(i_nodes) + index_row
        if n1v is not None:
            if int(c.phase) == num_of_phases:
                M[row, n1v] += z
            else:
                M[row, n1v] += 1
        if n2v is not None:
            if int(c.phase) == num_of_phases:
                M[row, n2v] += -z
            else:
                M[row, n2v] += -1
        index_row += 1
        return index_row

    def SCSI_add_VVT_tgn(self, M, v_nodes, i_nodes, c, index_row, v_graph_collapses, num_of_phases):
        if self.is_symbolic:
            e = c.sym_value
        else:
            e = c.value
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node3 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node4 + '_' + str(phase), v_graph_collapses)
            n3v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(phase), v_graph_collapses)
            n4v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(phase), v_graph_collapses)
            row = len(i_nodes) + index_row
            if n1v is not None:
                if phase == num_of_phases:
                    M[row, n1v] += -e*z
                else:
                    M[row, n1v] += -e
            if n2v is not None:
                if phase == num_of_phases:
                    M[row, n2v] += e*z
                else:
                    M[row, n2v] += e
            if n3v is not None:
                if phase == num_of_phases:
                    M[row, n3v] += z
                else:
                    M[row, n3v] += 1
            if n4v is not None:
                if phase == num_of_phases:
                    M[row, n4v] += -z
                else:
                    M[row, n4v] += -1
            index_row += 1
        return index_row

    def SCSI_add_QQT_tgn(self, M, v_nodes, i_nodes, c, index_col, i_graph_collapses, num_of_phases):
        if self.is_symbolic:
            f = c.sym_value
        else:
            f = c.value
        c_v = self.circuit.components[c.current_sensor]
        for phase in range(1, num_of_phases + 1):
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
            # TODO: fix after bug #8 fix
            n3i = self.SCSI_index_tgn(i_nodes, c_v.node2 + '_' + str(phase), i_graph_collapses)
            n4i = self.SCSI_index_tgn(i_nodes, c_v.shorted_node + '_' + str(phase), i_graph_collapses)
            col = len(v_nodes) + index_col
            if n1i is not None:
                if phase == num_of_phases:
                    M[n1i, col] += f*z
                else:
                    M[n1i, col] += f
            if n2i is not None:
                if phase == num_of_phases:
                    M[n2i, col] += -f*z
                else:
                    M[n2i, col] += -f
            if n3i is not None:
                if phase == num_of_phases:
                    M[n3i, col] += z
                else:
                    M[n3i, col] += 1
            if n4i is not None:
                if phase == num_of_phases:
                    M[n4i, col] += -z
                else:
                    M[n4i, col] += -1
            index_col += 1
        return index_col

    def SI_add_resistor_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses, num_of_phases):
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        y = 1 / val
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(phase), v_graph_collapses)
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
            if n1v is not None:
                if n1i is not None:
                    M[n1i, n1v] += +y
                if n2i is not None:
                    M[n2i, n1v] += -y
            if n2v is not None:
                if n1i is not None:
                    M[n1i, n2v] += -y
                if n2i is not None:
                    M[n2i, n2v] += +y

    def SI_add_current_source_tgn(self, S, i_nodes, c, i_graph_collapses, num_of_phases):
        if self.is_symbolic:
            val = c.sym_value
        else:
            # self.analysis_type is always "TF" for SwitchedAnalysis
            val = c.tf_val("rat")
        for phase in range(1, num_of_phases + 1):
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
            if n1i is not None:
                S[n1i, 0] += -val
            if n2i is not None:
                S[n2i, 0] += val

    def SI_add_VCT_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses, num_of_phases):
        if self.is_symbolic:
            g = c.sym_value
        else:
            g = c.value
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node3 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node4 + '_' + str(phase), v_graph_collapses)
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
            if n1v is not None:
                if n1i is not None:
                    M[n1i, n1v] += +g
                if n2i is not None:
                    M[n2i, n1v] += -g
            if n2v is not None:
                if n1i is not None:
                    M[n1i, n2v] += -g
                if n2i is not None:
                    M[n2i, n2v] += +g

    def SCSI_get_node_voltage(self, node, phase, force_z = False):
        """
        Returns node voltage of the specified node in specified phase.
        """
        phase_string = str(phase)
        formatted_node = str(node) + "_" + phase_string
        formatted_voltage = sympy.symbols("v(" + str(node) + ")_" + phase_string)
        try:
            value = self.solved_dict[formatted_voltage]
        except KeyError:
            for identity in self.node_voltage_identities:
                if formatted_node in identity:
                    if any(identity_node.startswith("0") for identity_node in identity):
                        value = 0
                    else:
                        for identity_node in identity:
                            try:
                                if str(identity_node).count("_") == 1:
                                    formatted_voltage = sympy.symbols("v(" + identity_node.split("_")[0]
                                                                      + ")_" + phase_string)
                                else:
                                    formatted_voltage = sympy.symbols("v(" + identity_node.split(")_")[0]
                                                                      + "))_" + phase_string)
                                value = self.solved_dict[formatted_voltage]
                            except KeyError:
                                pass
        if force_z:
            return value
        elif self.analysis_type == "tran":
            pass
            #return z_transform.iZT(value)
        else:
            return value

    def SCSI_component_voltage(self, name, phase):
        """
        Input - name of the component and the particular phase we want the value for.
        Returns the voltage of a component.
        """
        value_dict = {}
        c = self.circuit.components[name]
        biquad_input_key = sympy.symbols(f"v({name}-in)_{phase}")
        if c.node1 == "0":
            node1_value = 0
        else:
            node1_value = self.SCSI_get_node_voltage(c.node1, phase, force_z=True)
        if c.node2 == "0":
            node2_value = 0
        else:
            node2_value = self.SCSI_get_node_voltage(c.node2, phase, force_z=True)
        if c.type in ["a", "e", "g"]:
            if c.node3 == "0":
                node3_value = 0
            else:
                node3_value = self.SCSI_get_node_voltage(c.node3, phase, force_z=True)
            if c.node4 == "0":
                node4_value = 0
            else:
                node4_value = self.SCSI_get_node_voltage(c.node4, phase, force_z=True)
            value_dict[biquad_input_key] = sympy.cancel(node3_value - node4_value)
        if c.type == "f":
            c_v = self.circuit.components[c.current_sensor]
            if c_v.node2 == "0":
                cv_node_value = 0
            else:
                cv_node_value = self.SCSI_get_node_voltage(c_v.node2, phase, force_z=True)
            if c_v.shorted_node == "0": # TODO: fix after bug #8 fix
                cv_short_value = 0
            else:
                cv_short_value = self.SCSI_get_node_voltage(c_v.shorted_node, phase, force_z=True) # TODO: fix after bug #8 fix
            value_dict[biquad_input_key] = sympy.cancel(cv_node_value - cv_short_value)
        # if c.type in ["a", "e", "f", "g"]:
        #     value_dict[sympy.symbols(f"v({name}-out)_{phase}")] = sympy.cancel(node1_value - node2_value)
        # else:
        value_dict[sympy.symbols(f"v({name})_{phase}")] = sympy.cancel(node1_value - node2_value)
        # if self.analysis_type == "tran":
        #     for entry in value_dict:
        #         value_dict[entry] = z_transform.IZT(value_dict[entry])
        return value_dict

    def SCSI_component_charge(self, name, phase):
        """
        Input - name of the component and the particular phase we want the value for.
        Returns a charge (or a current if SI).
        """
        #TODO evaluate charges/currents for two_graph_node...somehow
        value_dict = {}
        if self.mode == "scideal":
            charge = sympy.symbols(f"q({name})_{phase}")
        else:
            charge = sympy.symbols(f"i({name})_{phase}")
        c = self.circuit.components[name]
        voltage_key = sympy.symbols(f"v({name})_{phase}")
        if self.is_symbolic:
            value = c.sym_value
        else:
            value = c.value
        if c.type in ["c", "r"]:
            if c.type == "r":
                value_dict[charge] = sympy.cancel(self.SCSI_component_voltage(name, phase)[voltage_key] / value)
            if c.type == "c":
                value_dict[charge] = sympy.cancel(self.SCSI_component_voltage(name, phase)[voltage_key] * value)
        elif c.type == "i":
            if self.is_symbolic:
                value = c.sym_value
            else:
                # self.analysis_type is always "TF" for SwitchedAnalysis
                value = c.tf_val("rat")
            value_dict[charge] = value
        elif c.type == "g":
            value_dict[charge] = sympy.cancel(self.SCSI_component_voltage(name, phase)[voltage_key] * value)
        elif c.type == "f":
            if self.mode == "scideal":
                control_charge = sympy.symbols(f"q({c.current_sensor})_{phase}")
            else:
                control_charge = sympy.symbols(f"i({c.current_sensor})_{phase}")
            value_dict[charge] = sympy.cancel(self.solved_dict[control_charge] * value)
        else:
            value_dict[charge] = charge
        # if self.analysis_type == "tran":
        #     for entry in value_dict:
        #         value_dict[entry] = z_transform.IZT(value_dict[entry])
        return value_dict

    def SCSI_component_values(self, name):
        """
        Input - name of a component.
        Returns a dictionary of all values of the component in all phases.
        """
        num_of_phases = self.phases[0]
        value_dict = {}
        for phase in range(1, num_of_phases + 1):
            value_dict.update(self.SCSI_component_voltage(name, phase))
            value_dict.update(self.SCSI_component_charge(name, phase))
            #TODO TypeError: 'Mul' object is not iterable
        return value_dict

    def SCSI_result_formatter(self, dict):
        """
        Input is a dictionary of values which aren't sorted by component and node quantities.
        Returns a dictionary of lists - every component and node quantity has
        more phase values.
        """
        #TODO delete unimportant data - OpAmp input etc.
        temp = {}
        for key in dict:
            if list(str(key)).count("_") == 1:
                name = str(key).split("_")[0]
            else:
                name = str(key).split(")_")[0]
                name = "".join([name, ')'])
            temp[name] = []
            for i in range(self.phases[0]):
                temp[name].append(0)
        for key in dict:
            if list(str(key)).count("_") == 1:
                name = str(key).split("_")[0]
                phase = int(str(key).split("_")[1])
            else:
                name = str(key).split(")_")[0]
                name = "".join([name, ')'])
                phase = int(str(key).split(")_")[1])
            temp[name][phase - 1] = dict[key]
        return temp