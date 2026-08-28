import copy
import sympy
from symcirc.utils import s

class TGNBuilder:
    def __init__(self, circuit, analysis_type, is_symbolic):
        self.analysis_type = analysis_type
        self.is_symbolic = is_symbolic
        self.circuit = circuit

    def build(self):
        symbols = []
        v_graph_collapses = []
        i_graph_collapses = []
        v_graph_nodes = []
        i_graph_nodes = []
        matrix_col_expand = 0
        for key in self.circuit.components:
            c = self.circuit.components[key]

            if c.type == "k":
                matrix_col_expand += 2

            if c.type in ["r", "l", "c"]:
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)

            if c.type == "v":
                # Modified topology to also represent current through the source
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                matrix_col_expand += 1

            if c.type == "i":
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)

            if c.type == "g":  # VCT/VCCS
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node3, v_graph_nodes)
                self.graph_append(c.node4, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self.graph_append(c.node3, i_graph_nodes)
                self.graph_append(c.node4, i_graph_nodes)

            if c.type == "e":  # VVT/VCVS
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node3, v_graph_nodes)
                self.graph_append(c.node4, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self.graph_append(c.node3, i_graph_nodes)
                self.graph_append(c.node4, i_graph_nodes)

                self._collapse(i_graph_collapses, c.node1, c.node2)

            if c.type == "f":  # CCT/CCCS
                node1 = c.node1
                node2 = c.node2
                node3 = c.node3
                node4 = c.node4
                self.graph_append(node1, v_graph_nodes)
                self.graph_append(node2, v_graph_nodes)
                self.graph_append(node3, v_graph_nodes)
                self.graph_append(node4, v_graph_nodes)
                self.graph_append(node1, i_graph_nodes)
                self.graph_append(node2, i_graph_nodes)
                self.graph_append(node3, i_graph_nodes)
                self.graph_append(node4, i_graph_nodes)

                self._collapse(v_graph_collapses, node3, node4)
                matrix_col_expand += 1

            if c.type == "h":  # CVT/CCVS
                node1 = c.node1
                node2 = c.node2
                node3 = c.node3
                node4 = c.node4
                self.graph_append(node1, v_graph_nodes)
                self.graph_append(node2, v_graph_nodes)
                self.graph_append(node3, v_graph_nodes)
                self.graph_append(node4, v_graph_nodes)
                self.graph_append(node1, i_graph_nodes)
                self.graph_append(node2, i_graph_nodes)
                self.graph_append(node3, i_graph_nodes)
                self.graph_append(node4, i_graph_nodes)

                matrix_col_expand += 1
                self._collapse(i_graph_collapses, node1, node2)
                self._collapse(v_graph_collapses, node3, node4)

            if c.type == "a":
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node3, v_graph_nodes)
                self.graph_append(c.node4, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self.graph_append(c.node3, i_graph_nodes)
                self.graph_append(c.node4, i_graph_nodes)

                self._collapse(i_graph_collapses, c.node1, c.node2)
                self._collapse(v_graph_collapses, c.node3, c.node4)

            if c.type == "w": # short
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self._collapse(i_graph_collapses, c.node1, c.node2)
                self._collapse(v_graph_collapses, c.node1, c.node2)

        """Collapse nodes based on collapse dictionaries"""
        for collapse_list in i_graph_collapses:
            i = 0
            tmp_i_graph_nodes = copy.copy(i_graph_nodes)
            for n in tmp_i_graph_nodes:
                if n in collapse_list:
                    if "0" in collapse_list:
                        i_graph_nodes.remove(n)
                    else:
                        i_graph_nodes[i] = min(collapse_list)
                i += 1

        for collapse_list in v_graph_collapses:
            i = 0
            tmp_v_graph_nodes = copy.copy(v_graph_nodes)
            for n in tmp_v_graph_nodes:
                if n in collapse_list:
                    if "0" in collapse_list:
                        v_graph_nodes.remove(n)
                    else:
                        v_graph_nodes[i] = min(collapse_list)
                i += 1

        node_voltage_identities = v_graph_collapses
        v_graph_nodes = list(set(v_graph_nodes))
        i_graph_nodes = list(set(i_graph_nodes))

        rows = len(i_graph_nodes)
        cols = len(v_graph_nodes)

        m_size = len(v_graph_nodes) + matrix_col_expand

        M = sympy.Matrix(sympy.zeros(m_size))
        S = sympy.Matrix(sympy.zeros(m_size, 1))
        index_row = 0
        index_col = 0
        symbols_to_append = []
        for key in self.circuit.components:
            c = self.circuit.components[key]
            if c.type in ["r", "c"]:
                self._add_basic_tgn(M, S, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
            if c.type == "l":
                if c.coupling is None:
                    self._add_basic_tgn(M, S, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
            if c.type == "v":
                self._add_voltage_source_tgn(M, S, v_graph_nodes, i_graph_nodes, c, index_row, index_col,
                                             i_graph_collapses, v_graph_collapses)
                symbols_to_append.append(sympy.Symbol(f"i({c.name})"))
                index_row += 1
                index_col += 1
            if c.type == "i":
                self._add_current_source_tgn(M, S, v_graph_nodes, i_graph_nodes, c, i_graph_collapses)
            if c.type == "g":
                self._add_VCT_tgn(M, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
            if c.type == "e":
                self._add_VVT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_row, i_graph_collapses, v_graph_collapses)
                index_row += 1
            if c.type == "f":
                # raise NotImplementedError("CCCS not supported in 'two_graph_node' analysis")
                self._add_CCT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_col, i_graph_collapses)
                symbols_to_append.append(sympy.Symbol(f"i({c.name})"))
                index_col += 1
            if c.type == "h":
                self._add_CVT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_col, index_row, i_graph_collapses,
                                  v_graph_collapses)
                symbols_to_append.append(sympy.Symbol(f"i({c.name}_in)"))
                index_col += 1
                index_row += 1
            if c.type == "a":
                pass
            if c.type == "s":
                raise NotImplementedError("Switch not supported in 'two_graph_node' analysis")

        if len(self.circuit.couplings) > 0:
            for coupling in self.circuit.couplings:
                ind1 = self.circuit.components[coupling.L1]
                ind2 = self.circuit.components[coupling.L2]
                if self.is_symbolic:
                    coef = coupling.sym_value
                else:
                    coef = coupling.value
                self._add_K_tgn(M, S, v_graph_nodes, i_graph_nodes, [ind1, ind2, coef], index_col, index_row,
                                i_graph_collapses, v_graph_collapses)
                symbols_to_append.append(sympy.Symbol(f"i({ind1.name})"))
                symbols_to_append.append(sympy.Symbol(f"i({ind2.name})"))
                index_col += 2
                index_row += 2

        equation_matrix = M.col_insert(m_size, S)

        for node in v_graph_nodes:
            symbols.append(sympy.Symbol(f"v({node})"))
        for symb in symbols_to_append:
            symbols.append(symb)

        # TODO: experiment with simplification inside matrix - seems like a huge performance upgrade in tgn method!
        '''for i in range(m_size ** 2):
            expr = equation_matrix[i]
            if expr != 0:
                equation_matrix[i] = sympy.cancel(expr)'''

        return equation_matrix, symbols, node_voltage_identities

    @staticmethod
    def _collapse(graph_collapses, node1, node2):
        node1_in = None
        node2_in = None
        for i in range(len(graph_collapses)):
            if node1 in graph_collapses[i]:
                node1_in = i
            if node2 in graph_collapses[i]:
                node2_in = i
        if (node1_in is None) and (node2_in is None): # collapsed nodes not present in any existing collapse list
            graph_collapses.append([node1, node2])  # set node2 to be collapsed into node1 on the graph
        elif (node1_in is not None) and (node2_in is not None):
            graph_collapses[node2_in] = list(set(graph_collapses[node2_in]) | set(graph_collapses[node1_in]))
            del graph_collapses[node1_in]
        elif node1_in is not None:
            graph_collapses[node1_in].append(node2)  # set node2 to be collapsed into node1 on the graph
        elif node2_in is not None:
            graph_collapses[node2_in].append(node1)  # set node1 to be collapsed into node2 on the graph

    def _add_CVT_tgn(self, M, v_nodes, i_nodes, c, index_col, index_row, i_graph_collapses, v_graph_collapses):
        if self.is_symbolic:
            r = c.sym_value
        else:
            r = c.value

        node1 = c.node1
        node2 = c.node2
        node3 = c.node3
        node4 = c.node4

        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        n3i = self.index_tgn(i_nodes, node3, i_graph_collapses)
        n4i = self.index_tgn(i_nodes, node4, i_graph_collapses)

        col = len(v_nodes) + index_col
        row = len(i_nodes) + index_row

        if n3i is not None:
            M[n3i, col] += 1
        if n4i is not None:
            M[n4i, col] += -1
        if n1v is not None:
            M[row, n1v] += 1
        if n2v is not None:
            M[row, n2v] += -1

        M[row, col] += -r

    def _add_CCT_tgn(self, M, v_nodes, i_nodes, c, index, i_graph_collapses):
        if self.is_symbolic:
            f_gain = c.sym_value
        else:
            f_gain = c.value

        node1 = c.node1
        node2 = c.node2
        node3 = c.node3
        node4 = c.node4


        n1i = self.index_tgn(i_nodes, node3, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node4, i_graph_collapses)
        n3i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n4i = self.index_tgn(i_nodes, node2, i_graph_collapses)
        col = len(v_nodes) + index

        if n1i is not None:
            M[n1i, col] += 1/f_gain
        if n2i is not None:
            M[n2i, col] += -1/f_gain
        if n3i is not None:
            M[n3i, col] += 1
        if n4i is not None:
            M[n4i, col] += -1

    def _add_VCT_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses):
        if self.is_symbolic:
            g = c.sym_value
        else:
            g = c.value
        node1 = c.node1
        node2 = c.node2
        node3 = c.node3
        node4 = c.node4
        n1v = self.index_tgn(v_nodes, node3, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node4, v_graph_collapses)
        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
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

    def _add_VVT_tgn(self, M, v_nodes, i_nodes, c, index, i_graph_collapses, v_graph_collapses):
        if self.is_symbolic:
            e = c.sym_value
        else:
            e = c.value
        node1 = c.node1
        node2 = c.node2
        node3 = c.node3
        node4 = c.node4
        n1v = self.index_tgn(v_nodes, node3, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node4, v_graph_collapses)
        n3v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n4v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        row = len(i_nodes)+index
        if n1v is not None:
            M[row, n1v] += -e
        if n2v is not None:
            M[row, n2v] += e
        if n3v is not None:
            M[row, n3v] += 1
        if n4v is not None:
            M[row, n4v] += -1

    def _add_voltage_source_tgn(self, M, S, v_nodes, i_nodes, c, index_row, index_col, i_graph_collapses, v_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        if self.is_symbolic:
            form = "sym"
        else:
            form = "rat"
        val = c.v(form, self.analysis_type)

        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)

        col = len(v_nodes) + index_col
        row = len(i_nodes) + index_row

        if n1i is not None:
            M[n1i, col] += 1
        if n2i is not None:
            M[n2i, col] += -1
        if n1v is not None:
            M[row, n1v] += 1
        if n2v is not None:
            M[row, n2v] += -1
        S[row, 0] += val

    def _add_current_source_tgn(self, M, S, v_nodes, i_nodes, c, i_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        if self.is_symbolic:
            form = "sym"
        else:
            form = "rat"
        val = c.i(form, self.analysis_type)

        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)

        if n1i is not None:
            S[n1i, 0] += -val
        if n2i is not None:
            S[n2i, 0] += val

    def _add_basic_tgn(self, M, S, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses):
        node1 = c.node1
        node2 = c.node2

        y = c.y(self.is_symbolic)
        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
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

        # Only process ic if there is one
        if self.analysis_type == "tran":
            if c.type in ["c", "l"] and c.init_cond is not None:
                node1 = c.node1
                node2 = c.node2
                n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
                n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)

                if c.type == "c":
                    # Capacitor Norton current: C * V(0)
                    if n1i is not None:
                        S[n1i, 0] += c.val(self.is_symbolic) * c.init_cond
                    if n2i is not None:
                        S[n2i, 0] += -c.val(self.is_symbolic) * c.init_cond

                elif c.type == "l":
                    if c.coupling is None:  # If theres no coupling handle IC by Norton current IC model
                        # Inductor Norton current: I(0) / s
                        if n1i is not None:
                            S[n1i, 0] += -(1 / s) * c.init_cond
                        if n2i is not None:
                            S[n2i, 0] += (1 / s) * c.init_cond
                    else:  # Else let the coupling stamp handle initial conditions
                        pass

    def _add_K_tgn(self, M, S, v_nodes, i_nodes, coupling, index_col, index_row, i_graph_collapses, v_graph_collapses):
        """
        Adds a coupled inductor block to TGN matrix.
        """
        ind1 = coupling[0]
        ind2 = coupling[1]

        if self.is_symbolic:
            L1 = ind1.sym_value
            L2 = ind2.sym_value
        else:
            L1 = ind1.value
            L2 = ind2.value

        coeff = coupling[2] * sympy.sqrt(L1 * L2)

        node1 = ind1.node1
        node2 = ind1.node2
        node3 = ind2.node1
        node4 = ind2.node2

        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        n3v = self.index_tgn(v_nodes, node3, v_graph_collapses)
        n4v = self.index_tgn(v_nodes, node4, v_graph_collapses)

        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
        n3i = self.index_tgn(i_nodes, node3, i_graph_collapses)
        n4i = self.index_tgn(i_nodes, node4, i_graph_collapses)

        col1 = len(v_nodes) + index_col
        col2 = len(v_nodes) + index_col + 1
        row1 = len(i_nodes) + index_row
        row2 = len(i_nodes) + index_row + 1

        # L1 KVL row
        if n1v is not None: M[row1, n1v] += 1
        if n2v is not None: M[row1, n2v] += -1
        M[row1, col1] += -s * L1  # self-inductor
        M[row1, col2] += -s * coeff  # mutual term

        # L2 KVL row
        if n3v is not None: M[row2, n3v] += 1
        if n4v is not None: M[row2, n4v] += -1
        M[row2, col2] += -s * L2  # self-inductor
        M[row2, col1] += -s * coeff  # mutual term

        # Node current contributions (KCL)
        if n1i is not None: M[n1i, col1] += 1
        if n2i is not None: M[n2i, col1] += -1
        if n3i is not None: M[n3i, col2] += 1
        if n4i is not None: M[n4i, col2] += -1

        if self.analysis_type == "tran":
            ic1 = ind1.init_cond
            ic2 = ind2.init_cond

            if ic1 is not None:
                S[row2, 0] += -coeff * ic1 - L2 * ic2
            if ic2 is not None:
                S[row1, 0] += -coeff * ic2 - L1 * ic1

    @staticmethod
    def index_tgn(nodes, node, collapses):
        try:
            return nodes.index(node)
        except ValueError:
            for collapse_list in collapses:
                if node in collapse_list:
                    if "0" in collapse_list:
                        return None
                    else:
                        return nodes.index(min(collapse_list))
                elif node == "0":
                    return None
    @staticmethod
    def graph_append(node, graph):
        if node == '0':
            pass
        elif node not in graph:
            graph.append(node)
        return graph


class TableauBuilder:
    def __init__(self, circuit, analysis_type, is_symbolic):
        self.analysis_type = analysis_type
        self.is_symbolic = is_symbolic
        self.circuit = circuit

    def build(self):
        port_count = self.circuit.count_ports()
        node_count = self.circuit.count_nodes()
        node_symbols = self.circuit.get_node_symbols()
        node_dict = self.circuit.get_node_dict()

        M = sympy.Matrix(sympy.zeros(2 * port_count + node_count))
        R = sympy.Matrix(sympy.zeros(port_count * 2 + node_count, 1))
        for i in range(port_count):
            M[i, i] = 1
        index = 0
        voltage_symbols = []
        current_symbols = []
        inductor_index = {}
        couplings = []
        for key in self.circuit.components:
            if self.circuit.components[key].type == "k":
                couplings.append(self.circuit.components[key])
            else:
                c = self.circuit.components[key]
                if c.type in ["r", "l", "c"]:
                    self._add_basic(M, R, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    if c.type == "l":
                        inductor_index[c.name] = index

                if c.type == "v":
                    self._add_voltage_source(M, R, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                if c.type == "i":
                    self._add_current_source(M, R, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                if c.type == "g":
                    self._add_VCT(M, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "e":
                    self._add_VVT(M, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "f":
                    self._add_CCT(M, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "h":
                    self._add_CVT(M, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "a":
                    self._add_A(M, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "s":
                    self._add_short(M, c, index, port_count, node_dict)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))

                index += 1
        for coupling in couplings:
            self._add_K(M, R, coupling, inductor_index, port_count)

        equation_matrix = M.col_insert(port_count * 2 + node_count, R)
        symbols = voltage_symbols + current_symbols + node_symbols
        return equation_matrix, symbols

    @staticmethod
    def _incidence_matrix_write(N1, N2, matrix, index, port_count, node_dict):
        if N1 == "0":
            pass
        else:
            node_pos = node_dict[N1]
            matrix[port_count*2 + node_pos, port_count + index] += 1
            matrix[index, port_count*2 + node_pos] -= 1
        if N2 == "0":
            pass
        else:
            node_pos = node_dict[N2]
            matrix[port_count*2 + node_pos, port_count + index] += -1
            matrix[index, port_count*2 + node_pos] -= -1
        return matrix

    def _add_basic(self, matrix, vi_vector, c, index, port_count, node_dict):
        N1 = c.node1
        N2 = c.node2
        y_b = 0
        z_b = 0

        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
            if val is None:
                val = c.sym_value

        if c.type == "r":
            y_b = 1
            z_b = -val
        elif c.type == "l":
            y_b = 1
            z_b = -s * val
            if self.analysis_type == "tran":
                vi_vector[port_count+index, 0] += -val * c.init_cond
        elif c.type == "c":
            y_b = s * val
            z_b = -1
            if self.analysis_type == "tran":
                vi_vector[port_count+index, 0] += val * c.init_cond

        matrix[port_count+index, index] += y_b
        matrix[port_count+index, port_count+index] += z_b

        self._incidence_matrix_write(N1, N2, matrix, index, port_count, node_dict)
        return matrix

    def _add_voltage_source(self, matrix, result, c, index, port_count, node_dict):
        N1 = c.node1
        N2 = c.node2
        if self.is_symbolic:
            form = "sym"
        else:
            form = "rat"
        val = c.val(form, self.analysis_type)

        matrix[port_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index, port_count, node_dict)
        result[port_count + index, 0] = val
        return matrix

    def _add_current_source(self, matrix, result, c, index, port_count, node_dict):
        N1 = c.node1
        N2 = c.node2
        if self.is_symbolic:
            form = "sym"
        else:
            form = "rat"
        val = c.val(form, self.analysis_type)
        matrix[port_count + index, port_count + index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index, port_count, node_dict)
        result[port_count + index, 0] = val
        return matrix

    def _add_VCT(self, matrix, c, index, port_count, node_dict):  # voltage to current transformer
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        matrix[port_count + index, port_count + index] = 1
        matrix[port_count + index + 1, index] = val
        matrix[port_count + index + 1, port_count + index + 1] = -1
        self._incidence_matrix_write(N3, N4, matrix, index, port_count, node_dict)
        self._incidence_matrix_write(N1, N2, matrix, index + 1, port_count, node_dict)
        return matrix

    def _add_VVT(self, matrix, c, index, port_count, node_dict):  # voltage to voltage transformer
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        matrix[port_count + index, port_count + index] = 1
        matrix[port_count + index + 1, index] = val
        matrix[port_count + index + 1, index + 1] = -1
        self._incidence_matrix_write(N3, N4, matrix, index, port_count, node_dict)
        self._incidence_matrix_write(N1, N2, matrix, index + 1, port_count, node_dict)
        return matrix

    def _add_CCT(self, matrix, c, index, port_count, node_dict):  # current to current transformer
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        matrix[port_count + index, index] = 1
        matrix[port_count + index + 1, port_count + index] = val
        matrix[port_count + index + 1, port_count + index + 1] = -1
        self._incidence_matrix_write(N3, N4, matrix, index, port_count, node_dict)
        self._incidence_matrix_write(N1, N2, matrix, index + 1, port_count, node_dict)
        return matrix

    def _add_CVT(self, matrix, c, index, port_count, node_dict):  # Current to voltage transformer
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4

        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        matrix[port_count + index, index] = 1
        matrix[port_count + index + 1, port_count + index] = val
        matrix[port_count + index + 1, index + 1] = -1
        self._incidence_matrix_write(N3, N4, matrix, index, port_count, node_dict)
        self._incidence_matrix_write(N1, N2, matrix, index + 1, port_count, node_dict)
        return matrix

    def _add_A(self, matrix, c, index, port_count, node_dict):  # Ideal OAMP
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        matrix[port_count + index, index] = 1
        matrix[port_count + index + 1, port_count + index] = 1
        self._incidence_matrix_write(N3, N4, matrix, index, port_count, node_dict)
        self._incidence_matrix_write(N1, N2, matrix, index + 1, port_count, node_dict)
        return matrix

    def _add_K(self, matrix, vi_vector, c, inductor_index, port_count):  # coupled inductors
        c_L1 = self.circuit.components[c.L1]
        c_L2 = self.circuit.components[c.L2]
        L1_index = inductor_index[c_L1.name]
        L2_index = inductor_index[c_L2.name]
        if self.is_symbolic:
            L1 = c_L1.sym_value
            L2 = c_L2.sym_value
            coupling_coeff = c.sym_value
            coeff = coupling_coeff * sympy.sqrt(L1 * L2)
        else:
            L1 = c_L1.value
            L2 = c_L2.value
            coupling_coeff = c.value
            coeff = coupling_coeff * sympy.sqrt(L1 * L2)

        matrix[port_count + L2_index, port_count + L1_index] += -s * coeff
        matrix[port_count + L1_index, port_count + L2_index] += -s * coeff

        if self.analysis_type == "tran":
            if c_L1.init_cond is not None:
                vi_vector[port_count + L2_index, 0] += -coeff * c_L1.init_cond
            if c_L2.init_cond is not None:
                vi_vector[port_count + L1_index, 0] += -coeff * c_L2.init_cond

        return coeff, c_L1, c_L2, L1_index, L2_index

    def _add_short(self, matrix, c, index, port_count, node_dict):
        N1 = c.node1
        N2 = c.node2
        matrix[port_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index, port_count, node_dict)
        return matrix