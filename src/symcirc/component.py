from typing import List, Set, Dict, Union, Any
import sympy

from symcirc.utils import s, j, pi
from symcirc.graph_signalflow import Branch
from symcirc.parse_utils import dc_value, ac_value, tran_value, value_enum


ANALYSIS_DEPENDENT = ["i", "v"]


class Component:
    """
    Parent component class

    :param name: component id
    :param node1: first node id
    :param node2: second node id
    :param sym_value: symbolic value
    :param value: numeric value
    """
    def __init__(self, name:str, node1:str, node2:str,
                 sym_value:Union[sympy.Symbol,sympy.Expr]=None, value:sympy.Rational=None, value_float:float=None):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.type = None

        self.sym_value = sym_value
        self.value = value
        self.value_float = value_float

    def nodes(self) -> Set[str]:
        """
        :return: set of component nodes
        """
        ret = set()
        if self.node1 is not None:
            ret.add(self.node1)
        if self.node2 is not None:
            ret.add(self.node2)
        return ret

    def set_value(self, value):
        rat, flt, sym, is_symbolic = value_enum(self.name, value)
        self.value = rat
        self.value_float = flt
        self.sym_value = sym

    def short(self):
        return Short(f"Short_{self.name}", self.node1, self.node2)



class ParallelAdmittance(Component):
    """
        Resistor component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param sym_value: symbolic value
        :param value: numeric value used in semisymbolic analyes
    """
    def __init__(self, component):

        name: str = component.name
        node1: str = component.node1
        node2: str = component.node2
        sym_value = component.y(symbolic=True)
        value = component.y(symbolic=False)
        value_float = component.y(symbolic=False, is_float=True)

        super().__init__(name, node1, node2, sym_value, value, value_float)

        self.netlist_keywords = ["R", "r"] # TODO: change to Z
        self.type = "r"
        self.subcomponents = {component.name: component}

    def add(self, component):
        self.name += f"|{component.name}"
        self.sym_value += component.y(symbolic=True)
        self.value += component.y(symbolic=False)
        self.value_float += component.y(symbolic=False, is_float=True)
        self.subcomponents[component.name] = component

    def remove(self, component):
        self.name = ""
        self.sym_value = 0
        self.value = 0
        self.value_float = 0
        self.subcomponents.pop(component.name)
        for component in self.subcomponents.values():
            if self.name != "":
                self.name += f"|{component.name}"
            else:
                self.name = component.name
            self.sym_value += component.y(symbolic=True)
            self.value += component.y(symbolic=False)
            self.value_float += component.y(symbolic=False, is_float=True)

    def z(self, symbolic=False, is_float=False):
        if symbolic:
            return 1 / self.sym_value
        elif is_float:
            return 1/ self.value_float
        else:
            return 1 / self.value

    def y(self, symbolic=False, is_float=False):
        if symbolic:
            return self.sym_value
        elif is_float:
            return self.value_float
        else:
            return self.value

    def branches(self):
        branches = []
        vnode1 = f"{self.node1}v"
        vnode2 = f"{self.node2}v"
        inode1 = f"{self.node1}i"
        inode2 = f"{self.node2}i"

        sym_adm = self.y(symbolic=True)
        num_adm = self.y(symbolic=False)
        sym_imp = self.z(symbolic=True)
        num_imp = self.z(symbolic=False)

        branches.append(Branch(self.name, inode1, vnode1, sym_imp, num_imp))
        branches.append(Branch(self.name, inode2, vnode2, sym_imp, num_imp))
        branches.append(Branch(self.name, vnode1, inode2, sym_adm, num_adm))
        branches.append(Branch(self.name, vnode2, inode1, sym_adm, num_adm))
        return branches

class SerialAdmittance(ParallelAdmittance):
    def add(self, component):
        self.name = f"({self.name})({component.name})"
        self.sym_value = 1/(self.z(symbolic=True) + component.z(symbolic=True))
        self.value = 1/(self.z(symbolic=False) + component.z(symbolic=False))
        self.value_float = 1/(self.z(symbolic=False, is_float=True) + component.z(symbolic=False, is_float=True))
        self.subcomponents[component.name] = component
        old_n_self = self.nodes()
        old_n_comp = component.nodes()
        all_nodes = set()
        all_nodes.update(old_n_self)
        all_nodes.update(old_n_comp)
        nodes = list(all_nodes - (old_n_comp & old_n_self))
        self.node1 = nodes[0]
        self.node2 = nodes[1]

    def remove(self, component):
        self.name = ""
        self.sym_value = 0
        self.value = 0
        self.value_float = 0
        self.subcomponents.pop(component.name)
        for component in self.subcomponents.values():
            if self.name != "":
                self.name += f"({component.name})"
            else:
                self.name = component.name
            self.sym_value = 1 / (self.z(symbolic=True) + component.z(symbolic=True))
            self.value = 1 / (self.z(symbolic=False) + component.z(symbolic=False))
            self.value_float = 1 / (self.z(symbolic=False, is_float=True) + component.z(symbolic=False, is_float=True))


class Port:
    def __init__(self, name: str, node1: str, node2: str, component: Component):
        self.name = name
        self.component = component
        self.node1 = node1
        self.node2 = node2
        self.type = "port"

    def nodes(self) -> Set[str]:
        """
        :return: set of component nodes
        """
        ret = set()
        if self.node1 is not None:
            ret.add(self.node1)
        if self.node2 is not None:
            ret.add(self.node2)
        return ret


class Coupling(Component):
    def __init__(self, name, L1, L2, sym_value, value, value_float):
        super().__init__(name, None, None, sym_value, value, value_float)
        self.L1 = L1
        self.L2 = L2
        self.sym_value = sym_value
        self.value = value
        self.value_float = value_float
        self.netlist_keywords = ["K", "k"]
        self.type = "k"

class Resistor(Component):
    """
        Resistor component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param sym_value: symbolic value
        :param value: numeric value used in semisymbolic analyes
    """
    def __init__(self, name:str, node1:str, node2:str, sym_value:Union[sympy.Symbol, sympy.Expr]=None,
                 value: sympy.Rational=None, value_float: float=None):
        if sym_value is None:
            sym_value = sympy.Symbol(name)
        super().__init__(name, node1, node2, sym_value, value, value_float)
        self.netlist_keywords = ["R", "r"]
        self.type = "r"

    def y(self, symbolic=False, is_float=False):
        if symbolic:
            return 1 / self.sym_value
        elif is_float:
            return 1 / self.value_float
        else:
            return 1 / self.value

    def z(self, symbolic=False, is_float=False):
        if symbolic:
            return self.sym_value
        elif is_float:
            return self.value_float
        else:
            return self.value

    def branches(self):
        branches = []
        vnode1 = f"{self.node1}v"
        vnode2 = f"{self.node2}v"
        inode1 = f"{self.node1}i"
        inode2 = f"{self.node2}i"

        sym_adm = self.y(symbolic=True)
        num_adm = self.y(symbolic=False)
        sym_imp = self.z(symbolic=True)
        num_imp = self.z(symbolic=False)

        branches.append(Branch(self.name, inode1, vnode1, sym_imp, num_imp))
        branches.append(Branch(self.name, inode2, vnode2, sym_imp, num_imp))
        branches.append(Branch(self.name, vnode1, inode2, sym_adm, num_adm))
        branches.append(Branch(self.name, vnode2, inode1, sym_adm, num_adm))
        return branches


class Capacitor(Component):
    """
        Capacitor component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param sym_value: symbolic value
        :param value: numeric value used in semisymbolic analyes
        :param init_cond: initial condition in volts
    """
    def __init__(self, name:str, node1:str, node2:str, sym_value:Union[sympy.Symbol, sympy.Expr]=None,
                value:sympy.Rational=None, value_float: float=None, init_cond:Union[sympy.Symbol, sympy.Expr]=0):
        if sym_value is None:
            sym_value = sympy.Symbol(name)
        super().__init__(name, node1, node2, sym_value, value, value_float)
        self.init_cond = init_cond
        self.netlist_keywords = ["C", "c"]
        self.type = "c"

    def y(self, symbolic=False, is_float=False):
        if symbolic:
            return s * self.sym_value
        elif is_float:
            return s * self.value_float
        else:
            return s * self.value

    def z(self, symbolic=False, is_float=False):
        if symbolic:
            return 1/ (s * self.sym_value)
        elif is_float:
            return 1/ (s * self.value_float)
        else:
            return 1/ (s * self.value)

    def branches(self):
        branches = []
        vnode1 = f"{self.node1}v"
        vnode2 = f"{self.node2}v"
        inode1 = f"{self.node1}i"
        inode2 = f"{self.node2}i"

        sym_adm = self.y(symbolic=True)
        num_adm = self.y(symbolic=False)
        sym_imp = self.z(symbolic=True)
        num_imp = self.z(symbolic=False)

        branches.append(Branch(self.name, inode1, vnode1, sym_imp, num_imp))
        branches.append(Branch(self.name, inode2, vnode2, sym_imp, num_imp))
        branches.append(Branch(self.name, vnode1, inode2, sym_adm, num_adm))
        branches.append(Branch(self.name, vnode2, inode1, sym_adm, num_adm))
        return branches


class Inductor(Component):
    """
        Inductor component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param sym_value: symbolic value
        :param value: numeric value used in semisymbolic analyes
        :param init_cond: initial condition in amps
        :param coupling: a coupling object defining the couple
    """
    def __init__(self, name:str, node1:str, node2:str, sym_value:Union[sympy.Symbol, sympy.Expr]=None,
                 value:sympy.Rational=None, value_float: float=None, init_cond:Union[sympy.Symbol, sympy.Expr]=0,
                 coupling:Union[Coupling, None]=None):
        if sym_value is None:
            sym_value = sympy.Symbol(name)
        super().__init__(name, node1, node2, sym_value, value, value_float)
        self.init_cond = init_cond
        self.coupling = coupling
        self.netlist_keywords = ["L", "l"]
        self.type = "l"

    def y(self, symbolic=False, is_float=False):
        if symbolic:
            return 1 / (self.sym_value * s)
        elif is_float:
            return 1 / (self.value_float * s)
        else:
            return 1 / (self.value * s)

    def z(self, symbolic=False, is_float=False):
        if symbolic:
            return self.sym_value * s
        elif is_float:
            return self.value_float * s
        else:
            return self.value * s

    def branches(self):
        branches = []
        vnode1 = f"{self.node1}v"
        vnode2 = f"{self.node2}v"
        inode1 = f"{self.node1}i"
        inode2 = f"{self.node2}i"

        sym_adm = self.y(symbolic=True)
        num_adm = self.y(symbolic=False)
        sym_imp = self.z(symbolic=True)
        num_imp = self.z(symbolic=False)

        branches.append(Branch(self.name, inode1, vnode1, sym_imp, num_imp))
        branches.append(Branch(self.name, inode2, vnode2, sym_imp, num_imp))
        branches.append(Branch(self.name, vnode1, inode2, sym_adm, num_adm))
        branches.append(Branch(self.name, vnode2, inode1, sym_adm, num_adm))
        return branches

class IndependentSource(Component):
    """
            Independent source component base class

            :param name: component id
            :param node1: first node id
            :param node2: second node id
            :param dc_num: numeric dc value used in semisymbolic dc analyes
            :param dc_sym: symbolic dc value
            :param ac_num: numeric ac value used in semisymbolic ac analysis
            :param ac_sym: symbolic ac value
            :param ac_phase: ac phase value
            :param tran_num: numeric transient value used in semisymbolic transient analysis
        """

    def __init__(self, name: str, node1: str, node2: str,
                 values: Dict[str, Any]):

        self.name = name
        self._dc_raw = values["dc_raw"]
        self._ac_raw = values["ac_raw"]
        self._tran_raw = values["tran_raw"]

        self._dc_num, self._dc_float, self._dc_sym = None, None, None
        self._ac_num, self._ac_float, self._ac_sym, self._ac_phase = None, None, None, None
        self._tran_num, self._tran_float = None, None

        self.set_dc_val(self._dc_raw)
        self.set_ac_val(self._ac_raw)
        self.set_tran_val(self._tran_raw)

        super().__init__(name, node1, node2, sym_value=self._dc_sym, value=self._dc_num)

    def set_dc_val(self, val):
        self._dc_raw = val
        self._dc_num, self._dc_float, self._dc_sym = dc_value(self.name, val)
        self._dc_raw = val
    def set_ac_val(self, val):
        self._ac_raw = val
        self._ac_num, self._ac_float, self._ac_sym, self._ac_phase = ac_value(self.name, val)
    def set_tran_val(self, val):
        self._tran_raw = val
        self._tran_num, self._tran_float = tran_value(self.name, val)

    def val(self, form="sym", analysis_type="dc"):
        # TODO: document options
        analysis_type = analysis_type.lower()
        if analysis_type == 'ac':
            return self.ac_val(form)
        elif analysis_type == 'dc':
            return self.dc_val(form)
        elif analysis_type == 'tf':
            return self.tf_val(form)
        elif analysis_type == 'tran':
            return self.tran_val(form)

    def dc_val(self, form="sym"):
        if form == "sym":
            mag = self._dc_sym
        elif form == "float":
            mag = self._dc_float
        else:
            mag = self._dc_num
        return mag

    def ac_val(self, form="sym"):
        has_phase = self._ac_phase != 0
        if form == "sym":
            if has_phase:
                mag = self._ac_num
            else:
                mag = self._ac_sym
        elif form == "float":
            mag = self._ac_float
        else:
            mag = self._ac_num

        if has_phase:
            ret = mag * sympy.exp(j * self._ac_phase)
        else:
            ret = mag
        return ret

    def tf_val(self, form="sym"):
        if form == "sym":
            mag = self._ac_sym
        elif form == "float":
            mag = self._ac_float
        else:
            mag = self._ac_num
        return mag

    def tran_val(self, form="sym"):
        if form == "float":
            params = self._tran_float
        else:
            params = self._tran_num

        if params[0] == "sin":
            offset = params[1]
            amp = params[2]
            freq = params[3]
            delay = params[4]
            damping = params[5]
            omega = 2 * pi * freq

            tran = amp * ((damping + s) * sympy.sin(delay) + omega * sympy.cos(delay)) / (
                    damping ** 2 + 2 * damping * s + omega ** 2 + s ** 2)
        else:
            tran = self.dc_val(form="num") / s
        return tran


class VoltageSource(IndependentSource):
    """
        Independent voltage source component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param values: a list of values defining the source behaviour
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.netlist_keywords = ["V", "v", "U", "u"]
        self.type = "v"

    def v(self, form="sym", analysis_type="dc"):
        return self.val(form, analysis_type)


class CurrentSource(IndependentSource):
    """
        Independent current source component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param values: a list of values defining the source behaviour
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.netlist_keywords = ["I", "i"]
        self.type = "i"

    def i(self, form="sym", analysis_type="dc"):
        return self.val(form, analysis_type)


class IdealOperationalAmplifier(Component):
    """
        Operational amplifier component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param node3: third node id
        :param node4: fourth node id
    """
    def __init__(self, name:str, node1:str, node2:str, node3:str, node4:str):
        super().__init__(name, node1, node2)
        self.node3 = node3
        self.node4 = node4
        self.netlist_keywords = ["A", "a"]
        self.type = "a"

    def nodes(self) -> Set[str]:
        """
        :return: set of component nodes
        """
        return {self.node1, self.node2, self.node3, self.node4}

    def generate_ports(self):
        self.port1 = Port(f"{self.name}_out", self.node1, self.node2, self)
        self.port2 = Port(f"{self.name}_in", self.node3, self.node4, self)

    def ports(self) -> List[Port]:
        return [self.port1, self.port2]


class CurrentControlledSource(Component):
    """
        Current controlled source component class

        :param name: component id
        :param type: component type id
        :param node1: first node id
        :param node2: second node id
        :param current_sensor: id of the element across which is the controlling current
        :param sym_value: symbolic value of the component
        :param value: numeric value of the component
    """
    def __init__(self, name:str, type:str, node1:str, node2:str, node3:str, node4:str,
                 sym_value:Union[sympy.Symbol, sympy.Expr], value: sympy.Rational, value_float: float):
        super().__init__(name, node1, node2, sym_value, value, value_float)
        self.node3 = node3
        self.node4 = node4
        self.type = type
        self.netlist_keywords = ["F", "f", "H", "h"]

    def nodes(self):
        return {self.node1, self.node2, self.node3, self.node4}

    def generate_ports(self):
        self.port1 = Port(f"{self.name}_out", self.node1, self.node2, self)
        self.port2 = Port(f"{self.name}_ctrl", self.node3, self.node4, self)

    def ports(self) -> List[Port]:
        return [self.port1, self.port2]


class VoltageControlledSource(Component):
    """
        Voltage controlled source component class

        :param name: component id
        :param type: component type id
        :param node1: first node id
        :param node2: second node id
        :param node3: third node id
        :param node4: fourth node id
        :param sym_value: symbolic value of the component
        :param value: numeric value of the component
    """
    def __init__(self, name:str, type:str, node1:str, node2:str, node3:str, node4:str,
                 sym_value:Union[sympy.Symbol, sympy.Expr], value:sympy.Rational, value_float: float):
        super().__init__(name, node1, node2, sym_value, value, value_float)
        self.node3 = node3
        self.node4 = node4
        self.type = type
        self.netlist_keywords = ["G", "g", "E", "e"]

    def nodes(self):
        return [self.node1, self.node2, self.node3, self.node4]

    def generate_ports(self):
        self.port1 = Port(f"{self.name}_out", self.node1, self.node2, self)
        self.port2 = Port(f"{self.name}_ctrl", self.node3, self.node4, self)

    def ports(self) -> List[Port]:
        return [self.port1, self.port2]

    def gain(self, symbolic=False, is_float=False):
        if symbolic:
            return self.sym_value
        if is_float:
            return self.value_float
        else:
            return self.value

    def branches(self):
        if self.type == "g":
            branches = []
            vnode3 = f"{self.node3}v"
            vnode4 = f"{self.node4}v"
            inode1 = f"{self.node1}i"
            inode2 = f"{self.node2}i"

            sym_gain = self.gain(symbolic=True)
            num_gain = self.gain(symbolic=False)
            if self.node3 != self.node1:
                branches.append(Branch(self.name, vnode3, inode1, -sym_gain, -num_gain))
            if self.node4 != self.node2:
                branches.append(Branch(self.name, vnode4, inode2, -sym_gain, -num_gain))
            if self.node3 != self.node2:
                branches.append(Branch(self.name, vnode3, inode2, sym_gain, num_gain))
            if self.node4 != self.node1:
                branches.append(Branch(self.name, vnode4, inode1, sym_gain, num_gain))
            return branches


class SubcktModel:
    """
        Voltage controlled source component class

        :param model_id: component id
        :param node_list: subcircuit interface node list
        :param param_dict: subcircuit parameters dictionary
    """
    def __init__(self, model_id:str, node_list:List[str], param_dict:Dict[str, Union[sympy.Symbol, sympy.Expr]]):
        self.model_id = model_id
        self.node_list = node_list
        self.param_dict = self._sanitize_param_dict(param_dict)
        self.elements = []

    def build_instance(self, operational_point) -> List[str]:
        return self.elements

    @staticmethod
    def _sanitize_param_dict(param_dict):
        sanitized_param_dict = {}
        for key in param_dict:
            sanitized_param_dict[key.lower()] = param_dict[key]
        param_dict = sanitized_param_dict
        return param_dict


class DiodeModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["a", "k"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op: Dict[str, Any] = None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        if "gd" in params:
            elements.append(f"rd a k 1/{param_dict['gd']}")
        else:
            elements.append("rd a k")

        return elements


class NPNModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["c", "b", "e"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op:Dict[str, Any]=None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        # Resistors
        if "gpi" in params:
            elements.append(f"rpi b e 1/{param_dict['gpi']}")
        else:
            elements.append(f"rpi b e")
        if "go" in params:
            elements.append(f"ro c e 1/{param_dict['go']}")
        elif "vaf" in params:
            elements.append(f"ro c e")

        #if "gx" in params:
        #    elements.append(f"rx b e 1/{param_dict['gx']}")

        # VCCS
        if "gm" in params:
            elements.append(f"gm c e b e {param_dict['gm']}")
        else:
            elements.append(f"gm c e b e")

        # Parasitic caps
        if "cmu" in params:
            elements.append(f"cmu b c {param_dict['cmu']}")
        elif "cje" in params:
            elements.append(f"cmu b c {param_dict['go']}*{param_dict['tr']}+{param_dict['cjc']}")

        if "cpi" in params:
            elements.append(f"cpi b c {param_dict['cpi']}")
        elif "cjc" in params:
            elements.append(f"cpi b e {param_dict['gm']}*{param_dict['tf']}+2*{param_dict['cje']}")

        return elements


class PNPModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["c", "b", "e"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op: Dict[str, Any] = None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        # Resistors
        if "gpi" in params:
            elements.append(f"rpi b e 1/{param_dict['gpi']}")
        else:
            elements.append(f"rpi b e")
        if "go" in params:
            elements.append(f"ro c e 1/{param_dict['go']}")
        elif "vaf" in params:
            elements.append(f"ro c e")
        if "gx" in params:
            elements.append(f"rx b e 1/{param_dict['gx']}")

        # VCCS
        if "gm" in params:
            elements.append(f"gm c e e b {param_dict['gm']}")
        else:
            elements.append(f"gm c e e b")


        # Parasitic caps
        if "cmu" in params:
            elements.append(f"cmu b c {param_dict['cmu']}")
        elif "cje" in params:
            elements.append(f"cmu b c {param_dict['go']}*{param_dict['tr']}+{param_dict['cjc']}")

        if "cpi" in params:
            elements.append(f"cpi b c {param_dict['cpi']}")
        elif "cjc" in params:
            elements.append(f"cpi b e {param_dict['gm']}*{param_dict['tf']}+2*{param_dict['cje']}")

        return elements


class NFETModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["d", "g", "s", "b"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op: Dict[str, Any] = None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        # VCCS
        if "gm" in params:
            elements.append(f"gm d s g s {param_dict['gm']}")
        else:
            elements.append("gm d s g s")

        if "gmb" in params:  # body-effect VCCS
            elements.append(f"gmb d s b s {param_dict['gmb']}")

        # Parasitic caps
        if "cbd" in params:  # parasitic cap bulk-drain
            elements.append(f"cbd b d {param_dict['cbd']}")
        if "cbs" in params:  # parasitic cap bulk-source
            elements.append(f"cbs b s {param_dict['cbs']}")
        if "cgb" in params:  # parasitic cap gate-bulk
            elements.append(f"cgb g b {param_dict['cgb']}")
        if "cgs" in params:  # parasitic cap gate-source
            elements.append(f"cgs g s {param_dict['cgs']}")
        if "cgd" in params:  # parasitic cap gate-drain
            elements.append(f"cgd g d {param_dict['cgd']}")
        if "cds" in params:  # parasitic cap drain-source
            elements.append(f"cds d s {param_dict['cds']}")

        # Parasitic res
        if "gds" in params:  # parasitic res drain-source
            elements.append(f"rds d s 1/{param_dict['gds']}")
        if "gbs" in params:  # parasitic res bulk-source
            elements.append(f"rbs b s 1/{param_dict['gbs']}")
        if "gbd" in params:  # parasitic res bulk-drain
            elements.append(f"rbd b d 1/{param_dict['gbd']}")

        return elements



class PFETModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["d", "g", "s", "b"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op: Dict[str, Any] = None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        # VCCS
        if "gm" in params:  # main VCCS
            elements.append(f"gm d s g s {param_dict['gm']}")
        else:
            elements.append("gm d s g s")

        if "gmb" in params: # body-effect VCCS
            elements.append(f"gmb d s b s {param_dict['gmb']}")

        # Parasitic caps
        if "cbd" in params:  # parasitic cap bulk-drain
            elements.append(f"cbd b d {param_dict['cbd']}")
        if "cbs" in params:  # parasitic cap bulk-source
            elements.append(f"cbs b s {param_dict['cbs']}")
        if "cgb" in params:  # parasitic cap gate-bulk
            elements.append(f"cgb g b {param_dict['cgb']}")
        if "cgs" in params:  # parasitic cap gate-source
            elements.append(f"cgs g s {param_dict['cgs']}")
        if "cgd" in params:  # parasitic cap gate-drain
            elements.append(f"cgd g d {param_dict['cgd']}")
        if "cds" in params:  # parasitic cap drain-source
            elements.append(f"cds d s {param_dict['cds']}")

        # Parasitic res
        if "gds" in params:  # parasitic res drain-source
            elements.append(f"rds d s 1/{param_dict['gds']}")
        if "gbs" in params:  # parasitic res bulk-source
            elements.append(f"rbs b s 1/{param_dict['gbs']}")
        if "gbd" in params:  # parasitic res bulk-drain
            elements.append(f"rbd b d 1/{param_dict['gbd']}")

        return elements


class PeriodicSwitch(Component):
    # periodic switch used for SC/SI analysis
    def __init__(self, name, node1, node2, phase):
        super().__init__(name, node1, node2)
        self.phase = phase
        self.netlist_keywords = ["S", "s"]


class Short(Component):
    def __init__(self, name, node1, node2):
        super().__init__(name, node1, node2)
        self.netlist_keywords = ["W", "w"]
        self.type = "w"
