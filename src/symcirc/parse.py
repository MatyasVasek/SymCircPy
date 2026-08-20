from symcirc.component import *
from symcirc.parse_utils import *

def preparse_source(words):
    """
    Locates the DC, AC, and TRAN signatures within a source (V/I) netlist
    line's word list

    :param words: the full split netlist line
    """
    lower_words = [w.lower() for w in words]
    dc_sig = None
    if len(words) < 3:
        dc_sig = None
    elif "dc" in lower_words:
        dc_index = lower_words.index("dc")
        dc_sig = words[dc_index+1]
    elif len(words) == 4:
        dc_sig = words[3]

    ac_sig = None
    try:
        if "ac" in lower_words:
            ac_index = lower_words.index("ac")
            mag = words[ac_index+1]
            try:
                phase = words[ac_index+2] if words[ac_index+2] not in RESERVED else None
            except IndexError:
                phase = None
            ac_sig = [mag, phase]
    except IndexError:
        ac_sig = None

    tran_sig = None
    if "sin" in lower_words:
        index = lower_words.index("sin")
        tran_sig = words[index:index + 6]

    return dc_sig, ac_sig, tran_sig

def parse_subcircuits(netlist, operating_point):
    subckt_models = {}
    in_model = False
    subckt_model_id = ""
    current_model = None
    parsed_netlist = []

    for line in netlist[1:]:
        words = line.split()
        keyword = words[0].lower()
        if words in [[], "\n", " "]:
            pass

        elif words[0][0] not in NETLIST_KEYCHARS:
            raise NotImplementedError(f"Keyword/Element '{words[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")

        elif keyword == ".subckt":
            in_model = True
            loading_nodes = True
            node_list = []
            param_dict = {}
            subckt_model_id = words[1]
            for w in words[2:]:
                if w in ["PARAMS:", "params:"]:
                    loading_nodes = False
                elif loading_nodes:
                    node_list.append(w)
                else:
                    key, val = w.split("=")
                    param_dict[key], _ = convert_units(val)
            current_model = SubcktModel(subckt_model_id, node_list, param_dict)

        elif keyword == ".model":
            model_id = words[1]
            model_type = words[2].lower()
            # strip optional spice model brackets
            words[3] = words[3].replace("(", "")
            words[-1] = words[-1].replace(")", "")
            # parse arguments

            param_dict = {"nr": 1}
            for w in words[3:]:
                key, val = w.split("=")
                param_dict[key], _ = convert_units(val)
            if model_type in ["npn"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"NPN model not implemented outside of AC or TF analysis.")
                model = NPNModelAC(model_id, param_dict)
                subckt_models[model_id] = model
            elif model_type in ["pnp"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"PNP model not implemented outside of AC or TF analysis.")
                model = PNPModelAC(model_id, param_dict)
                subckt_models[model_id] = model
            elif model_type in ["nmos", "njf"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"NFET model not implemented outside of AC or TF analysis.")
                model = NFETModelAC(model_id, param_dict)
                subckt_models[model_id] = model
            elif model_type in ["pmos", "pjf"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"PFET model not implemented outside of AC or TF analysis.")
                model = PFETModelAC(model_id, param_dict)
                subckt_models[model_id] = model
            elif model_type in ["d"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"Diode model not implemented outside of AC or TF analysis.")
                model = DiodeModelAC(model_id, param_dict)
                subckt_models[model_id] = model

        elif words[0][0] == ".":
            if keyword == ".ends":
                in_model = False
                subckt_models[subckt_model_id] = current_model
            elif keyword == ".end":
                break
            else:
                raise NotImplementedError(f"Keyword/Element '{words[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")

        elif in_model:
            current_model.elements.append(line)
        else:
            parsed_netlist.append(line)

    final_netlist = unpack(parsed_netlist, subckt_models, operating_point)
    return final_netlist


def unpack(parsed_netlist, subckt_models, operating_point):
    """
    Identifies all subcircuits, unpacks them and returns a unpacked netlist
    Elements inside a subcircuit inherit it's name in the following format: ElementName_SubcircuitName
    :param list parsed_netlist: each cell contains one netlist line in string format
    :return: list unpacked_netlist: each cell contains one netlist line in string format
    """
    final_netlist = []
    for line in parsed_netlist:
        words = line.split()
        name = words[0]
        if name[0].lower() in ["x", "q", "m", "d", "j"]:
            loading_params = False
            nodes = []
            params = {}
            model = None
            for word in words[1:]:
                if word in ["params:", "PARAMS:"]:
                    if model:
                        loading_params = True
                    else:
                        raise NotImplementedError(f"Model of element '{name}' is not present in the netlist.")

                elif not loading_params:
                    try:
                        model = subckt_models[word]
                        params = model.param_dict
                    except KeyError:
                        nodes.append(word)

                else:
                    tmp = word.split("=")
                    if len(tmp) != 2:
                        raise SyntaxError(f"Parameter '{word}' is not in the correct format.")
                    params[tmp[0]], _ = convert_units(tmp[1])
            node_index = 0
            node_dict = {}
            for node in nodes:
                node_dict[model.node_list[node_index]] = node
                node_index += 1

            # TODO: at this stage operating point has to be in format {"gm:Q1": 0.24}
            op = {}
            if operating_point is not None:
                for element in operating_point:
                    parsed_element = element.split(":")
                    if parsed_element[1] == name:
                        op[parsed_element[0]] = operating_point[element]
            else:
                op = None

            elements = model.build_instance(op)
            for elem in elements:
                split_elem = elem.split(" ")
                if split_elem[0][0] not in NETLIST_KEYCHARS: # filter unsupported elements and syntax errors
                    raise NotImplementedError(
                        f"Keyword/Element '{split_elem[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")
                elif split_elem[0] in ["k", "K"]: # parse couplings
                    split_elem[1] = f"{split_elem[1]}_{words[0]}"
                    split_elem[2] = f"{split_elem[2]}_{words[0]}"
                    if split_elem[3][0] == "{":
                        try:
                            split_elem[3] = str(params[split_elem[3][1:-1]])
                        except KeyError:
                            split_elem[3] = str(model.param_dict[split_elem[3][1:-1]])
                else: # parse the rest
                    if split_elem[0][0] in ["f", "F", "h", "H"]: # correct subcircuit CC(C/V)S control voltage source (current sensor) name
                        split_elem[3] = f"{split_elem[3]}_{words[0]}"
                    node_count = nodes_per_element(split_elem[0][0]) # get expected node count
                    for i in range(1, node_count+1): # scan and connect subcircuit nodes with outside circuit nodes
                        if split_elem[i] in ["0", 0]:
                            pass
                        else:
                            try:
                                split_elem[i] = node_dict[split_elem[i]]
                            except KeyError:
                                split_elem[i] = f"{split_elem[i]}_{words[0]}"
                    index = node_count+1
                    for e in split_elem[node_count+1:]: # substitute parameters
                        if e == "":
                            continue

                        name = words[0]
                        if "IC" not in params: # Needed to differentiate between OP of different semiconductor elements
                            params["IC"] = sympy.Symbol(f"IC_{name}")
                        if "ID" not in params: # Needed to differentiate between OP of different semiconductor elements
                            params["ID"] = sympy.Symbol(f"ID_{name}")

                        local = sympy.abc._clash
                        split_items = e.split("=")
                        if len(split_items) == 1:
                            tmp_elem = split_items[0]
                            param_prefix = ""
                        elif len(split_items) == 2:
                            tmp_elem = split_items[1]
                            param_prefix = split_items[0]+'='
                        else:
                            raise SyntaxError(f"Unrecognized syntax in netlist line: {elem}, unrecognized part: {e}")

                        tmp_elem = tmp_elem.replace("{", "")
                        tmp_elem = tmp_elem.replace("}", "")
                        local.update(params)
                        #tmp_elem, _ = convert_units(tmp_elem, local_dict=local)
                        enum_res = value_enum("tmp", tmp_elem, local_dict=local)
                        tmp_elem = enum_res[0]
                        string_tmp = str(tmp_elem).replace(" ", "")
                        split_elem[index] = param_prefix + string_tmp

                        index += 1

                device_name = words[0]
                split_elem = append_dev_name(split_elem, device_name)
                final_netlist.append(" ".join(split_elem))
        else:
            final_netlist.append(line)
    return final_netlist

def append_dev_name(split_elem, device_name):
    if len(split_elem) == 3:
        if split_elem[0] == "rpi":
            split_elem.append(f"1/gpi_{device_name}")
        elif split_elem[0] == "rd":
            split_elem.append(f"1/gd_{device_name}")
        elif split_elem[0] == "ro":
            split_elem.append(f"1/go_{device_name}")

    split_elem[0] = f"{split_elem[0]}_{device_name}"
    return split_elem

def append_lines_plus(netlist_lines):
    new_netlist_lines = []
    for line in netlist_lines:
        if len(line) > 0:
            if line[0] == '+':
                new_netlist_lines[-1] += line[1:]
            else:
                new_netlist_lines.append(line)
        else:
            new_netlist_lines.append(line)
    return new_netlist_lines

def preparse(netlist):
    netlist_lines = netlist.splitlines()
    sanitized_netlist_lines = append_lines_plus(netlist_lines)
    preparsed_netlist_lines = [sanitized_netlist_lines[0]]
    for line in sanitized_netlist_lines[1:]:
        split_line = line.split()
        if split_line:
            first_char = split_line[0][0]
            if first_char in ["m", "M"]:  # add 'PARAMS:' between mos model name and 'W' and 'L' params
                if split_line[6] != "PARAMS:":
                    split_line.insert(6, "PARAMS:")
                preparsed_line = " ".join(split_line)
                preparsed_netlist_lines.append(preparsed_line)
            elif first_char == "*":  # Strip commentaries
                pass
            else:
                preparsed_line = " ".join(split_line)
                preparsed_netlist_lines.append(preparsed_line)
    return preparsed_netlist_lines

def parse(netlist, operating_point=None):
    """
    Translates
    :param str netlist: netlist in a string format
    :return: list data: data contains four items: \n
    * :node_dict: dictionary of circuit nodes
    * :code_count: amount of nodes
    * :matrix_size: matrix size needed to fit all components in
    * :components: list of components
    \n
    Input example: \n
    Circuit AC6
    V1 a 0 dc 0 ac 1 0 sin 0 1 14k 0 0
    R1 a b R1
    L b 0 L1
    R2 b 0 1k
    """
    data = {}
    parsed_netlist = preparse(netlist)
    components = {}
    count = 0
    c = None
    nodes = []
    independent_sources = []
    basic_components = []
    controlled_sources = []
    operational_amplifiers = []
    couplings = []
    SCSI_components = []
    add_short = {}
    parsed_netlist = parse_subcircuits(parsed_netlist, operating_point)

    for line in parsed_netlist:
        words = line.split()
        name = words[0]
        if name[0] == "*":
            continue
        if name[0] not in ['k', 'K']:
            node1 = words[1]
            node2 = words[2]
            symbolic = False
            if node1 not in nodes:
                nodes.append(node1)
            if node2 not in nodes:
                nodes.append(node2)

        # identify variant of component and assign symbolic value

        if name[0] in ["i", "I"]:
            dc_sig, ac_sig, tran_sig = preparse_source(words)

            values = {"dc_raw": dc_sig,
                      "ac_raw": ac_sig,
                      "tran_raw": tran_sig}

            c = CurrentSource(name, node1, node2, values)
            independent_sources.append(c)

        elif name[0] in ["v", "V", "u", "U"]:
            dc_sig, ac_sig, tran_sig = preparse_source(words)

            values = {"dc_raw": dc_sig,
                      "ac_raw": ac_sig,
                      "tran_raw": tran_sig}

            c = VoltageSource(name, node1, node2, values)
            independent_sources.append(c)

        elif name[0] in ["r", "R"]:
            if len(words) < 4:
                val = words[0]
            else:
                val = words[3]
            val_rat, val_flt, val_sym, is_symbolic = value_enum(words[0], val)
            c = Resistor(name, node1, node2, sym_value=val_sym, value=val_rat, value_float=val_flt)
            basic_components.append(c)

        elif name[0] in ["c", "C"]:
            if len(words) < 4:
                val = words[0]
            else:
                val = words[3]
            val_rat, val_flt, val_sym, is_symbolic = value_enum(words[0], val)
            try:
                ic_rat, ic_flt, _, _ = value_enum(words[0], words[4][3:])
            except IndexError:
                ic_rat, ic_flt = 0, 0

            c = Capacitor(name, node1, node2, sym_value=val_sym, init_cond=ic_rat, value=val_rat, value_float=val_flt)
            basic_components.append(c)

        elif name[0] in ["l", "L"]:
            if len(words) < 4:
                val = words[0]
            else:
                val = words[3]
            val_rat, val_flt, val_sym, is_symbolic = value_enum(words[0], val)
            try:
                ic_rat, ic_flt, _, _ = value_enum(words[0], words[4][3:])
            except IndexError:
                ic_rat, ic_flt = 0, 0

            c = Inductor(name, node1, node2, sym_value=val_sym, init_cond=ic_rat, value=val_rat, value_float=val_flt)
            basic_components.append(c)

        elif name[0] in ["a", "A"]:
            node3 = words[3]
            node4 = words[4]
            if node3 not in nodes:
                nodes.append(node3)
            if node4 not in nodes:
                nodes.append(node4)
            c = IdealOperationalAmplifier(name, node1, node2, node3, node4)
            operational_amplifiers.append(c)

        elif name[0] in ["e", "E"]:  # VVT (VCVS)
            variant = "e"
            node3 = words[3]
            node4 = words[4]
            if node3 not in nodes:
                nodes.append(node3)
            if node4 not in nodes:
                nodes.append(node4)

            try:
                val = words[5]
            except IndexError:
                val = words[0]

            val_rat, val_flt, val_sym, is_symbolic = value_enum(words[0], val)

            c = VoltageControlledSource(name, variant, node1, node2, node3, node4, value=val_rat,
                                        value_float=val_flt, sym_value=val_sym)
            controlled_sources.append(c)

        elif name[0] in ["g", "G"]:  # VCT (VCCS)
            variant = "g"
            node3 = words[3]
            node4 = words[4]
            if node3 not in nodes:
                nodes.append(node3)
            if node4 not in nodes:
                nodes.append(node4)
            try:
                val = words[5]
            except IndexError:
                val = words[0]

            val_rat, val_flt, val_sym, is_symbolic = value_enum(words[0], val)

            c = VoltageControlledSource(name, variant, node1, node2, node3, node4, value=val_rat,
                                        value_float=val_flt, sym_value=val_sym)
            controlled_sources.append(c)

        elif name[0] in ["f", "F"]:  # CCT (CCCS)
            variant = "f"
            v_c = words[3]
            try:
                add_short[v_c].append(name)
                node3 = f"*ctrl{v_c}{add_short[v_c][-2]}"
                node4 = f"*ctrl{v_c}{add_short[v_c][-1]}"
            except:
                add_short[v_c] = [name]
                node3 = None
                node4 = f"*ctrl{v_c}{add_short[v_c][0]}"

            try:
                val = words[4]
            except IndexError:
                val = words[0]

            val_rat, val_flt, val_sym, is_symbolic = value_enum(words[0], val)

            c = CurrentControlledSource(name, variant, node1, node2, node3, node4, value=val_rat,
                                        value_float=val_flt, sym_value=val_sym)
            controlled_sources.append(c)

        elif name[0] in ["h", "H"]:  # CVT (CCVS)
            variant = "h"
            v_c = words[3]
            try:
                add_short[v_c].append(name)
                node3 = f"*ctrl{v_c}{add_short[v_c][-2]}"
                node4 = f"*ctrl{v_c}{add_short[v_c][-1]}"
            except:
                add_short[v_c] = [name]
                node3 = None
                node4 = f"*ctrl{v_c}{add_short[v_c][0]}"

            try:
                val = words[4]
            except IndexError:
                val = words[0]

            val_rat, val_flt, val_sym, is_symbolic = value_enum(words[0], val)

            c = CurrentControlledSource(name, variant, node1, node2, node3, node4, value=val_rat,
                                        value_float=val_flt, sym_value=val_sym)
            controlled_sources.append(c)

        elif name[0] in ["k", "K"]:  # coupled inductors
            val_rat, val_flt, val_sym, is_symbolic = value_enum(words[0], words[3])
            L1 = words[1]
            L2 = words[2]
            c = Coupling(name, L1, L2, val_sym, value=val_rat, value_float=val_flt)
            couplings.append(c)

        elif name[0] in ["s", "S"]:  # periodic switch used for SC/SI analysis
            phase = words[3]
            c = PeriodicSwitch(name, node1, node2, phase=phase)
            SCSI_components.append(c)

        components[c.name] = c
        count += 1

    for key in add_short:
        c_s = components[key]
        new_node = None

        shorted_node = c_s.node1
        for c_name in add_short[key]:
            c = components[c_name]
            if c.node3 is None:
                c.node3 = shorted_node
            new_node = f"*ctrl{c_s.name}{c_name}"
            nodes.append(new_node)

        c_s.node1 = new_node

    for couple in couplings:
        L1 = components[couple.L1]
        L2 = components[couple.L2]
        L1.coupling = couple
        L2.coupling = couple

    data["components"] = components
    data["basic_components"] = basic_components
    data["independent_sources"] = independent_sources
    data["controlled_sources"] = controlled_sources
    data["oamps"] = operational_amplifiers
    data["couplings"] = couplings
    data["SCSI_components"] = SCSI_components
    return components, couplings
