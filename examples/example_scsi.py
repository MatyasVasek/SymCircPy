from symcirc.analysis import Circuit, TF
from symcirc.analysis_scsi import SwitchedAnalysis
from symcirc.utils import pprint

netlist = """
C1 node1 0 10n 
C2 nodeout 0 10n 
S1 node1 nodeout 2 
S2 nodein node1 1 
V1 nodein 0 DC 1 AC 1 {0*Pi/180}
"""

circuit = Circuit(netlist)

print(circuit.components)
switched_analysis = SwitchedAnalysis(circuit, phases='P=2', mode="scideal", method="two_graph_node", symbolic=True)
pprint(switched_analysis.eqn_matrix)
results = switched_analysis.get_all_results()

for i in results:
    print(f"{i}: {results[i]}")

# Results should be this:
""" 
v(1)[1] = V1
v(1)[2] = C1*V1/(C1*z+C2*z-C2)*z^(1/2)
v(in)[1] = V1
v(in)[2] = V1
v(out)[1] = C1*V1/(C1*z+C2*z-C2)
v(out)[2] = C1*V1/(C1*z+C2*z-C2)*z^(1/2)
v(C1)[1] = V1
v(C1)[2] = C1*V1/(C1*z+C2*z-C2)*z^(1/2)
q(C1)[1] = C1*V1
q(C1)[2] = V1*C1^2/(C1*z+C2*z-C2)*z^(1/2)
v(C2)[1] = C1*V1/(C1*z+C2*z-C2)
v(C2)[2] = C1*V1/(C1*z+C2*z-C2)*z^(1/2)
q(C2)[1] = C2*C1*V1/(C1*z+C2*z-C2)
q(C2)[2] = C2*C1*V1/(C1*z+C2*z-C2)*z^(1/2)
v(S1)[1] = V1*(C1*z+C2*z-C1-C2)/(C1*z+C2*z-C2)
v(S1)[2] = 0
q(S1)[1] = 0
q(S1)[2] = C2*(z-1)*C1*V1/(C1*z+C2*z-C2)/z^(1/2)
v(S2)[1] = 0
v(S2)[2] = -C1*V1/(C1*z+C2*z-C2)*z^(1/2)+V1
q(S2)[1] = C1*V1*(C1*z+C2*z-C1-C2)/(C1*z+C2*z-C2)
q(S2)[2] = 0
v(V1)[1] = V1
v(V1)[2] = V1
q(V1)[1] = -C1*V1*(C1*z+C2*z-C1-C2)/(C1*z+C2*z-C2)
q(V1)[2] = 0
"""
