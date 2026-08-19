from mrfmsim import formula
from mrfmsim import ExperimentGroup, Node
from mrfmsim.experiment.stdelements import STANDARD_NODES, STANDARD_COMPONENTS_FORCE

# developing new cermittd-force experiment simulator



node_objects = [
    Node("rel_dpol td_sat mw", formula.rel_dpol_sat_td_mw, output="rel_dpol"),
    Node("phase one pulse", formula.phase_one_pulse, output="phase_delay"),
    Node("positiveFeedback_efficiency_halfcycle", formula.positiveFeedback_efficiency_halfcycle, output="positivefeedback_eff"),
    Node("positiveFeedback_efficiency_onetime", formula.positiveFeedback_efficiency_onetime, output="positivefeedback_eff"),
    Node(
        "force td",
        formula.force_td,
        inputs=["Bzx", "positivefeedback_eff", "mz_eq", "spin_density", "grid_voxel",'Q',"k_c"],
        output="amp_spin",
        doc="Calculate force account for the negative sign in the approximation.",
    ),
]

CermitTD_force_mul_edges = [
    ["Bz", "B_tot"],
    ["B_tot", ["mz_eq", "B_offset"]],
    [["B_offset"], "rel_dpol td_sat mw"],
    ["rel_dpol td_sat mw", "positiveFeedback_efficiency_halfcycle"],
    [["mz_eq", "Bzx", "positiveFeedback_efficiency_halfcycle"], "force td"],
]


CermitTD_force_one_edges = [
    ["Bz", "B_tot"],
    ["B_tot", ["mz_eq", "B_offset"]],
    [["B_offset"], ["rel_dpol td_sat mw","phase one pulse"]],
    [["rel_dpol td_sat mw","phase one pulse"], "positiveFeedback_efficiency_onetime"],
    [["mz_eq", "Bzx", "positiveFeedback_efficiency_onetime"], "force td"],
]

experiment_recipes = {
    "CermitTD_force_multipulse": {
        "grouped_edges": CermitTD_force_mul_edges,
        "doc": "Time-dependent CERMIT experiment for a large tip. Saturating spin with external microwave, detecting force signal",
    },
    "CermitTD_force_onepulse": {
        "grouped_edges": CermitTD_force_one_edges,
        "doc": "Time-dependent CERMIT experiment for a large tip. Saturating spin with external microwave, detecting force signal assuming saturation only happens once per cantilever cycle",
    },
}

docstring = """\
Simulates a Cornell-style frequency shift magnetic
resonance force microscope experiment considering the time-dependent
nature of the saturation, averaged over multiple pulses and with
small-step approximation.
"""

CermitTDGroup_force = ExperimentGroup(
    name="CermitTDGroup_force",
    node_objects=list(STANDARD_NODES) + node_objects,
    experiment_recipes=experiment_recipes,
    experiment_defaults={"components": STANDARD_COMPONENTS_FORCE},
    doc=docstring,
)