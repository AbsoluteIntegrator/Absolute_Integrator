import numpy as np

# dictionary describing options available to tune this algorithm
options = {
    "peak_size": {"purpose": "Estimate of the peak size, in pixels.  If 'auto', attempts to determine automatically.  Otherwise, this should be an integer.",
                  "default": "auto",
                  "type": "int",
                  "has_auto": True},
    "refine_positions": {"purpose": "TODO",
                         "default": False,
                         "type": "bool"},
    "progress_object": {"purpose": "Object used to present a progress bar to the user.  For definition, see UI_interface folder.",
                        "default": None},
}

def run(data):
    # TODO: need to actually implement this peak finder.
    return np.zeros((4,2))