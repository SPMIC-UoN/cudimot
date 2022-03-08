"""
CUDIMOT GUI: Model parameters setup
"""
import os

import wx
import wx.grid

from . import OptionError
from .widgets import TabPage, ParameterList, NumberChooser, NumberList

class ModelParameters(TabPage):
    """
    Tab page containing definitions of model parameters
    """

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Model parameters", idx, n, name="params")

        self.section("Variable parameters")
        self.ptables = []
        self.ptables.append(ParameterList(self, variable=True))
        self.pack("", self.ptables[0], span=3, expand=True)
        self.section("Common fixed parameters")
        self.ptables.append(ParameterList(self, variable=False))
        self.pack("", self.ptables[1], span=3, expand=True)
        self.section("Voxelwise fixed parameters")
        self.ptables.append(ParameterList(self, variable=False))
        self.pack("", self.ptables[2], span=3, expand=True)
        self.add_next_prev_btn()

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.Fit()
        self.Layout()
        self.Bind(wx.EVT_SIZE, self._resize_tables)

    def _resize_tables(self, event):
        width = self.GetClientSize()[0] - 10
        for table in self.ptables:
            table.ResizeCols(width)
        event.Skip()

    def reset(self, projdir):
        self.ptables[0].ClearParameters()
        self.ptables[1].ClearParameters()
        self.ptables[2].ClearParameters()

    def config(self):
        config = {
            "params" : self.ptables[0].GetParameters(),
            "common_fixed_params" : self.ptables[1].GetParameters(),
            "voxelwise_fixed_params" : self.ptables[2].GetParameters(),
        }
        for idx, ptype in enumerate(["param", "cfp", "vfp"]):
            params = self.ptables[idx].GetParameters()
            config[f"n_{ptype}"] = len(params)
            if params:
                for key in params[0].keys():
                    config[f"{ptype}_{key}s"] = ",".join([str(p[key]) for p in params])
            else:
                config[f"{ptype}_names"] = ""
                config[f"{ptype}_sizes"] = ""
        return config

    def load(self, projdir):
        params_file = os.path.join(projdir, "modelparameters.h")
        priors_file = os.path.join(projdir, "modelpriors")
        param_names = self._parse_params(params_file, "NPARAMS")
        param_inits, param_bounds, param_priors = self._parse_bounds_inits(priors_file, len(param_names))
        self.ptables[0].ClearParameters()
        for name, init, bounds, prior in zip(param_names, param_inits, param_bounds, param_priors):
            self.ptables[0].AddParameter(name=name, size=1, init=init, lbound=bounds[0], ubound=bounds[1], prior=prior)

        cfp_names = self._parse_params(params_file, "NCFP")
        param_sizes_file = os.path.join(projdir, "modelparameters.cc")
        cfp_sizes = self._parse_param_sizes(param_sizes_file, "CFP", len(cfp_names))
        self.ptables[1].ClearParameters()
        for param_name, param_size in zip(cfp_names, cfp_sizes):
            self.ptables[1].AddParameter(name=param_name, size=param_size)

        fixp_names = self._parse_params(params_file, "NFIXP")
        fixp_sizes = self._parse_param_sizes(param_sizes_file, "FixP", len(fixp_names))
        self.ptables[2].ClearParameters()
        for param_name, param_size in zip(fixp_names, fixp_sizes):
            self.ptables[2].AddParameter(name=param_name, size=param_size)

    def _parse_params(self, fname, label):
        nparams = self.config_from_line_regex("nparams", fname, 
                                              f"#define\s+{label}\s+(\d+)")
        param_names = self.config_from_line_regex("param_names", fname, 
                                              f"#define\s+{label}\s+\d+\s*//\s*(.+)")

        try:
            nparams = int(nparams)
        except ValueError:
            print(f"Unable to parse number of parameters for {label}: {nparams}")
            nparams = 0

        if nparams:
            param_names = [n.strip() for n in param_names.split(",")]
            if len(param_names) != nparams:
                print(f"Wrong number of parameter names: {nparams} vs {param_names}")
            for idx in range(len(param_names), nparams):
                param_names.append("param%i" % (idx+1))

            return param_names[:nparams]
        else:
            return []

    def _parse_param_sizes(self, fname, label, nparams):
        sizes = self.config_from_line_regex("param_sizes", fname, 
                                              f"int\s+MODEL::{label}_size\s*\[\]\s*=\s*{{([\d,]*)}};")
        print(sizes)
        try:
            if not sizes.strip():
                sizes = []
            else:
                sizes = [int(s) for s in sizes.split(",")]
            print(sizes)
            if len(sizes) != nparams:
                print(f"Wrong number of parameter names: {nparams} vs {sizes}")
            for idx in range(len(sizes), nparams):
                sizes.append(1)
            return sizes
        except ValueError:
            print(f"Unable to parse parameter sizes for {label}: {sizes}")
            return [1] * nparams

    def _parse_bounds_inits(self, fname, nparams):
        inits = self.config_from_line_regex("p_init", fname, "P_init\s*\[([\d,\.\-\s]+)\]")
        try:
            print(inits)
            inits = [float(s) for s in inits.split(",")]
            print(inits)
            if len(inits) != nparams:
                print(f"Wrong number of parameter bounds: {nparams} vs {inits}")
            for idx in range(len(inits), nparams):
                inits.append(0.0)
        except ValueError:
            print(f"Unable to parse parameter inits: {inits}")
            inits = [0.0] * nparams

        param_bounds = []
        for idx in range(nparams):
            bounds = self.config_from_line_regex("bounds", fname, "bounds\s*\[\s*%i\s*\]\s*=\s*\(([\d,\.\-\s]+)\)" % idx)
            try:
                bounds = [float(s) for s in bounds.split(",")]
                if len(bounds) != 2:
                    print(f"Unable to parse parameter bounds: {bounds}")
                    param_bounds.append((-10000, 10000))
                else:
                    param_bounds.append(bounds)
            except ValueError:
                print(f"Unable to parse parameter bounds for parameter %i: {bounds}" % idx)
                param_bounds.append((-10000, 10000))

        param_priors = []
        for idx in range(nparams):
            prior = self.config_from_line_regex("priors", fname, "priors\s*\[\s*%i\s*\]\s*=\s*(.+)" % idx)
            if not prior:
                prior = "None"
            param_priors.append(prior)

        return inits, param_bounds, param_priors
