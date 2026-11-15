import itertools


class Trait:
    NAME = "trait"
    DESCRIPTION = "Trait description."
    INPUT = {}
    OUTPUT = {}
    DEFAULT_PARAMS = {}
    SUPPORTED_CURVE_TYPES = ("elliptic",)

    def __init__(self):
        pass

    def __call__(self, curve, **params):
        curve_type = getattr(curve, "curve_type", lambda: "elliptic")()
        if curve_type not in self.SUPPORTED_CURVE_TYPES:
            return dict.fromkeys(self.OUTPUT)
        return self.compute(curve, params)

    def compute(self, curve, params):
        raise NotImplementedError("Compute method for trait not implemented")

    def params(self):
        return dict(self.DEFAULT_PARAMS)

    def params_iter(self):
        for params_values in itertools.product(*self.DEFAULT_PARAMS.values()):
            yield dict(zip(self.DEFAULT_PARAMS.keys(), params_values))

    def numeric_outputs(self):
        outputs = []

        for output, info in self.OUTPUT.items():
            if info[0] == int or info[0] == float:
                outputs.append(output)

        return outputs

    def nonnumeric_outputs(self):
        outputs = []

        for output, info in self.OUTPUT.items():
            if not (info[0] == int or info[0] == float):
                outputs.append(output)

        return outputs

    def outputs(self):
        return list(self.OUTPUT.keys())
