"""Bounded computations with machine-readable applicability and failure states."""

from dissect.traits.trait import Trait
from dissect.utils.curve_analysis import AnalysisUnavailable


class DetailedTrait(Trait):
    SUPPORTED_CURVE_TYPES = ("hyperelliptic", "non_hyperelliptic")
    GENUS = None
    MODEL_TYPES = None
    INPUT = {
        "timeout": (int, "Wall-clock limit in seconds"),
        "max_field_size": (int, "Maximum base-field cardinality"),
        "max_genus": (int, "Maximum genus"),
    }
    DEFAULT_PARAMS = {"timeout": [30], "max_field_size": [4096], "max_genus": [12]}
    OUTPUT = {
        "status": (str, "ok, partial, unsupported, resource_limit, timeout, or not_applicable"),
        "reason": (str, "Explanation of any limitation"),
    }

    def __call__(self, curve, **params):
        options = {key: values[0] for key, values in self.DEFAULT_PARAMS.items()}
        unknown = set(params) - set(options)
        if unknown:
            raise ValueError("Unknown trait parameters: " + ", ".join(sorted(unknown)))
        options.update(params)
        for key, value in options.items():
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError("Trait parameter %s must be a positive integer" % key)
        result = dict.fromkeys(self.OUTPUT)
        family = getattr(curve, "curve_type", lambda: "elliptic")()
        if family not in self.SUPPORTED_CURVE_TYPES:
            return dict(result, status="not_applicable", reason="Requires a supported higher-genus model")
        if self.GENUS is not None and curve.genus() != self.GENUS:
            return dict(result, status="not_applicable", reason="Requires genus %d" % self.GENUS)
        if self.MODEL_TYPES is not None and curve.model_type() not in self.MODEL_TYPES:
            return dict(result, status="not_applicable", reason="Requires one of these models: " + ", ".join(self.MODEL_TYPES))
        from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt

        try:
            alarm(options["timeout"])
            values = self.compute(curve, options)
            result.update(status="ok", reason=None)
            result.update(values)
        except AnalysisUnavailable as error:
            result.update(status=error.status, reason=str(error))
        except AlarmInterrupt:
            result.update(status="timeout", reason="Exceeded the configured wall-clock limit")
        except NotImplementedError as error:
            result.update(status="unsupported", reason=str(error) or "Sage backend does not implement this computation")
        finally:
            cancel_alarm()
        return result
