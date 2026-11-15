from typing import List

from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, function_field, semigroup_profile


class WeierstrassGapsTrait(DetailedTrait):
    NAME = "weierstrass_gaps"
    DESCRIPTION = "Generic gaps, all closed Weierstrass places and numerical semigroups, including nonclassical curves."
    INPUT = dict(DetailedTrait.INPUT, max_places=(int, "Maximum place records returned"), rational_samples=(int, "Maximum additional rational places"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_places=[256], rational_samples=[16])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "generic_gaps": (List[int], "Generic geometric gap sequence"),
        "classical": (int, "Generic gaps are 1 through genus"),
        "weierstrass_places": (list, "Closed places, degrees, gaps, gap weights and semigroup generators"),
        "weierstrass_closed_points": (int, "Total number of closed Weierstrass points"),
        "weierstrass_geometric_points": (int, "Sum of closed-point degrees"),
        "rational_weierstrass_points": (int, "Number of degree-one Weierstrass places"),
        "rational_samples": (list, "Gap profiles at additional base-field rational places"),
        "returned_locus_complete": (int, "Whether the full computed locus fits max_places"),
        "weighted_gap_excess": (int, "Degree-weighted excess above generic gaps"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        L = function_field(curve)
        generic = [int(n) for n in L.gaps()]
        places = sorted(L.weierstrass_places(), key=lambda P: (int(P.degree()), str(P)))

        def record(place):
            profile = semigroup_profile(place.gaps())
            profile.update(place=str(place), degree=int(place.degree()),
                           generic_gap_excess=sum(profile["gaps"])-sum(generic))
            return profile

        records = [record(P) for P in places]
        special = set(places)
        ordinary = [P for P in L.places(1) if P not in special][:params["rational_samples"]]
        return {
            "generic_gaps": generic, "classical": int(generic == list(range(1, curve.genus()+1))),
            "weierstrass_places": records[:params["max_places"]],
            "weierstrass_closed_points": len(places),
            "weierstrass_geometric_points": sum(int(P.degree()) for P in places),
            "rational_weierstrass_points": sum(P.degree() == 1 for P in places),
            "rational_samples": [record(P) for P in ordinary],
            "returned_locus_complete": int(len(places) <= params["max_places"]),
            "weighted_gap_excess": sum(r["degree"]*r["generic_gap_excess"] for r in records),
            "status": "ok" if len(places) <= params["max_places"] else "partial",
            "reason": None if len(places) <= params["max_places"] else "Place records truncated by max_places; totals cover the full locus",
        }
