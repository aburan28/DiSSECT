from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, AnalysisUnavailable
from dissect.utils.curve_automorphisms import extension_field, hyperelliptic_odd, hyperelliptic_even, plane_automorphisms, trigonal_automorphisms


class AutomorphismStructureTrait(DetailedTrait):
    NAME = "automorphism_structure"
    DESCRIPTION = "Finite-field automorphism groups, canonical actions and tame cyclic quotient genera."
    INPUT = dict(DetailedTrait.INPUT, extension_degree=(int, "Constant-field extension in which automorphisms are sought"),
                 max_candidates=(int, "Projective or affine search budget"), max_transformations=(int, "Maximum explicit transformations returned"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, extension_degree=[1], max_candidates=[100000], max_transformations=[128])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "search_field_size": (int, "Cardinality of the field over which automorphisms are defined"),
        "search_field_extension_degree": (int, "Relative extension degree over the input field"),
        "search_field_modulus": (str, "Power-basis modulus for serialized extension-field elements"),
        "source_field_generator_image": (list, "Image of the input constant-field generator in the search field basis"),
        "base_field_order": (int, "Exact order over the search field, or null if incomplete"),
        "computed_subgroup_order": (int, "Certified subgroup order, a lower bound if incomplete"),
        "search_complete": (int, "Whether all automorphisms over this finite field were enumerated"),
        "cover_stabilizer_complete": (int, "Whether the chosen trigonal-cover stabilizer was exhausted"),
        "group_structure": (str, "GAP structure description of the computed group"),
        "abelian": (int, "Whether the computed group is abelian"),
        "center_order": (int, "Order of its center"),
        "element_order_distribution": (list, "Exact order distribution in the computed group"),
        "reduced_automorphism_order": (int, "Image order on the hyperelliptic quotient"),
        "branch_stabilizer_order": (int, "PGL2 branch stabilizers, including transformations without rational lifts"),
        "transformations": (list, "Verified coordinate transformations and regular differential actions"),
        "transformations_complete": (int, "Whether all computed transformations fit the output bound"),
        "canonical_action_faithful": (int, "0 for characteristic-two hyperelliptic deck involutions"),
        "tame_cyclic_quotient_genera": (list, "Invariant-differential dimensions only for groups prime to p"),
        "geometric_automorphism_order": (int, "Null: a finite extension search does not certify the algebraic-closure group"),
        "hyperelliptic_quotient_genus": (int, "Genus zero for the deck involution, including wild double covers"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        if int(curve.q())**params["extension_degree"] > params["max_field_size"]:
            raise AnalysisUnavailable("Search extension exceeds max_field_size", "resource_limit")
        field, embedding = extension_field(curve.field(), params["extension_degree"])
        if curve.model_type() == "CyclicTrigonal":
            operation = trigonal_automorphisms
        elif curve.is_hyperelliptic():
            operation = hyperelliptic_even if field.characteristic() == 2 else hyperelliptic_odd
        else:
            operation = plane_automorphisms
        result = operation(curve, field, embedding, params)
        from dissect.utils.curve_analysis import field_value

        generator_image = field_value(embedding(curve.field().gen()))
        if not isinstance(generator_image, list):
            generator_image = [generator_image]
        result.update(
            search_field_size=int(field.cardinality()), search_field_extension_degree=params["extension_degree"],
            search_field_modulus=str(field.modulus()),
            source_field_generator_image=generator_image,
            transformations_complete=int(len(result["transformations"]) == result["computed_subgroup_order"]),
            geometric_automorphism_order=None,
            hyperelliptic_quotient_genus=0 if curve.is_hyperelliptic() else None,
            status="ok" if result["search_complete"] else "partial",
        )
        return result
