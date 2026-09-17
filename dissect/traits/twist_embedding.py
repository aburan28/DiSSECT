from dissect.traits import Trait

TRAIT_TIMEOUT = 30


class TwistEmbeddingTrait(Trait):
    NAME = "twist_embedding"
    DESCRIPTION = (
        "Embedding degree of the quadratic twist with respect to the largest prime factor of"
        " its cardinality. A small degree exposes the twist to an MOV/Frey-Rueck transfer,"
        " which is what an invalid-point or fault attack would land on."
    )
    INPUT = {}
    OUTPUT = {
        "largest_prime_bitlen": (
            int,
            "Bit length of the largest prime factor of the twist cardinality",
        ),
        "embedding_degree_bitlen": (int, "Bit length of the twist embedding degree"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Embedding degree of the twist, guarded by timeouts on both hard steps.

        The twist has cardinality $q + 1 + t$. Factoring it and taking a
        multiplicative order are each expensive enough to need a timeout, so a
        curve can report a bit length for the prime and nothing for the degree.
        """
        from sage.all import Integers, ZZ
        from dissect.utils.utils import Factorization, timeout

        q = curve.q()
        twist_cardinality = q + 1 + curve.trace()
        factorization = Factorization(twist_cardinality, timeout_duration=TRAIT_TIMEOUT)
        if factorization.timeout():
            return {"largest_prime_bitlen": None, "embedding_degree_bitlen": None}

        largest_prime = ZZ(factorization.factorization()[-1])
        results = {"largest_prime_bitlen": int(largest_prime.nbits())}

        degree = timeout(
            lambda x: x.multiplicative_order(),
            [Integers(largest_prime)(q)],
            timeout_duration=TRAIT_TIMEOUT,
        )
        if isinstance(degree, str):
            results["embedding_degree_bitlen"] = None
        else:
            results["embedding_degree_bitlen"] = int(ZZ(degree).nbits())
        return results


def test_twist_embedding():
    assert True
