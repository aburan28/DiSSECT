from dissect.traits import Trait


def hamming_weight(value):
    """Number of set bits."""
    return bin(value).count("1")


def naf_weight(value):
    """Number of nonzero digits in the non-adjacent form of `value`.

    NAF is the signed binary representation with no two adjacent nonzero
    digits, and it is minimal in nonzero digits among signed representations.
    A field prime chosen for fast reduction -- Solinas, pseudo-Mersenne -- is
    a short sum of signed powers of two and so has a tiny NAF weight, where a
    prime with no such structure sits near bitlength/3.
    """
    weight, remaining = 0, value
    while remaining:
        if remaining & 1:
            digit = 2 - (remaining & 3)  # +1 or -1, forcing the next bit to zero
            remaining -= digit
            weight += 1
        remaining >>= 1
    return weight


class PrimeShapeTrait(Trait):
    NAME = "prime_shape"
    DESCRIPTION = (
        "Sparsity of the field size in binary and in non-adjacent form. A prime picked for fast"
        " reduction is a short signed sum of powers of two and has a very small NAF weight;"
        " one with no such structure has roughly bitlength/3."
    )
    INPUT = {}
    OUTPUT = {
        "hamming_weight": (int, "Set bits in $q$"),
        "naf_weight": (int, "Nonzero digits in the non-adjacent form of $q$"),
        "naf_density": (float, "NAF weight divided by bit length"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Sparsity of the field size.

        This describes the field, not the curve, so within any one simulated
        pool -- which varies the seed over a fixed prime -- it is constant and
        contributes nothing to a standard-versus-simulated comparison. Its use
        is screening across the database.
        """
        q = int(curve.q())
        weight = naf_weight(q)
        return {
            "hamming_weight": hamming_weight(q),
            "naf_weight": weight,
            "naf_density": weight / q.bit_length(),
        }


def test_hamming_weight():
    assert hamming_weight(0) == 0
    assert hamming_weight(1) == 1
    assert hamming_weight(2**64 - 1) == 64


def test_naf_weight_of_powers_and_neighbours():
    assert naf_weight(0) == 0
    assert naf_weight(2**100) == 1
    # 2^100 - 1 is 100 set bits in binary but two signed digits in NAF.
    assert hamming_weight(2**100 - 1) == 100
    assert naf_weight(2**100 - 1) == 2


def test_naf_weight_no_adjacent_nonzeros():
    # 0b1011 = 11 = 8 + 4 - 1, three digits, not four.
    assert naf_weight(11) == 3
    assert naf_weight(7) == 2


def test_solinas_prime_is_sparse():
    # NIST P-256: p = 2^256 - 2^224 + 2^192 + 2^96 - 1, five signed terms.
    p256 = 2**256 - 2**224 + 2**192 + 2**96 - 1
    assert naf_weight(p256) == 5
    # Curve25519: p = 2^255 - 19.
    assert naf_weight(2**255 - 19) == 4


def test_unstructured_prime_is_dense():
    # A Brainpool prime is generated from pi and has no reduction structure.
    bp256 = 76884956397045344220809746629001649093037950200943055203735601445031516197751
    assert naf_weight(bp256) > bp256.bit_length() // 4
