from .trait import Trait
from .cofactor import CofactorTrait
from .discriminant import DiscriminantTrait
from .twist_order import TwistOrderTrait
from .kn_factorization import KNFactorizationTrait
from .torsion_extension import TorsionExtensionTrait
from .conductor import ConductorTrait
from .embedding import EmbeddingTrait
from .class_number import ClassNumberTrait
from .small_prime_order import SmallPrimeOrderTrait
from .division_polynomials import DivisionPolynomialsTrait
from .volcano import VolcanoTrait
from .isogeny_extension import IsogenyExtensionTrait
from .trace_factorization import TraceFactorizationTrait
from .isogeny_neighbors import IsogenyNeighborsTrait
from .q_torsion import QTorsionTrait
from .hamming_x import HammingXTrait
from .square_4p1 import Square4P1Trait
from .pow_distance import PowDistanceTrait
from .multiples_x import MultiplesXTrait
from .x962_invariant import X962InvariantTrait
from .brainpool_overlap import BrainpoolOverlapTrait
from .weierstrass import WeierstrassTrait
from .extension_degree import ExtensionDegreeTrait
from .subfield_curve import SubfieldCurveTrait
from .automorphisms import AutomorphismsTrait
from .twist_embedding import TwistEmbeddingTrait
from .sato_tate import SatoTateTrait
from .prime_shape import PrimeShapeTrait
from .group_structure import GroupStructureTrait
from .montgomery_form import MontgomeryFormTrait
from .cm_discriminant_size import CMDiscriminantSizeTrait
from .field_representation import FieldRepresentationTrait
from .cheon import CheonTrait
from .canonical_geometry import CanonicalGeometryTrait
from .hyperelliptic_branch import HyperellipticBranchTrait
from .frobenius_polynomial import FrobeniusPolynomialTrait
from .jacobian_decomposition import JacobianDecompositionTrait
from .hasse_witt import HasseWittTrait
from .newton_polygon import NewtonPolygonTrait
from .weierstrass_gaps import WeierstrassGapsTrait
from .automorphism_structure import AutomorphismStructureTrait
from .moduli_invariants import ModuliInvariantsTrait
from .gonality_and_quotients import GonalityAndQuotientsTrait
from .cover_ramification import CoverRamificationTrait
from .genus2_rosenhain import Genus2RosenhainTrait
from .genus2_elliptic_quotients import Genus2EllipticQuotientsTrait
from .genus3_bitangents import Genus3BitangentsTrait
from .genus3_flexes import Genus3FlexesTrait
from .genus4_canonical_quadric import Genus4CanonicalQuadricTrait
from .genus4_tritangents import Genus4TritangentsTrait

TRAITS = dict(
    map(
        lambda x: (x.NAME, x()),
        [
            CofactorTrait,
            DiscriminantTrait,
            TwistOrderTrait,
            KNFactorizationTrait,
            TorsionExtensionTrait,
            ConductorTrait,
            EmbeddingTrait,
            ClassNumberTrait,
            SmallPrimeOrderTrait,
            DivisionPolynomialsTrait,
            VolcanoTrait,
            IsogenyExtensionTrait,
            TraceFactorizationTrait,
            IsogenyNeighborsTrait,
            QTorsionTrait,
            HammingXTrait,
            Square4P1Trait,
            PowDistanceTrait,
            MultiplesXTrait,
            X962InvariantTrait,
            BrainpoolOverlapTrait,
            WeierstrassTrait,
            ExtensionDegreeTrait,
            SubfieldCurveTrait,
            AutomorphismsTrait,
            TwistEmbeddingTrait,
            SatoTateTrait,
            PrimeShapeTrait,
            GroupStructureTrait,
            MontgomeryFormTrait,
            CMDiscriminantSizeTrait,
            FieldRepresentationTrait,
            CheonTrait,
            CanonicalGeometryTrait,
            HyperellipticBranchTrait,
            FrobeniusPolynomialTrait,
            JacobianDecompositionTrait,
            HasseWittTrait,
            NewtonPolygonTrait,
            WeierstrassGapsTrait,
            AutomorphismStructureTrait,
            ModuliInvariantsTrait,
            GonalityAndQuotientsTrait,
            CoverRamificationTrait,
            Genus2RosenhainTrait,
            Genus2EllipticQuotientsTrait,
            Genus3BitangentsTrait,
            Genus3FlexesTrait,
            Genus4CanonicalQuadricTrait,
            Genus4TritangentsTrait,
        ],
    )
)
