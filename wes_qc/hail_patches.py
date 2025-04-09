"""
This file contains modified versions of Hail functions
used when the upstream Hail release contains some issues
"""

import hail as hl
from hail.matrixtable import MatrixTable
from hail.table import Table
from hail.typecheck import oneof, typecheck


"""
The modified version of Hail split_multi_hts() function.
Original code is here: https://hail.is/docs/0.2/_modules/hail/methods/statgen.html#split_multi_hts
"""


@typecheck(
    ds=oneof(Table, MatrixTable),
    keep_star=bool,
    left_aligned=bool,
    vep_root=str,
    permit_shuffle=bool,
    recalculate_gq=bool,
)
def split_multi_hts(
    ds, keep_star=False, left_aligned=False, vep_root="vep", *, permit_shuffle=False, recalculate_gq=False
):
    """Split multiallelic variants for datasets that contain one or more fields
    from a standard high-throughput sequencing entry schema.

    .. code-block:: text

      struct {
        GT: call,
        AD: array<int32>,
        DP: int32,
        GQ: int32,
        PL: array<int32>,
        PGT: call,
        PID: str
      }

    For other entry fields, write your own splitting logic using
    :meth:`.MatrixTable.annotate_entries`.

    Examples
    --------

    >>> hl.split_multi_hts(dataset).write('output/split.mt')

    Warning
    -------
    This method assumes `ds` contains at most one non-split variant per locus. This assumption permits the
    most efficient implementation of the splitting algorithm. If your queries involving `split_multi_hts`
    crash with errors about out-of-order keys, this assumption may be violated. Otherwise, this
    warning likely does not apply to your dataset.

    If each locus in `ds` contains one multiallelic variant and one or more biallelic variants, you
    can filter to the multiallelic variants, split those, and then combine the split variants with
    the original biallelic variants.

    For example, the following code splits a dataset `mt` which contains a mixture of split and
    non-split variants.

    >>> bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    >>> bi = bi.annotate_rows(a_index=1, was_split=False)
    >>> multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    >>> split = hl.split_multi_hts(multi)
    >>> mt = split.union_rows(bi)

    Notes
    -----

    We will explain by example. Consider a hypothetical 3-allelic
    variant:

    .. code-block:: text

      A   C,T 0/2:7,2,6:15:45:99,50,99,0,45,99

    :func:`.split_multi_hts` will create two biallelic variants (one for each
    alternate allele) at the same position

    .. code-block:: text

      A   C   0/0:13,2:15:45:0,45,99
      A   T   0/1:9,6:15:50:50,0,99

    Each multiallelic `GT` or `PGT` field is downcoded once for each alternate allele. A
    call for an alternate allele maps to 1 in the biallelic variant
    corresponding to itself and 0 otherwise. For example, in the example above,
    0/2 maps to 0/0 and 0/1. The genotype 1/2 maps to 0/1 and 0/1.

    The biallelic alt `AD` entry is just the multiallelic `AD` entry
    corresponding to the alternate allele. The ref AD entry is the sum of the
    other multiallelic entries.

    The biallelic `DP` is the same as the multiallelic `DP`.

    The biallelic `PL` entry for a genotype g is the minimum over `PL` entries
    for multiallelic genotypes that downcode to g. For example, the `PL` for (A,
    T) at 0/1 is the minimum of the PLs for 0/1 (50) and 1/2 (45), and thus 45.

    Fixing an alternate allele and biallelic variant, downcoding gives a map
    from multiallelic to biallelic alleles and genotypes. The biallelic `AD` entry
    for an allele is just the sum of the multiallelic `AD` entries for alleles
    that map to that allele. Similarly, the biallelic `PL` entry for a genotype is
    the minimum over multiallelic `PL` entries for genotypes that map to that
    genotype.

    `GQ` is recomputed from `PL` if `PL` is provided and is not
    missing. If not, it is copied from the original GQ.

    Here is a second example for a het non-ref

    .. code-block:: text

      A   C,T 1/2:2,8,6:16:45:99,50,99,45,0,99

    splits as

    .. code-block:: text

      A   C   0/1:8,8:16:45:45,0,99
      A   T   0/1:10,6:16:50:50,0,99

    **VCF Info Fields**

    Hail does not split fields in the info field. This means that if a
    multiallelic site with `info.AC` value ``[10, 2]`` is split, each split
    site will contain the same array ``[10, 2]``. The provided allele index
    field `a_index` can be used to select the value corresponding to the split
    allele's position:

    >>> split_ds = hl.split_multi_hts(dataset)
    >>> split_ds = split_ds.filter_rows(split_ds.info.AC[split_ds.a_index - 1] < 10,
    ...                                 keep = False)

    VCFs split by Hail and exported to new VCFs may be
    incompatible with other tools, if action is not taken
    first. Since the "Number" of the arrays in split multiallelic
    sites no longer matches the structure on import ("A" for 1 per
    allele, for example), Hail will export these fields with
    number ".".

    If the desired output is one value per site, then it is
    possible to use annotate_variants_expr to remap these
    values. Here is an example:

    >>> split_ds = hl.split_multi_hts(dataset)
    >>> split_ds = split_ds.annotate_rows(info = split_ds.info.annotate(AC = split_ds.info.AC[split_ds.a_index - 1]))
    >>> hl.export_vcf(split_ds, 'output/export.vcf') # doctest: +SKIP

    The info field AC in *data/export.vcf* will have ``Number=1``.

    **New Fields**

    :func:`.split_multi_hts` adds the following fields:

     - `was_split` (*bool*) -- ``True`` if this variant was originally
       multiallelic, otherwise ``False``.

     - `a_index` (*int*) -- The original index of this alternate allele in the
       multiallelic representation (NB: 1 is the first alternate allele or the
       only alternate allele in a biallelic variant). For example, 1:100:A:T,C
       splits into two variants: 1:100:A:T with ``a_index = 1`` and 1:100:A:C
       with ``a_index = 2``.

    See Also
    --------
    :func:`.split_multi`

    Parameters
    ----------
    ds : :class:`.MatrixTable` or :class:`.Table`
        An unsplit dataset.
    keep_star : :obj:`bool`
        Do not filter out * alleles.
    left_aligned : :obj:`bool`
        If ``True``, variants are assumed to be left
        aligned and have unique loci. This avoids a shuffle. If the assumption
        is violated, an error is generated.
    vep_root : :class:`str`
        Top-level location of vep data. All variable-length VEP fields
        (intergenic_consequences, motif_feature_consequences,
        regulatory_feature_consequences, and transcript_consequences)
        will be split properly (i.e. a_index corresponding to the VEP allele_num).
    permit_shuffle : :obj:`bool`
        If ``True``, permit a data shuffle to sort out-of-order split results.
        This will only be required if input data has duplicate loci, one of
        which contains more than one alternate allele.
    recalculate_gq : :obj:`bool`
        If ``True``, recalcualte genotype quality scores based on PL

    Returns
    -------
    :class:`.MatrixTable` or :class:`.Table`
        A biallelic variant dataset.

    """

    split = hl.split_multi(ds, keep_star=keep_star, left_aligned=left_aligned, permit_shuffle=permit_shuffle)

    row_fields = set(ds.row)
    update_rows_expression = {}
    if vep_root in row_fields:
        update_rows_expression[vep_root] = split[vep_root].annotate(
            **{
                x: split[vep_root][x].filter(lambda csq: csq.allele_num == split.a_index)
                for x in (
                    "intergenic_consequences",
                    "motif_feature_consequences",
                    "regulatory_feature_consequences",
                    "transcript_consequences",
                )
            }
        )

    if isinstance(ds, Table):
        return split.annotate(**update_rows_expression).drop("old_locus", "old_alleles")

    split = split.annotate_rows(**update_rows_expression)
    entry_fields = ds.entry

    expected_field_types = {
        "GT": hl.tcall,
        "AD": hl.tarray(hl.tint),
        "DP": hl.tint,
        "GQ": hl.tint,
        "PL": hl.tarray(hl.tint),
        "PGT": hl.tcall,
        "PID": hl.tstr,
    }

    bad_fields = []
    for field in entry_fields:
        if field in expected_field_types and entry_fields[field].dtype != expected_field_types[field]:
            bad_fields.append((field, entry_fields[field].dtype, expected_field_types[field]))

    if bad_fields:
        msg = "\n  ".join([f"'{x[0]}'\tfound: {x[1]}\texpected: {x[2]}" for x in bad_fields])
        raise TypeError("'split_multi_hts': Found invalid types for the following fields:\n  " + msg)

    update_entries_expression = {}
    if "GT" in entry_fields:
        update_entries_expression["GT"] = hl.downcode(split.GT, split.a_index)
    if "DP" in entry_fields:
        update_entries_expression["DP"] = split.DP
    if "AD" in entry_fields:
        update_entries_expression["AD"] = hl.or_missing(
            hl.is_defined(split.AD), [hl.sum(split.AD) - split.AD[split.a_index], split.AD[split.a_index]]
        )
    if "PL" in entry_fields:
        pl = hl.or_missing(
            hl.is_defined(split.PL),
            (
                hl.range(0, 3).map(
                    lambda i: hl.min(
                        (
                            hl.range(0, hl.triangle(split.old_alleles.length()))
                            .filter(
                                lambda j: hl.downcode(
                                    hl.unphased_diploid_gt_index_call(j), split.a_index
                                ).unphased_diploid_gt_index()
                                == i
                            )
                            .map(lambda j: split.PL[j])
                        )
                    )
                )
            ),
        )
        update_entries_expression["PL"] = pl
        if "GQ" in entry_fields:
            if recalculate_gq:  # Recalculating GQ - the inherited Hail behavior
                update_entries_expression["GQ"] = hl.or_else(hl.gq_from_pl(pl), split.GQ)
            else:  # Keeping GQ the same - bcftools behavior
                update_entries_expression["GQ"] = split.GQ

    elif "GQ" in entry_fields:
        update_entries_expression["GQ"] = split.GQ

    if "PGT" in entry_fields:
        update_entries_expression["PGT"] = hl.downcode(split.PGT, split.a_index)
    if "PID" in entry_fields:
        update_entries_expression["PID"] = split.PID
    return split.annotate_entries(**update_entries_expression).drop("old_locus", "old_alleles")
