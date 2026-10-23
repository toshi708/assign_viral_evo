import numpy as np          

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
  differences = 0
  total = 0

  for i in range(len(reference_sequence)):
    reference = reference_sequence[i]
    distant = distant_sequence[i]

    if reference != "-" and distant != "-":
      total += 1

      if reference != distant:
        differences += 1

  if total == 0:
    return 0.0

  proportion = differences / total

  if proportion >= 0.75:
    return float("inf")

  return -0.75 * np.log(1 - (4 / 3) * proportion) * total


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()
