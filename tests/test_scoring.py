"""Tests for mapping scoring choices to blastp arguments."""

from frankensearch import blast, scoring


def test_identity_gapped_args():
    assert scoring.blast_args("identity", ungapped=False) == [
        "-matrix", "IDENTITY", "-gapopen", "15", "-gapextend", "2",
    ]


def test_identity_ungapped_args():
    assert scoring.blast_args("identity", ungapped=True) == ["-matrix", "IDENTITY", "-ungapped"]


def test_other_matrix_uses_blast_default_gaps():
    assert scoring.blast_args("blosum62", ungapped=False) == ["-matrix", "BLOSUM62"]


def test_matrix_name_is_uppercased():
    assert scoring.matrix_blast_name("pam30") == "PAM30"


def test_gap_description():
    assert scoring.gap_description("identity", ungapped=False) == "open 15, extend 2"
    assert scoring.gap_description("identity", ungapped=True) == "ungapped"
    # Non-identity matrices report blastp's actual default gap costs, tagged as such.
    defaults = {
        "pam30": "open 9, extend 1 (matrix default)",
        "blosum45": "open 14, extend 2 (matrix default)",
        "blosum62": "open 11, extend 1 (matrix default)",
    }
    for matrix, expected in defaults.items():
        assert scoring.gap_description(matrix, ungapped=False) == expected


def test_remote_args_identity_falls_back_to_pam30():
    args, warning = scoring.remote_blast_args("identity", ungapped=False)
    assert args == ["-matrix", "PAM30"]
    assert warning and "PAM30" in warning


def test_remote_args_other_matrix_no_warning():
    args, warning = scoring.remote_blast_args("blosum62", ungapped=False)
    assert args == ["-matrix", "BLOSUM62"]
    assert warning is None


def test_remote_args_ungapped():
    args, _ = scoring.remote_blast_args("blosum62", ungapped=True)
    assert args == ["-matrix", "BLOSUM62", "-ungapped"]


def test_selftest_identity_matrix_loads():
    # Runs the real blastp/makeblastdb available in the dev environment.
    ok, detail = blast.selftest_matrix("identity")
    assert ok, detail
